import math
import numpy as np
import numba
import pathlib
import tomllib
from enum import IntEnum 
import output
import graph
eD2Q9=[[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]#D2Q9模型的离散速度
eD2Q9=np.array(eD2Q9)
reD2Q9=[0,3,4,1,2,7,8,5,6]#D2Q9模型的反向索引
reD2Q9=np.array(reD2Q9)
wwD2Q9=[4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.]#D2Q9模型的权重
wwD2Q9=np.array(wwD2Q9)
QD2Q9=9#D2Q9模型的离散速度数量

class controller:
    def __init__(self,config_path,image_path,
                 D=2,Q=9,G=np.zeros(shape=(2))):
        # 读取物理量和格子配置
        config_path = pathlib.Path(config_path)
        with config_path.open("rb") as f:
            config = tomllib.load(f)
        #两类，格子参数
        lattice = config["Lattice"]#格子参数
        physics = config["Physics"]#物理参数

        #格子基本信息：
        self.Nx, self.Ny = lattice["Nx"], lattice["Ny"]
        self.U0 = np.array([lattice["U0x"], lattice["U0y"]])
        
        #物理基本信息
        self.Up = np.array([physics["Ux"], physics["Uy"]])#物理流速
        self.Th = physics["Th"]
        self.Tc = physics["Tc"]
        self.P_total = physics["P"]
        # 物性相关
        # 计算固体和流体的热扩散率
        self.Rho_s=physics["Rho_s"]
        self.Rho_l=physics["Rho_l"]
        self.Chi_s = physics["kappa_s"]/(physics["Rho_s"]*physics["C_s"])
        self.Chi_l = physics["kappa_l"]/(physics["Rho_l"]*physics["C_l"])
        self.Nu_l = physics["Nu_l"]#参考Tc时数值，运动粘度
        self.beta = physics["beta"]


        # 读取图像并构造三段式通道掩码
        self.block = graph.read_image(str(image_path))
        self.block = graph.clip_image(self.block, self.Nx, self.Ny)
        self.block = np.round(self.block)
        empty = np.zeros((round(self.Nx*0.25), self.Ny), dtype=self.block.dtype)
        self.block = np.concatenate([empty, self.block, empty], axis=0)/255
        self.Chi = np.array(self.block, dtype=np.float64, copy=True)*(self.Chi_s-self.Chi_l) + self.Chi_l
        self.source_delta_T=np.array(self.block, dtype=np.float64, copy=True)
        self.source_delta_T=self.source_delta_T/np.sum(self.source_delta_T)*self.P_total/(physics["Rho_s"]*physics["C_s"])

        # 三段式通道对应的格子数和物理长度
        self.Nx, self.Ny = self.block.shape[0], self.Ny
        #print(self.Nx,self.Ny)
        self.X = physics["X"]*1.5
        
        self.delta_T = self.Th-self.Tc

        #构建物理量到格子量过渡工具
        self.deltax=self.X/self.Nx
        self.deltat=self.deltax*(self.U0[0]/self.Up[0])

        #格子量
        self.re = np.array([0,3,4,1,2,7,8,5,6])#相反方向索引
        self.ww = np.array([4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36,])# 各个方向的权重（）

        self.block = np.round(self.block)
        

        self.Rho = np.ones((self.Nx, self.Ny,), dtype = np.float64)# 每个格点的总密度
        self.U = np.ones((self.Nx, self.Ny, D), dtype = np.float64) # 速度矢量
        self.UU = np.ones((self.Nx, self.Ny,), dtype = np.float64) #速度标量值
        self.U_pre = np.ones((self.Nx, self.Ny, D), dtype = np.float64)#定义一个前一步的速度量来计算相邻两步之间的误差
        self.f = np.zeros((self.Nx, self.Ny, Q), dtype = np.float64)#速度分布函数
        self.f_col=np.zeros(shape=self.f.shape)#碰撞使用的临时速度分布变量

        self.T = np.ones((self.Nx,self.Ny),dtype=np.float64)#每个格点的总热量（归一化后）
        self.g = np.zeros((self.Nx, self.Ny, Q), dtype = np.float64)#温度分布函数
        self.g_col=np.zeros(shape=self.g.shape)#碰撞使用的临时温度分布变量

        self.col_stat=np.zeros(shape=self.g.shape)

        self.Rho0=1
        self.T0=0.0#取入口水温为低温
        self.G=G

        #物性量纲调整

        self.Chi = self.Chi * (self.deltat / (self.deltax ** 2))#修正
        self.source_delta_T = self.source_delta_T*self.deltat/(self.deltax**3 * self.delta_T)
        self.Nu_l=self.Nu_l / (self.deltax**2/self.deltat)#算出格子运动粘度
        self.beta = self.beta*self.delta_T
        self.Chi_l = self.Chi_l * (self.deltat / (self.deltax ** 2))
        self.Chi_s = self.Chi_s * (self.deltat / (self.deltax ** 2))
        #
        self.U_tau=3*self.Nu_l+0.5#算出U_tau暂不考虑温度引起的运动粘度变化，待修正

        # 可靠性判据的历史状态
        self.U_pre=None
        self.T_pre=None
        self._reliability_stable_checks = 0

        #打印指标，用于参考
        self.Re=np.sqrt(np.sum(self.U0**2))*self.Nx/self.Nu_l
        self.Pr=self.Nu_l/self.Chi_l
        self.Ra=((9.8*(self.deltat**2/self.deltax)
                 *self.beta
                 *1
                 *(self.Nx**3))/(
                 self.Chi_l
                 *self.Nu_l
                ))
        print("Re:",self.Re)
        print("Pr:",self.Pr)
        print("Ra:",self.Ra)

        for i in range(self.Nx):
            for j in range(self.Ny):
                if(self.block[i][j] == 1):#如果是墙壁格子，密度为0，速度为0
                    self.Rho[i][j] = 1.0#调整尺度
                    self.U[i][j] = np.zeros(shape=(2))#无速度
                    self.T[i][j]=0.0
                else:
                    self.Rho[i][j] = self.Rho0
                    self.U[i][j] = self.U0
                    self.T[i][j]= 0.0
                #以平衡态分布函数作为初始的密度分布函数
                for m in range(Q):
                    self.f[i][j][m] =f_eq(m, self.Rho[i][j], self.U[i][j])#按平衡态分布
                    self.g[i][j][m] =g_eq(m, self.T[i][j], self.U[i][j])

    def run_one_step(self):
        self.U_pre=self.U.copy()
        self.T_pre=self.T.copy()
        self.f, self.Rho, self.U, self.UU, self.T = evolution(NX=self.Nx, NY=self.Ny, 
                                                f=self.f,f_col=self.f_col,Rho=self.Rho,
                                                U=self.U,UU=self.UU,U_tau=self.U_tau,
                                                g=self.g, g_col=self.g_col, T=self.T,
                                                block= self.block, Chi=self.Chi, source_delta_T=self.source_delta_T,
                                                Rho0=self.Rho0,U0=self.U0,T0=self.T0,G=self.G,
                                                col_stat=self.col_stat)

    def report_reliability(self, velocity_tol=1e-6, temperature_tol=1e-6,peak_tol=1e-6,
                           required_checks=3):
        # 采集当前速度场、温度场和出口带热量
        source_mask = self.source_delta_T > 0.0
        if not np.any(source_mask):
            source_peak = 0
        else:
            source_peak = np.max(self.T[source_mask])
        
        # 首次调用只建立比较基准
        first_check = self.T_pre is None
        if first_check:
            velocity_error = np.inf
            velocity_error_per = np.inf

            temperature_error = np.inf
            temperature_error_per = np.inf
            
            peak_error = np.inf
            peak_error_per = np.inf

        if not first_check:
            #1.计算流场收敛程度
            velocity_error = np.sqrt(np.mean((self.U_pre-self.U)**2))
            velocity_error_per=velocity_error/max(1e-12,np.abs(np.sqrt(np.mean(self.U**2))),np.abs(np.sqrt(np.mean(self.U_pre**2))))
            #2.计算温度场收敛程度
            temperature_error = np.sqrt(np.mean((self.T_pre-self.T)**2))
            temperature_error_per=temperature_error/max(1e-12,np.abs(np.sqrt(np.mean(self.T**2))),np.abs(np.sqrt(np.mean(self.T_pre**2))))
            #3.最高温度收敛程度
            if not np.any(source_mask):
                peak_error = 0
                peak_error_per = 0
            else:
                peak_error = abs(source_peak-np.max(self.T_pre[source_mask]))
                peak_error_per = peak_error/max(1e-12,np.max(self.T[source_mask]),np.max(self.T_pre[source_mask]))

        # 要求速度和三个热稳态指标连续多次同时达标
        stable = (velocity_error_per < velocity_tol and temperature_error_per < temperature_tol and
                  peak_error_per < peak_tol )
        
        self._reliability_stable_checks = (
            self._reliability_stable_checks+1 if stable else 0
        )
        converged = self._reliability_stable_checks >= required_checks
        return {
            "converged": bool(converged),
            "stable_checks": int(self._reliability_stable_checks),
            "velocity_field_relative_change": float(velocity_error_per),
            "temperature_field_relative_change": float(temperature_error_per),
            "source_peak_temperature_relative_change": float(peak_error_per),
        }

    def report_performance(self):
        #统计最高温
        source_mask = self.source_delta_T > 0.0
        if not np.any(source_mask):
            heat_Peak_l=0
        else: 
            heat_Peak_l=np.max(self.T[source_mask])
        heat_Peak_h=self.Tc+heat_Peak_l*(self.Th-self.Tc)

        #统计入口处压差
        avg_Income_Rho=np.mean(self.Rho[0][1:-1:1])
        P_l=1/3*(avg_Income_Rho-1)#末尾参考为1
        P_p=P_l*self.Rho_l*(self.deltax**2)/(self.deltat**2)

        #统计入口体积流量
        V_l=np.sum(self.U[0, 1:-1, 0])
        V_p=V_l*(self.deltax**3)/(self.deltat)

        return {
            "max_source_temperature_k": float(heat_Peak_h),
            "pressure_drop_pa": float(P_p),
            "volume_flow_rate_m3_s": float(V_p),
        }

    def final_report(self, reliability):
        report = {
            "performance": self.report_performance(),
            "reliability": self.report_reliability(),
        }
        #print(report)
        return report

    def write(self,n):
        output.writetecplot(self.Nx, self.Ny, self.Rho, self.T, self.U, n)

class col_status(IntEnum):
    Finished=0
    Open_Boundary=1
    Bouncy_Boundary=2
    Block_Liquid=3
    Block_Solid=4

@numba.njit
def f_eq(m, RHO, U,e=eD2Q9,ww=wwD2Q9,re=reD2Q9):#计算平衡态分布函数，参数表为 方向索引，密度，宏观速度（x轴），宏观速度（y轴）
    eu = e[m][0] * U[0] + e[m][1] * U[1]
    uv = U[0] * U[0] + U[1] * U[1]
    return ww[m] * RHO * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv)#公式可见视频

@numba.njit
def g_eq(m,T,U,e=eD2Q9,ww=wwD2Q9,re=reD2Q9):#温度分布函数的平衡态计算公式，但是返回的是温度分布函数的平衡态值（温度）
	eu=e[m][0]*U[0]+e[m][1]*U[1]
	x=ww[m]*T*(1.0+3*eu)
	return x

#NX,NY为网格的行列数，Q为离散速度方向的数量，f为密度分布函数，f_col为碰撞后的密度分布函数，Rho为密度，U为速度，UU为速度的模长，U_pre为上一步的速度，tau为弛豫时间，G为外力,Cs为固体比热容，Cl为流体比热容（均为单位质量）
@numba.njit
def evolution(NX, NY,
              Rho,
              f, f_col, 
              U, UU, U_tau,
              g, g_col, 
              T,
              block, source_delta_T, Chi,
              col_stat,
              Rho0, U0=[0,0], T0=0.5, Cs=1, Cl=1, G=[0,0], Q=QD2Q9 ,e=eD2Q9, ww=wwD2Q9, re=reD2Q9):
    #step1:求碰撞之后的密度分布函数
    for i in range(NX):
        for j in range(NY):
            alpha=Chi[i][j]*1/Rho[i][j]
            T_tau=0.5+3*alpha#实际上τ=0.5+3*α*Δt/（Δx）
            if(block[i][j] == 1):#如果是墙壁格子，不进行操作
                
                #正常更新
                for m in range(Q):
                    g_col[i][j][m] = (g[i][j][m] +
                                    (g_eq(m,T[i][j],[0,0])-g[i][j][m])/T_tau
                                    )+ww[m]*source_delta_T[i][j]*Rho[i][j]#虽然这里Rho[i][j]是多余的，但是处于规范还是加了
                    #计算热量碰撞
                    f_col[i][j][m] = f[i][j][m]#密度不动
            else:
                for m in range(Q):
                    g_col[i][j][m] = (g[i][j][m] +
                                    (g_eq(m,T[i][j],U[i][j])-g[i][j][m])/T_tau
                                    )#计算热量碰撞
                    
                    f_col[i][j][m] = (f[i][j][m] + 
                                    (f_eq(m, Rho[i][j], U[i][j]) - f[i][j][m]) / U_tau + 
                                    ww[m]*(1-1/(2*U_tau))*3*(e[m][0]*G[0]+e[m][1]*G[1])#注：浮力项不需要考虑，因为我们考虑的是水平截面，竖直方向的不用统计。这里的G是梯度力……
                                    )#计算速度碰撞
                    
                    
    #step2各个量第一轮迁移。同时完成迁移情况分类
    for i in range(NX):
        for j in range(NY):
            col_stat[i][j]=col_status.Finished.value#默认改为0
            is_open_section = j > 0 and j < NY-1

            # 1. 入口格点：给定速度和温度，全方向重构
            if(i == 0 and is_open_section):
                col_stat[i][j]=col_status.Open_Boundary.value
                ref_i = 1
                ref_Rho = Rho[ref_i][j]
                ref_U = U[ref_i][j]
                ref_T = T[ref_i][j]
                if(block[i][j] == 1):#固体密度不动，只重构温度分布
                    for m in range(Q):
                        f[i][j][m]=f_col[i][j][m]
                        g_neq=g_col[ref_i][j][m]-g_eq(m,ref_T,ref_U)
                        g[i][j][m]=g_eq(m,T0,np.array([0.0,0.0]))+g_neq
                else:
                    for m in range(Q):
                        f_neq=f_col[ref_i][j][m]-f_eq(m,ref_Rho,ref_U)
                        f[i][j][m]=f_eq(m,ref_Rho,U0)+f_neq
                        g_neq=g_col[ref_i][j][m]-g_eq(m,ref_T,ref_U)
                        g[i][j][m]=g_eq(m,T0,U0)+g_neq
                continue

            # 2. 出口格点：给定密度，速度和温度沿用内部参考值
            if(i == NX-1 and is_open_section):
                col_stat[i][j]=col_status.Open_Boundary.value
                ref_i = NX-2
                ref_Rho = Rho[ref_i][j]
                ref_U = U[ref_i][j]
                ref_T = T[ref_i][j]
                if(block[i][j] == 1):#固体密度不动，只重构温度分布
                    for m in range(Q):
                        f[i][j][m]=f_col[i][j][m]
                        g_neq=g_col[ref_i][j][m]-g_eq(m,ref_T,ref_U)
                        g[i][j][m]=g_eq(m,ref_T,np.array([0.0,0.0]))+g_neq
                else:
                    for m in range(Q):
                        f_neq=f_col[ref_i][j][m]-f_eq(m,ref_Rho,ref_U)
                        f[i][j][m]=f_eq(m,Rho0,ref_U)+f_neq
                        g_neq=g_col[ref_i][j][m]-g_eq(m,ref_T,ref_U)
                        g[i][j][m]=g_eq(m,ref_T,ref_U)+g_neq
                continue

            current_is_solid = block[i][j] == 1
            for m in range(Q):
                # pull-streaming 来源坐标
                src_i = i - e[m][0]
                src_j = j - e[m][1]
                source_is_outside = (src_i < 0 or src_i >= NX or
                                     src_j < 0 or src_j >= NY)

                # 3. 该方向碰到上下绝热边界，或位于四个角点
                if(source_is_outside):
                    col_stat[i][j]=col_status.Bouncy_Boundary.value
                    if(current_is_solid):
                        f[i][j][m]=f_col[i][j][m]
                    else:
                        f[i][j][m]=f_col[i][j][re[m]]
                    g[i][j][m]=g_col[i][j][re[m]]
                    continue

                # 4a. 当前格点是固体：密度不动，温度分布正常传递
                if(current_is_solid):
                    if(block[src_i][src_j] == 0):
                        col_stat[i][j]=col_status.Block_Solid.value
                    f[i][j][m]=f_col[i][j][m]
                    g[i][j][m]=g_col[src_i][src_j][m]
                    continue

                # 4b. 流体从固体方向拉取：动量反弹，热量跨界传递
                if(block[src_i][src_j] == 1):
                    col_stat[i][j]=col_status.Block_Liquid.value
                    f[i][j][m]=f_col[i][j][re[m]]
                    g[i][j][m]=g_col[src_i][src_j][m]
                    continue

                # 5. 正常迁移
                f[i][j][m]=f_col[src_i][src_j][m]
                g[i][j][m]=g_col[src_i][src_j][m]
        
    #用迁移后的密度分布函数求宏观量
    for i in range(NX):
        for j in range(NY):
            Rho[i][j] = 0.#准备叠加
            U[i][j][0] =0.#准备叠加
            U[i][j][1] =0.
            T[i][j]=0.
            for m in range(Q):
                Rho[i][j] += f[i][j][m]
                U[i][j][0] += f[i][j][m] * e[m][0]
                U[i][j][1] += f[i][j][m] * e[m][1]
                T[i][j]+=g[i][j][m]
            if(Rho[i][j] > 0):
                U[i][j][0] /= Rho[i][j]#取加权
                U[i][j][1] /= Rho[i][j]#取加权
                UU[i][j] = math.sqrt(U[i][j][0] * U[i][j][0] + U[i][j][1] * U[i][j][1])#算标量
            if block[i][j] == 1:
                U[i][j][0] = 0.0
                U[i][j][1] = 0.0
                UU[i][j] = 0.0
    return f, Rho, U, UU, T

def error(U, U_pre):#未使用
    U = np.array(U)
    U_pre = np.array(U_pre)
    error1 = np.sum((U - U_pre) ** 2)
    error2 = np.sum(U_pre ** 2)
    return error1 / (error2 + 1e-8)
