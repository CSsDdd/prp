import math
import numpy as np
import numba
from enum import IntEnum 
import output
eD2Q9=[[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]#D2Q9模型的离散速度
eD2Q9=np.array(eD2Q9)
reD2Q9=[0,3,4,1,2,7,8,5,6]#D2Q9模型的反向索引
reD2Q9=np.array(reD2Q9)
wwD2Q9=[4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.]#D2Q9模型的权重
wwD2Q9=np.array(wwD2Q9)
QD2Q9=9#D2Q9模型的离散速度数量

class controller:
    def __init__(self,NX,NY,U0,
                 block,Chi,
                 X,U,Re,Pr,beta,Ra,
                 D=2,Q=9,G=np.zeros(shape=(2))):
        self.block=block
        self.Chi=Chi

        self.re = np.array([0,3,4,1,2,7,8,5,6])#相反方向索引
        self.ww = np.array([4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36,])# 各个方向的权重（）

        self.NX=NX
        self.NY=NY

        self.deltax=X/NX
        self.deltat=self.deltax*(U0[0]/U[0])

        self.Chi = self.Chi * (self.deltat / (self.deltax ** 2))#修正

        viu = np.sqrt(np.sum(U0*U0))*NY/Re#根据雷诺数计算粘度
        self.U_tau = 3*viu + 0.5#根据运动粘度算出速度松弛系数
        self.T_tau = 3*viu/Pr + 0.5 #根据运动粘度和普朗克数算出温度松弛系数
        self.delta_T = Ra*(viu**2)/((NY**3)*9.8*beta*Pr)#根据普朗克数算出deltaT（开尔文，关联到归一化后的T的1.0）

        self.Rho = np.ones((NX, NY,), dtype = np.float64)# 每个格点的总密度
        self.U = np.ones((NX, NY, D), dtype = np.float64) # 速度矢量
        self.UU = np.ones((NX, NY,), dtype = np.float64) #速度标量值
        self.U_pre = np.ones((NX, NY, D), dtype = np.float64)#定义一个前一步的速度量来计算相邻两步之间的误差
        self.f = np.zeros((NX, NY, Q), dtype = np.float64)#速度分布函数
        self.f_col=np.zeros(shape=self.f.shape)#碰撞使用的临时速度分布变量

        self.T = np.ones((NX,NY),dtype=np.float64)#每个格点的总热量（归一化后）
        self.g = np.zeros((NX, NY, Q), dtype = np.float64)#温度分布函数
        self.g_col=np.zeros(shape=self.g.shape)#碰撞使用的临时温度分布变量

        self.col_stat=np.zeros(shape=self.g.shape)

        self.U0=U0
        self.Rho0=1
        self.Th=1#热的
        self.T0=0.5#默认水温
        self.Tc=0#冷的
        self.G=G

        for i in range(NX):
            for j in range(NY):
                if(block[i][j] == 1):#如果是墙壁格子，密度为0，速度为0
                    self.Rho[i][j] = 1.0#调整尺度
                    self.U[i][j] = np.zeros(shape=(2))#无速度
                    self.T[i][j]=1
                else:
                    self.Rho[i][j] = self.Rho0
                    self.U[i][j] = self.U0
                    self.T[i][j]= self.T0
                #以平衡态分布函数作为初始的密度分布函数
                for m in range(Q):
                    self.f[i][j][m] =f_eq(m, self.Rho[i][j], self.U[i][j])#按平衡态分布
                    self.g[i][j][m] =g_eq(m, self.T[i][j], self.U[i][j])
            #self.U[i][0] = np.zeros(shape=(2))#无速度
            #self.T[i][0]=self.Tc
            #self.U[i][NY-1] = np.zeros(shape=(2))#无速度
            #self.T[i][NY-1]=self.Tc

    def run_one_step(self):
        self.f, self.Rho, self.U, self.UU, self.T = evolution(NX=self.NX, NY=self.NY, 
                                                f=self.f,f_col=self.f_col,Rho=self.Rho,
                                                U=self.U,UU=self.UU,U_tau=self.U_tau,
                                                g=self.g, g_col=self.g_col, T=self.T,
                                                block= self.block, Chi=self.Chi,
                                                Rho0=self.Rho0,U0=self.U0,T0=self.T0,G=self.G,
                                                col_stat=self.col_stat)
        
    def report_err(self):
        return error(self.U,self.U_pre)

    def write(self,n):
        output.writetecplot(self.NX, self.NY, self.Rho, self.T, self.U, n)

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
              block, Chi,
              col_stat, 
              Rho0, U0=[0,0], T0=0.5, Cs=1, Cl=1, G=[0,0], Q=QD2Q9 ,e=eD2Q9, ww=wwD2Q9, re=reD2Q9):
    #step1:求碰撞之后的密度分布函数
    for i in range(NX):
        for j in range(NY):
            alpha=Chi[i][j]*1/Rho[i][j]
            T_tau=0.5+3*alpha#实际上τ=0.5+3*α*Δt/（Δx）²

            if(block[i][j] == 1):#如果是墙壁格子，不进行操作
                for m in range(Q):
                    g_col[i][j][m] = (g[i][j][m] +
                                    (g_eq(m,T[i][j],[0,0])-g[i][j][m])/T_tau
                                    )#计算热量碰撞
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
            if(block[i][j] == 1):#如果是墙壁格子，不进行操作
                for m in range(Q):
                    ip = i - e[m][0]#计算i方向预期坐标
                    jp = j - e[m][1]#计算j方向预期坐标
                    if(ip<0 or ip>=NX):
                        col_stat[i][j]=col_status.Open_Boundary.value
                        f[i][j][m]=f_col[i][j][m]

                        g_neq=g_col[i+e[m][0]][j+e[m][1]][m]-g_eq(m,T[i+e[m][0]][j+e[m][1]],U[i+e[m][0]][j+e[m][1]])#同样存在不连续风险。
                        g[i][j][m]=g_eq(m,T0,U0)+g_neq#这里外面的速度U修正，用里面的估计，稍微提高一点连续性
                        continue
                    elif(jp < 0 or jp >= NY):
                        col_stat[i][j]=col_status.Bouncy_Boundary.value
                        f[i][j][m]=f_col[i][j][m]
                        g[i][j][m]=g_col[i][j][re[m]]
                        continue
                    elif(block[ip][jp] == 0):#如果是来自流体的
                        col_stat[i][j]=col_status.Block_Solid.value
                        f[i][j][m] = f_col[i][j][m]
                        g[i][j][m] = g_col[ip][jp][m]
                        continue
                    else:
                        f[i][j][m] = f_col[i][j][m]#固体密度不动（默认没速度）
                        g[i][j][m] = g_col[ip][jp][m]
            else:
                for m in range(Q):
                    ip = i - e[m][0]#计算i方向预期坐标
                    jp = j - e[m][1]#计算j方向预期坐标
                    if(ip<0 or ip>=NX):
                        col_stat[i][j]=col_status.Open_Boundary.value
                        f_neq=f_col[i+e[m][0]][j+e[m][1]][m]-f_eq(m,Rho[i+e[m][0]][j+e[m][1]],U[i+e[m][0]][j+e[m][1]])#密度会不匹配吗？
                        f[i][j][m]=f_eq(m,Rho0,U0)+f_neq#这里外面的速度U修正，用里面的估计，稍微提高一点连续性

                        g_neq=g_col[i+e[m][0]][j+e[m][1]][m]-g_eq(m,T[i+e[m][0]][j+e[m][1]],U[i+e[m][0]][j+e[m][1]])#同样存在不连续风险。
                        g[i][j][m]=g_eq(m,T0,U0)+g_neq#这里外面的速度U修正，用里面的估计，稍微提高一点连续性
                        continue#写switchcase习惯了
                    elif(jp < 0 or jp >= NY ):
                        col_stat[i][j]=col_status.Bouncy_Boundary.value
                        f[i][j][m]=f_col[i][j][re[m]]
                        g[i][j][m]=g_col[i][j][re[m]]
                        continue
                    elif(block[ip][jp] == 1):#如果碰到边界或者墙壁
                        col_stat[i][j]=col_status.Block_Liquid.value
                        f[i][j][m] = f_col[i][j][re[m]]
                        g[i][j][m] = g_col[ip][jp][m]
                        continue
                    else:
                        f[i][j][m] = f_col[ip][jp][m]#不碰撞就迁移
                        g[i][j][m] = g_col[ip][jp][m]
                        continue
        
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

def error(U, U_pre):
    U = np.array(U)
    U_pre = np.array(U_pre)
    error1 = np.sum((U - U_pre) ** 2)
    error2 = np.sum(U_pre ** 2)
    return error1 / (error2 + 1e-8)