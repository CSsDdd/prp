import math
import numpy as np
eD2Q9=[[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]#D2Q9模型的离散速度
reD2Q9=[0,3,4,1,2,7,8,5,6]#D2Q9模型的反向索引
wwD2Q9=[4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.]#D2Q9模型的权重
QD2Q9=9#D2Q9模型的离散速度数量
def feq(m, RHO, U, V,e=eD2Q9,ww=wwD2Q9,re=reD2Q9):#计算平衡态分布函数，参数表为 方向索引，密度，宏观速度（x轴），宏观速度（y轴）
    eu = e[m][0] * U + e[m][1] * V
    uv = U * U + V * V
    return ww[m] * RHO * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv)#公式可见视频

#NX,NY为网格的行列数，Q为离散速度方向的数量，f为密度分布函数，f_col为碰撞后的密度分布函数，Rho为密度，U为速度，UU为速度的模长，U_pre为上一步的速度，tau为弛豫时间，G为外力
def evolution(NX, NY, f, Rho, U, UU, tau, block ,G=0, Q=QD2Q9 ,e=eD2Q9, ww=wwD2Q9, re=reD2Q9):
    f_col = [[[0. for m in range(Q)] for j in range(NY)] for i in range(NX)]
    #求碰撞之后的密度分布函数
    for i in range(NX):
        for j in range(NY):
            for m in range(Q):
                f_col[i][j][m] = f[i][j][m] + (feq(m, Rho[i][j], U[i][j][0], U[i][j][1], e, ww, re) - f[i][j][m]) / tau + 3.0 * ww[m] * G * e[m][0];#计算碰撞
    #密度函数量的迁移加边界条件
    for i in range(NX):
        for j in range(NY):
            for m in range(Q):
                #i方向用余数的方法来获得周期性
                ip = (i - e[m][0] + NX) % NX
                jp = j - e[m][1]#计算j方向预期坐标
                if(jp < 0 or jp >= NY or block[i][j] == 1):#如果碰到边界或者墙壁
                    f[i][j][m] = f_col[i][j][re[m]]#如果碰到墙壁，碰撞后会反方向速度回到同一点
                else:
                    f[i][j][m] = f_col[ip][jp][m]#不碰撞就迁移
    #用迁移后的密度分布函数求宏观量
    for i in range(NX):
        for j in range(NY):
            Rho[i][j] = 0.#准备叠加
            U[i][j] = [0.,0.]#准备叠加
            for m in range(Q):
                Rho[i][j] += f[i][j][m]
                U[i][j][0] += f[i][j][m] * e[m][0]
                U[i][j][1] += f[i][j][m] * e[m][1]
            if(Rho[i][j] > 0):
                U[i][j][0] /= Rho[i][j]#取加权
                U[i][j][1] /= Rho[i][j]#取加权
                UU[i][j] = math.sqrt(U[i][j][0] * U[i][j][0] + U[i][j][1] * U[i][j][1])#算标量
           #else:
            #    U[i][j] = [0.,0.]
    return f, Rho, U, UU

def error(U, U_pre):
    U = np.array(U)
    U_pre = np.array(U_pre)
    error1 = np.sum((U - U_pre) ** 2)
    error2 = np.sum(U_pre ** 2)
    return error1 / (error2 + 1e-8)