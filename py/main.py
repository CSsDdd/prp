import simulation
import output
import graph
import numpy as np
import matplotlib.pyplot as plt
NX = 20#格子 x 轴 数目 （格点数）
NY = 20#格子 y 轴 数目（格点数）
Q = 9  # D2Q9 有九个方向
Rho0 = 1.0 # 密度
U0 = [0.,0.] # 速度向量 （vx,vy）
Vis = 0.6/3 #运动粘度 
G = 1e-4 # 压力梯度
tau = 3 * Vis + 0.5 # 这里根据运动粘度算出 松弛时间
e = np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]])
#定义一个相反的方向量作为反弹边界的方向判定
re = np.array([0,3,4,1,2,7,8,5,6])#相反方向索引
ww = np.array([4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36,])# 各个方向的权重（）
Rho = np.ones((NX,NY), dtype = np.float64)# 每个格点的总密度
U = np.ones((NX, NY, 2), dtype = np.float64) # 速度矢量
UU = np.ones((NX, NY), dtype = np.float64) #速度标量值
#定义一个前一步的速度量来计算相邻两步之间的误差
U_pre = np.ones((NX, NY, 2), dtype = np.float64)
#密度分布函数
f = np.zeros((NX, NY, Q), dtype = np.float64)
#初始化
def initial():
    for i in range(NX):
        for j in range(NY):
            Rho[i][j] = Rho0
            U[i][j] = U0
            #以平衡态分布函数作为初始的密度分布函数
            for m in range(Q):
                f[i][j][m] = simulation.feq(m, Rho[i][j], U[i][j][0], U[i][j][1])#按平衡态分布
#定义一个墙壁的格子，1代表墙壁，0代表流体
block = np.zeros((NX, NY), dtype = np.int32)
#在中间放一个障碍物
for i in range(8, 12):
    for j in range(8, 12):
        block[i][j] = 1
initial()
max_iter = 10000
for n in range(max_iter):
    U_pre = np.copy(U)
    f, Rho, U, UU = simulation.evolution(NX, NY, f, Rho, U, UU, tau, block ,G, Q ,e, ww, re)
    if(n % 1000 == 0):
        print('Iter: {}, Error: {}'.format(n, simulation.error(U, U_pre)))
        output.writetecplot(NX, NY, Rho, U, n)
        x = np.arange(0, NX, 1)
        y = np.arange(0, NY, 1)
        xx, yy = np.meshgrid(x, y)
        fig, ax = plt.subplots()
        ax.contourf(yy, xx, UU, color = "k")
        ax.quiver(yy, xx, U[:,:,0], U[:,:,1])
        fig.savefig('velocity_field_{:08d}.png'.format(n))
    if(simulation.error(U, U_pre) < 1e-6):
        print('Converged at Iter: {}, Error: {}'.format(n, simulation.error(U, U_pre)))
        output.writetecplot(NX, NY, Rho, U, n)
        break

    