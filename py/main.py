import simulation
import output
import graph
import numpy as np
import matplotlib.pyplot as plt
import cv2

NX = 600#格子 x 轴 数目（格点数）
NY = 200#格子 y 轴 数目 （格点数）
Q = 9  # D2Q9 有九个方向
Rho0 = 1.0 # 密度
G = [0,0] # 压力梯度
G=np.array(G)
e = np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]])
#定义一个相反的方向量作为反弹边界的方向判定
re = np.array([0,3,4,1,2,7,8,5,6])#相反方向索引
ww = np.array([4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36,])# 各个方向的权重（）
Rho = np.ones((NX, NY,), dtype = np.float64)# 每个格点的总密度
U = np.ones((NX, NY, 2), dtype = np.float64) # 速度矢量
UU = np.ones((NX, NY,), dtype = np.float64) #速度标量值
#定义一个前一步的速度量来计算相邻两步之间的误差
U_pre = np.ones((NX, NY, 2), dtype = np.float64)
#密度分布函数
f = np.zeros((NX, NY, Q), dtype = np.float64)
#初始化
def init(U0,block,Rho0):
    for i in range(NX):
        for j in range(NY):
            if(block[i][j] == 1):#如果是墙壁格子，密度为0，速度为0
                Rho[i][j] = 0.0
                U[i][j] = [0.0, 0.0]
            else:
                Rho[i][j] = Rho0
                U[i][j] = 0
            #以平衡态分布函数作为初始的密度分布函数
            for m in range(Q):
                f[i][j][m] = simulation.feq(m, Rho[i][j], U[i][j][0], U[i][j][1])#按平衡态分布

#定义一个墙壁的格子，1代表墙壁，0代表流体
if __name__ == "__main__":
    output.clear_output()#清理之前的输出文件
    block = np.zeros((NX, NY), dtype = np.int32)
    #在中间放一个障碍物
    '''for i in range(80, 160):#这里填写X轴就好了
        for j in range(40, 160):#这里填写Y轴就好了
            block[i][j] = 1#这里交换了顺序'''
    #block=graph.generate_random_img(size=(NX,NY))
    block=graph.read_image(r"D:\coding\projects\prp\2.png")#读取图片（这里我是手画的）
    block=graph.clip_image(block,NX,NY)#这是转换成指定大小
    block=block/255#转换成掩码（（（

    Re = float(input("enter Re"))#输入雷诺数
    U1=float(input("enter Ux"))#输入初始速度x
    U2=float(input("enter Uy"))#输入初始速度y
    U0 = [U1,U2]
    U0=np.array(U0)
    viu = np.sqrt(U1*U1+U2*U2)*NY/Re#根据雷诺数计算粘度
    tau = 3*viu + 0.5#根据粘度
    print("Re: {}, U: {}, tau: {}".format(Re, U0, tau))

    init(U0, block, Rho0)
    f_col=np.zeros(shape=f.shape)
    max_iter = 3000
    imgs=[]
    for n in range(max_iter):
        U_pre = np.copy(U)
        f, Rho, U, UU = simulation.evolution(NX, NY, f,f_col, Rho, U, UU, tau, block ,Rho0, U0, G, Q ,e, ww, re)
        if(n % 5 == 0):#每5步输出一次结果
            print('Iter: {}, Error: {}'.format(n, simulation.error(U, U_pre)))
            output.writetecplot(NX, NY, Rho, U, n)
            #pic=output.draw_velocity_field(NX, NY, U, UU, Rho, n, save_fig=False)
            #imgs.append(pic)
        if(simulation.error(U, U_pre) < 2e-5):#误差足够小，认为收敛
            print('Converged at Iter: {}, Error: {}'.format(n, simulation.error(U, U_pre)))
            output.writetecplot(NX, NY, Rho, U, n)
            #pic=output.draw_velocity_field(NX, NY, U, UU, Rho,n, save_fig=False)
            #imgs.append(pic)
            break
    #将结果做成视频
    '''output_video = 'lbm_simulation.mp4'
    height, width, layers = imgs[0].shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(output_video, fourcc, 5, (width, height))
    for img in imgs:
        video.write(cv2.cvtColor(img, cv2.COLOR_RGBA2BGR))
    video.release()'''
    

    