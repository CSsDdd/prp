import simulation
import output
import graph
import numpy as np
import matplotlib.pyplot as plt
import cv2
import tomllib
import pathlib

#定义一个墙壁的格子，1代表墙壁，0代表流体
if __name__ == "__main__":
    output.clear_output()#清理之前的输出文件
    
    configpth=pathlib.Path("config.toml")
    if configpth.exists():
        with open("config.toml", "rb") as f:
            config = tomllib.load(f)
    else:
        print("请给出物理量配置！")
        exit(0)
    Nx = config["Lattice"]["Nx"]
    Ny = config["Lattice"]["Ny"]
    U0x= config["Lattice"]["U0x"]#输入初始速度x
    U0y= config["Lattice"]["U0y"]#输入初始速度y
    U0 = [U0x,U0y]
    U0=np.array(U0)

    X = config["Physics"]["X"]#输入真实物理长度（m）
    Re = config["Physics"]["Re"]#输入雷诺数
    Pr= config["Physics"]["Pr"]#输入普朗克数
    beta = config["Physics"]["beta"]#输入热膨胀系数
    Ra = config["Physics"]["Ra"]#输入瑞利数
    Ux= config["Physics"]["Ux"]#输入初始速度x
    Uy= config["Physics"]["Uy"]#输入初始速度y
    U = [Ux,Uy]
    U=np.array(U)
    kappa_s = config["Physics"]["kappa_s"]#输入固体热传导系数
    C_s = config["Physics"]["C_s"]#输入固体比热容
    Rho_s = config["Physics"]["Rho_s"]#输入固体密度
    Chi_s=kappa_s/(Rho_s*C_s)#默认归一化密度为1时，使用时必须处理密度修正！

    kappa_l = config["Physics"]["kappa_l"]#输入流体热传导系数
    C_l = config["Physics"]["C_l"]#输入流体比热容
    Rho_l = config["Physics"]["Rho_l"]#输入流体密度
    Chi_l=kappa_l/(Rho_l*C_l)#默认归一化密度为1时，使用时必须处理密度修正！

    block = np.zeros((Nx, Ny), dtype = np.int32)
    #在中间放一个障碍物
    '''for i in range(80, 160):#这里填写X轴就好了
        for j in range(40, 160):#这里填写Y轴就好了
            block[i][j] = 1#这里交换了顺序'''
    #block=graph.generate_random_img(size=(NX,NY))
    block=graph.read_image(r"D:\coding\projects\prp\cell_random_img2.png")#读取图片
    block=graph.clip_image(block,Nx,Ny)#这是转换成指定大小
    block=np.concatenate([np.zeros((Nx,Ny),dtype = block.dtype),block,np.zeros((Nx,Ny),dtype = block.dtype)],axis=0)
    print(block.shape)
    block=block/255#转换成掩码（（（

    Chi=block*(Chi_s-Chi_l)
    Chi+=Chi_l

    controller=simulation.controller(NX=Nx*3,NY=Ny,U0=U0,
                                     block=block,Chi=Chi,
                                     X=X*3,U=U,Re=Re,Pr=Pr,beta=beta,Ra=Ra,
                                     D=2,Q=9)
    
    max_iter = 10000
    for n in range(max_iter):
        controller.run_one_step()        
        if(n % 5 == 0):#每5步输出一次结果
            err=controller.report_err()
            print('Iter: {}, Error: {}'.format(n, err))
            controller.write(n)#只传个编号就好了
            if(err < 1e-6):#误差足够小，认为收敛
                print('Converged at Iter: {}, Error: {}'.format(n, err))
                controller.write(n)
                break
    #将结果做成视频

    