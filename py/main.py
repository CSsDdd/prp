import simulation
import output
import graph
import numpy as np
import matplotlib.pyplot as plt
import cv2

NX = 200*3#格子 x 轴 数目（格点数）#记得乘以3保证是3的整数倍！
NY = 200#格子 y 轴 数目 （格点数）
#定义一个墙壁的格子，1代表墙壁，0代表流体
if __name__ == "__main__":
    output.clear_output()#清理之前的输出文件
    block = np.zeros((NX//3, NY), dtype = np.int32)
    #在中间放一个障碍物
    '''for i in range(80, 160):#这里填写X轴就好了
        for j in range(40, 160):#这里填写Y轴就好了
            block[i][j] = 1#这里交换了顺序'''
    #block=graph.generate_random_img(size=(NX,NY))
    block=graph.read_image(r"D:\coding\projects\prp\cell_random_img.jpg")#读取图片（这里我是手画的）
    block=graph.clip_image(block,NX//3,NY)#这是转换成指定大小
    block=np.concatenate([np.zeros((NX//3,NY),dtype = block.dtype),block,np.zeros((NX//3,NY),dtype = block.dtype)],axis=0)
    print(block.shape)
    block=block/255#转换成掩码（（（
    
    Re = float(input("enter Re"))#输入雷诺数
    Pr= float(input("enter Pr"))#输入普朗克数
    beta = float(input("enter β"))#输入热膨胀系数
    Ra = float(input("enter Ra"))#输入瑞利数
    U1=float(input("enter Ux"))#输入初始速度x
    U2=float(input("enter Uy"))#输入初始速度y
    U0 = [U1,U2]
    U0=np.array(U0)
    controller=simulation.controller(NX=NX,NY=NY,block=block,
                                     Re=Re,Pr=Pr,beta=beta,Ra=Ra,U0=U0,
                                     g=9.8,D=2,Q=9)
    
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

    