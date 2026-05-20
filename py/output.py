import numpy as np
import matplotlib.pyplot as plt
import os 
def clear_output():
    for filename in os.listdir('.'):
        if filename.startswith('output_') and filename.endswith('.dat'):
            os.remove(filename)
        if filename.startswith('velocity_field_') and filename.endswith('.png'):
            os.remove(filename)

def writetecplot(NX, NY, rho, T, U, n):
    filename = 'output_{:08d}.dat'.format(n)
    with open(filename, 'w') as f:
        f.write('TITLE = "LBM Simulation"\n')
        f.write('VARIABLES = "X", "Y", "rho","T" ,"u", "v"\n')
        f.write('ZONE T="Zone {}", I={}, J={}, F=POINT\n'.format(n//1000, NX, NY))
        for j in range(NY):
            for i in range(NX):
                if(not rho[i][j]==rho[i][j]):
                    print("nan in rho")
                if(not T[i][j]==T[i][j]):
                    print("nan in T")
                if(not(U[i][j][0]==U[i][j][0]) or not(U[i][j][1]==U[i][j][1])):
                    print("nan in U")
                f.write('{} {} {} {} {} {}\n'.format(i, j, rho[i][j], T[i][j], U[i][j][0], U[i][j][1]))

def draw_velocity_field(NY, NX, U, UU, Rho, n, save_fig=False):
    x = np.arange(0, NY, 1)
    y = np.arange(0, NX, 1)
    xx, yy = np.meshgrid(y, x)#符合矩阵行为,生成NY*NX的坐标矩阵
    fig, ax = plt.subplots()
    ax.contourf(xx, yy, Rho, color = "k")
    ax.quiver(xx, yy, U[:,:,0], U[:,:,1])
    if(save_fig):
        fig.savefig('velocity_field_{:08d}.png'.format(n))
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())
    return arr

def draw_Rho_field(NY, NX, Rho, n, save_fig=False):
    x = np.arange(0, NY, 1)
    y = np.arange(0, NX, 1)
    xx, yy = np.meshgrid(y, x)#符合矩阵行为,生成NY*NX的坐标矩阵
    fig, ax = plt.subplots()
    ax.contourf(xx, yy, Rho, color = "k")
    if(save_fig):
        fig.savefig('velocity_field_{:08d}.png'.format(n))
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())
    return arr