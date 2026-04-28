import numpy as np
import matplotlib.pyplot as plt
def writetecplot(NX, NY, rho, U, n):
    filename = 'output_{:08d}.dat'.format(n)
    with open(filename, 'w') as f:
        f.write('TITLE = "LBM Simulation"\n')
        f.write('VARIABLES = "X", "Y", "rho", "u", "v"\n')
        f.write('ZONE T="Zone {}", I={}, J={}, F=POINT\n'.format(n//1000, NX, NY))
        for j in range(NY):
            for i in range(NX):
                f.write('{} {} {} {} {}\n'.format(i, j, rho[j][i], U[j][i][0], U[j][i][1]))#ij调换以和c++版本保持一致（虽然我不是很清楚）

def draw_velocity_field(NX, NY, U, UU, n, save_fig=False):
    x = np.arange(0, NX, 1)
    y = np.arange(0, NY, 1)
    xx, yy = np.meshgrid(x, y)
    fig, ax = plt.subplots()
    ax.contourf(yy, xx, UU, color = "k")
    ax.quiver(yy, xx, U[:,:,0], U[:,:,1])
    if(save_fig):
        fig.savefig('velocity_field_{:08d}.png'.format(n))
    fig.canvas.draw()
    arr = np.asarray(fig.canvas.buffer_rgba())
    return arr
