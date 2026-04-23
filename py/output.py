def writetecplot(NX, NY, rho, U, n):
    filename = 'output_{:08d}.dat'.format(n)
    with open(filename, 'w') as f:
        f.write('TITLE = "LBM Simulation"\n')
        f.write('VARIABLES = "X", "Y", "rho", "u", "v"\n')
        f.write('ZONE T="Zone {}", I={}, J={}, F=POINT\n'.format(n//1000, NX, NY))
        for j in range(NY):
            for i in range(NX):
                f.write('{} {} {} {} {}\n'.format(i, j, rho[j][i], U[j][i][0], U[j][i][1]))#ij调换以和c++版本保持一致（虽然我不是很清楚）