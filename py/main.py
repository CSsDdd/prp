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
    controller=simulation.controller(image_path=r"D:\coding\projects\prp\cell_random_img2.png",config_path=r"D:\coding\projects\prp\config3.toml")
    
    max_iter = 10000
    for n in range(max_iter):
        if(n % 5 == 0):#每5步输出一次结果
            err=controller.report_thermal_err(field_tol=1e-6, peak_tol=1e-6,
                           outlet_tol=1e-4, required_checks=3)
            print('Iter: {}, feild_Error: {},peak_Error: {},outlet_Error: {}'.format(n, err["field_rms_error"],err["source_peak_error"],err["outlet_heat_error"]))
            controller.write(n)#只传个编号就好了
            if(err["converged"]):#误差足够小，认为收敛
                print('Converged at Iter: {}, Error: {}'.format(n, err))
                controller.write(n)
                break
        controller.run_one_step()