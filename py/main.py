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
    controller=simulation.controller(image_path=r"pics\cell_random_img2.png",config_path=r"configs\config.toml")
    
    max_iter = 40000
    for n in range(max_iter):
        if(n % 100 == 0):#每200步输出一次结果
            controller.write(n)
        controller.run_one_step()
        reliability=controller.report_reliability(
            velocity_tol=1e-4, temperature_tol=1e-4, peak_tol=1e-4,
            required_checks=3)
        if(n % 100 == 0):
            print('Iter: {}, velocity_Change: {}, temperature_Change: {}, '
                  'peak_Change: {}'.format(
                      n, reliability["velocity_field_relative_change"],
                      reliability["temperature_field_relative_change"],
                      reliability["source_peak_temperature_relative_change"]))

        if(reliability["converged"]):#误差足够小，认为收敛
            if(not n%100 == 0):
                controller.write(n)
            print('Converged at Iter: {}, Reliability: {}'.format(n, reliability))
            break
        
    controller.final_report(reliability)
