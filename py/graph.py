import cv2
import numpy as np
import numba
def read_image(path):
    img = cv2.imread(path)
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mask = cv2.inRange(gray_img, 0, 128)
    return mask

def convert_img_to_mask(img: np.ndarray):
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mask = cv2.inRange(gray_img, 0, 128)
    return mask

def generate_random_img(seed=42, size=(20, 20), save_path='random_img.png'):
    np.random.seed(seed)
    img_tensor = np.random.randint(0, 2, size=size)*255
    cv2.imwrite(save_path, img_tensor)
    return img_tensor

def clip_image(image,NX,NY):
    return cv2.resize(image,dsize=(NX,NY)).T#cv2是反的，

@numba.njit(parallel=True)
def perlin_noise(size,grid_size=(10,10),seed=42):#生成柏林噪声
    #计算网格大小
    x_cnt=size[0]//grid_size[0]+1#x方向多少随机的向量
    y_cnt=size[1]//grid_size[1]+1#y方向多少随机的向量
    #随机生成向量
    angles=np.random.uniform(0,2*np.pi,(x_cnt,y_cnt))#每个网格点的随机向量的角度
    x_dir = np.cos(angles)
    y_dir = np.sin(angles)
    #插值
    #定义插值方法
    fade=lambda t : 6*t**5 - 15*t**4 + 10*t**3 
    res=np.zeros(size)
    for i in range(0,size[0]):
        for j in range(0,size[1]):
            id_x=i//grid_size[0]
            id_y=j//grid_size[1]
            dx=i%grid_size[0]/grid_size[0]
            dy=j%grid_size[1]/grid_size[1]
            n_lt=(x_dir[id_x][id_y]*(dx)+y_dir[id_x][id_y]*(dy))
            n_rt=(x_dir[id_x+1][id_y]*(dx-1)+y_dir[id_x+1][id_y]*(dy))
            n_lb=(x_dir[id_x][id_y+1]*(dx)+y_dir[id_x][id_y+1]*(dy-1))
            n_rb=(x_dir[id_x+1][id_y+1]*(dx-1)+y_dir[id_x+1][id_y+1]*(dy-1))
            u=fade(dx)
            v=fade(dy)
            res[i][j]=v*(u*n_rb+(1-u)*n_lb)+(1-v)*(u*n_rt+(1-u)*n_lt)#插值
    res=(res-res.min())/(res.max()-res.min())
    return res



def gri_sin_perlin(max_X_power=3,max_Y_power=3,size=(400,400),save_path='random_img.png'):
    index=np.zeros(size)
    weight=0.7
    k=1+0.5
    for p in range(0,max_X_power):
        rw=weight*np.random.uniform(0,1)
        for i in range(0,size[0]):#正弦打底
            for j in range(0,size[1]):#正弦打底
                index[i][j]+=np.sin(j*2*np.pi*k/size[1])*rw#
        k=(k-0.5)*2+0.5
        weight*=0.7
    
    weight=1.
    k=1+0.5
    for p in range(0,max_Y_power):
        rw=weight*np.random.uniform(0,1)
        for i in range(0,size[0]):#正弦打底
            for j in range(0,size[1]):#正弦打底
                index[i][j]+=np.sin(i*2*np.pi*k/size[0])*rw#
        k=(k-0.5)*2+0.5
        weight*=0.7
    index=(index-index.min())/(index.max()-index.min())
    index+=perlin_noise(size=(size[0], size[1]), grid_size=(40, 40), seed=42)*0.3
    #归一化
    index=(index-index.min())/(index.max()-index.min())
    #映射
    index=index*255*0.7+255*0.3#垫了一个底让他不那么容易堵住
    cv2.imwrite(save_path,index)
    return index

gri_sin_perlin(max_X_power=3,max_Y_power=3)
#help(noise)