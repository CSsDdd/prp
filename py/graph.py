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

#这里chance 是各个图形的概率 ，里面每个数组（其实改成tuple会更可读），分别是概率，和对应的形状 ， graph_size是图形半径相对grid大小，size_float是半径大小的波动范围
def cell_noise(size,grid_size=(40,40),seed=42,chance=[[0.2,"circle"],[0.2,"square"],[0.2,"diamond"]],graph_size=0.4/2,size_float=0.1,save_path="cell_random_img.jpg"):
    #np.random.seed(seed)
    map=np.ones(shape=size)
    x_cnt=size[0]//grid_size[0]
    y_cnt=size[1]//grid_size[1]
    gszx=int(grid_size[0]*graph_size)
    gszy=int(grid_size[1]*graph_size)
    print(gszx,gszy)
    for i in range(1,x_cnt-1):#边缘不取
        for j in range(1,y_cnt-1):#边缘不取特征点
            posx=i*grid_size[0]+int(np.random.uniform(0,grid_size[0]))
            posy=j*grid_size[1]+int(np.random.uniform(0,grid_size[1]))
            #print(posx,posy)#随机特征点位置
            op=np.random.uniform(0.,1.)
            tpe="null"
            for k in chance:
                if(op<=k[0]):
                    tpe=k[1]
                    break
                else:
                    op-=k[0]
            sz_index=np.random.uniform(1.-size_float,1.+size_float)
            szx=int(gszx*sz_index)
            szy=int(gszy*sz_index)
            for dx in range(0,szx+1):#考虑到取余问题，多动一位，毕竟是有进一步的距离判断的，问题不大
                for dy in range(0,szy+1):
                    check=None#check判断是否需要着色。可以以此为入手点，进行进一步丰富多种形状！tpe参数描述形状，这里可以加一些其他的比如方形（我这里写的是圆形）
                    if(tpe=="circle"):
                        check=lambda dx,dy: True if np.sqrt(dx**2+dy**2)<=min(szx,szy) else False#
                    elif(tpe=="square"):
                        check=lambda dx,dy: True
                    elif(tpe=="diamond"):
                        check=lambda dx,dy: True if (abs(dx)+abs(dy))<=min(szx,szy) else False
                    elif(tpe=="null"):
                        check=lambda dx,dy: False
                    if check(dx,dy):
                        map[posx+dx][posy+dy]=0
                    if check(-dx,dy):
                        map[posx-dx][posy+dy]=0
                    if check(dx,-dy):
                        map[posx+dx][posy-dy]=0
                    if check(-dx,-dy):
                        map[posx-dx][posy-dy]=0
    map=map*255
    cv2.imwrite(save_path,map)
    return map
cell_noise((400,400))


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

#gri_sin_perlin(max_X_power=3,max_Y_power=3)
#help(noise)