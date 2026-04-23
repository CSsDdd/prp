import cv2
import numpy as np
def read_image(path):
    img = cv2.imread(path)
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mask = cv2.inRange(gray_img, 0, 127)
    return mask

def convert_img_to_mask(img: np.ndarray):
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mask = cv2.inRange(gray_img, 0, 127)
    return mask

def generate_random_img(seed=42, size=(20, 20), save_path='random_img.png'):
    np.random.seed(seed)
    img_tensor = np.random.randint(0, 2, size=size)*255
    cv2.imwrite(save_path, img_tensor)
    return img_tensor

def clip_image(image,size=(40, 40)):
    return cv2.resize(image,size)