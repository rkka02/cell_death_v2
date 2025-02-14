import torch
import numpy as np
import cv2

# f : tcfile[index]
# data_resolution : file.data_resolution
def resize_tomogram_2d(f, data_resolution, target_resolution=0.1632, mode='mip'):
    fz= data_resolution[0]/target_resolution
    fy= data_resolution[1]/target_resolution
    fx= data_resolution[2]/target_resolution
    
    temp = f
    if mode=='mip':
        temp = np.max(temp, axis=0)
    elif mode=='qpi':
        temp = np.sum(temp, axis=0)
        
    if fx!=1 and fy!=1:
        temp = cv2.resize(temp, dsize=(1+int(temp.shape[0]*fy), 1+int(temp.shape[0]*fx)))

    return temp

# f : tcfile[index]
# data_resolution : file.data_resolution
def resize_tomogram_3d(f, data_resolution, target_resolution=0.1632):
    temp = torch.from_numpy(f)
    temp = temp.unsqueeze(0).unsqueeze(0)
    
    fz= data_resolution[0]/target_resolution
    fy= data_resolution[1]/target_resolution
    fx= data_resolution[2]/target_resolution

    temp = torch.nn.functional.interpolate(
        temp, scale_factor=(fz, fy, fx), mode='trilinear', align_corners=False
    )

    temp = temp.squeeze(0).squeeze(0)
    temp = np.array(temp)
    
    return temp

# f : tcfile[index]
# data_resolution : file.data_resolution
def resize_tomogram_mip(f, data_resolution, target_resolution=0.1632, mode='mip'):
    fy= data_resolution[0]/target_resolution
    fx= data_resolution[1]/target_resolution
    
    temp = f
    if fx!=1 and fy!=1:
        temp = cv2.resize(temp, dsize=(1+int(temp.shape[0]*fy), 1+int(temp.shape[0]*fx)))

    return temp