from TCFile import TCFile
import numpy as np
import pandas as pd
import napari
import analysis
from skimage import morphology, filters, measure
from tqdm import tqdm
from pathlib import Path
import h5py
import matlab.engine

def get_file_names(directory):
    return [f.name for f in Path(directory).iterdir() if f.is_file()]

if __name__=='__main__':

    # Example usage:

    file_paths = []
    properties = ['area', 'area_convex', 'axis_major_length', 'euler_number', 'extent',
                'inertia_tensor_eigvals', 'solidity'] 
    t_list = [0, 12, 24, 36]

    directory = ['C:/rkka_Projects/Cell/Data/Apoptosis']

    for dir in directory:
        file_names = get_file_names(dir)
        for fn in file_names:
            file = dir + '/' + fn
            
            tcfile = TCFile(file, '3D')
            print('loaded ' + file)

            # write file
            with h5py.File(dir + '/timelapsed_' + fn + '.h5', 'w') as f:

                # for t in range(tcfile.length):
                for t_idx, t in tqdm(enumerate(t_list)): 
                    a = tcfile[t]
                    # set focal plane range
                    start_layer = 32
                    end_layer = 37 + 1

                    mask = np.zeros(a[start_layer:end_layer].shape)
                    mask = mask + a[start_layer:end_layer]
                    thr = filters.threshold_multiotsu(mask, 4)
                    mask[mask<thr[1]] = 0

                    mask = filters.gaussian(mask)
                    
                    cutoff = 1/10

                    for i in tqdm(range(mask.shape[0])):
                        a_fourier = np.fft.fftshift(np.fft.fft2(mask[i]))
                        H, W = a_fourier.shape
                        a_fourier[H//2:H//2 + 1, W//2:W//2 + 1]=0
                        a_fourier[0:int(H*cutoff), :] = 0
                        a_fourier[H-int(H*cutoff):-1, :] = 0
                        a_fourier[:, 0:int(W*cutoff)] = 0
                        a_fourier[:, W-int(W*cutoff):-1] = 0
                        a_mod = np.fft.ifft2(a_fourier)
                        a_mod = np.abs(a_mod)

                        # 푸리에 필터 적용한거에 다시 threshold 걸어주기
                        thr = filters.threshold_multiotsu(a_mod, 4)
                        a_mod[a_mod<thr[0]] = 0

                        # gaussian 적용후 sobel 걸어주면 좀 더 잘됨
                        # a_mod = filters.gaussian(a_mod)
                        a_mod = filters.sobel(a_mod)
                        a_mod = morphology.area_closing(a_mod, area_threshold=10000)

                        # 다시 background 제거
                        a_mod_thr = filters.threshold_multiotsu(a_mod, 3)
                        a_mod[a_mod<a_mod_thr[0]] = 0
                        a_mod = morphology.area_closing(a_mod, area_threshold=50000)

                        mask[i] = a_mod

                    mask = mask.astype(bool)

                    aa = np.zeros(a[start_layer:end_layer].shape)
                    aa = aa + a[start_layer:end_layer]
                    aa[mask==0] = 0

                    labels = measure.label(aa)
                    props = measure.regionprops_table(labels, properties=properties)
                    props = pd.DataFrame(props)

                    # data save
                    f.create_dataset('labels/' + str(t_idx*6).zfill(2) + 'h', data=labels)
                    
                    props.to_csv(dir + '/' + fn + '(' + str(t) + ')' + '.csv', index=False)
