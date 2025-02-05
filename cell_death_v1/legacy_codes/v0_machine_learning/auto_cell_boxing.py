from TCFile import TCFile
import numpy as np
import pandas as pd
from skimage import morphology, filters, measure
from tqdm import tqdm
from pathlib import Path
import h5py
import copy

def get_file_names(directory):
    return [f.name for f in Path(directory).iterdir() if f.is_file()]

if __name__=='__main__':

    # Example usage:

    file_paths = []
    properties = ['centroid', 'area', 'area_convex', 'axis_major_length'] 
    t_list = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36] # 0은 이미 포함

    box_side = 400
    mode_list = ['Apoptosis', 'Necrosis', 'Necroptosis']

    Apoptosis_list = ['compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P12.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P15.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P19.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P20.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P26.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P28.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P29.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P34.TCF',
                'compressed_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P35.TCF'
                ]
    Apoptosis_area_threshold_list = [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000]

    Necrosis_list = ['230721.132040.Necrosis_NaOH.001.Group1.A1.T001P01.TCF',
                 '230725.092931.HeLa_NaOH.002.Group1.A1.T002P01.TCF',
                 '230725.094439.HeLa_NaOH.002.Group1.A1.T002P01.TCF',
                 '230725.095615.HeLa_NaOH.003.Group1.A1.T003P01.TCF',
                 '230725.100700.HeLa_NaOH.004.Group1.A1.T004P01.TCF',
                 '230725.102033.HeLa_NaOH.005.Group1.A1.T005P01.TCF'
                ]
    Necrosis_area_threshold_list = [10000, 10000, 10000, 10000, 10000, 10000]

    Necroptosis_list = ['compressed_230510.174730.HeLa_Hoechst.001.Group1.A2.T001P02.TCF',
                'compressed_230510.174730.HeLa_Hoechst.001.Group1.A2.T001P03.TCF',
                'compressed_230510.174730.HeLa_Hoechst.001.Group1.A2.T001P17.TCF',
                'compressed_230510.174730.HeLa_Hoechst.001.Group2.A1.T001P02.TCF',
                'compressed_230510.174730.HeLa_Hoechst.001.Group2.A1.T001P10.TCF',
                'compressed_230510.174730.HeLa_Hoechst.001.Group2.A1.T001P19.TCF',
                'compressed_230510.174730.HeLa_Hoechst.001.Group2.A1.T001P20.TCF'
                ]
    Necroptosis_area_threshold_list = [30000, 20000, 20000, 20000, 10000, 20000, 20000]

    file_list_list = [Apoptosis_list, Necrosis_list, Necroptosis_list]
    area_threshold_list_list = [Apoptosis_area_threshold_list, Necrosis_area_threshold_list, Necroptosis_area_threshold_list]
    start_layer_list = [27, 33, 28]
    end_layer_list = [33, 39, 34]

    for m_idx, mode in enumerate(mode_list):
        directory = 'C:/rkka_Projects/Cell/Data/' + mode + '/'
        file_list = file_list_list[m_idx]
        area_threshold_list = area_threshold_list_list[m_idx]
        start_layer = start_layer_list[m_idx]
        end_layer = end_layer_list[m_idx]

        for f_idx, file in enumerate(file_list):
            
            area_threshold = area_threshold_list[f_idx]

            tcfile = TCFile(directory + file, '3D')
            print('file loaded')

            # write file
            with h5py.File(directory + '/timelapsed_boxed_' + file + '.h5', 'w') as f:
                
                ###############################
                # base box 만들기
                a = tcfile[0]

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

                properties = ['centroid', 'area'] 

                labels = measure.label(aa)
                props = measure.regionprops_table(labels, properties=properties)
                props = pd.DataFrame(props)
                props = props.drop(props[props['area']<area_threshold].index)

                # Calculate Centroids
                centroids = np.zeros(np.array(props[['centroid-1', 'centroid-2']]).shape, dtype=int)
                for i, c in enumerate(np.array(props[['centroid-1', 'centroid-2']])):
                    centroids[i][0] = int(c[0])
                    centroids[i][1] = int(c[1])

                # Construct Bounding Box

                for k, c in enumerate(centroids):
                    box = np.zeros((labels.shape[0],box_side,box_side)) # (6, 300, 300)

                    h_start = c[0]-box_side//2
                    h_end = c[0]+box_side//2
                    w_start = c[1]-box_side//2
                    w_end = c[1]+box_side//2

                    if h_start <= 0:
                        delta_h = copy.deepcopy(h_start)
                        h_start = 0
                        h_end = h_end + np.abs(delta_h)
                    if h_end >= labels.shape[1]:
                        delta_h = copy.deepcopy(h_end-labels.shape[1])
                        h_end = labels.shape[1]
                        h_start = h_start - np.abs(delta_h)
                    if w_start <= 0:
                        delta_w = copy.deepcopy(w_start)
                        w_start = 0
                        w_end = w_end + np.abs(delta_w)
                    if w_end >= labels.shape[2]:
                        delta_w = copy.deepcopy(w_end-labels.shape[1])
                        w_end = labels.shape[1]
                        w_start = w_start - np.abs(delta_w)
                                
                    box = a[start_layer:end_layer, h_start:h_end, w_start:w_end]
                        
                    # data save
                    f.create_dataset('time/00h/box/' + str(k).zfill(3), data=box)

                ######################################

                # for t in range(tcfile.length):
                for t_idx, t in tqdm(enumerate(t_list)): 
                    a = tcfile[t]
                
                    # Construct Bounding Box
                    for k, c in enumerate(centroids):
                        box = np.zeros((labels.shape[0],box_side,box_side)) # (6, 300, 300)

                        h_start = c[0]-box_side//2
                        h_end = c[0]+box_side//2
                        w_start = c[1]-box_side//2
                        w_end = c[1]+box_side//2

                        if h_start <= 0:
                            delta_h = copy.deepcopy(h_start)
                            h_start = 0
                            h_end = h_end + np.abs(delta_h)
                        if h_end >= labels.shape[1]:
                            delta_h = copy.deepcopy(h_end-labels.shape[1])
                            h_end = labels.shape[1]
                            h_start = h_start - np.abs(delta_h)
                        if w_start <= 0:
                            delta_w = copy.deepcopy(w_start)
                            w_start = 0
                            w_end = w_end + np.abs(delta_w)
                        if w_end >= labels.shape[2]:
                            delta_w = copy.deepcopy(w_end-labels.shape[1])
                            w_end = labels.shape[1]
                            w_start = w_start - np.abs(delta_w)
                                    
                        box = a[start_layer:end_layer, h_start:h_end, w_start:w_end]
                            

                        # data save
                        f.create_dataset('time/' + str(t).zfill(2) + 'h/' + 'box/' + str(k).zfill(3), data=box)