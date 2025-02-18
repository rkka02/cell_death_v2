import numpy as np

def check_patch(patch, threshold=0.1):
    # remove bad quality patch
    patch[patch<45] = 0
    patch[patch>=45] = 1
    proportion = np.count_nonzero(patch==1)/(np.count_nonzero(patch==0)+np.count_nonzero(patch==1))
    if proportion > threshold:
        return True
    else:
        return False