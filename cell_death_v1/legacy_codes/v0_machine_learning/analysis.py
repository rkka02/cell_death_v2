from TCFile import TCFile
import numpy as np
import h5py

from skimage import morphology
from skimage import filters
from segment_anything import SamAutomaticMaskGenerator, sam_model_registry

def rescale(data, min=1.33, max=1.4):
    dmax = np.max(data)
    dmin = np.min(data)
    data = min + (data - dmin)/(dmax-dmin)
    return data

def get_thresholds(data, classes=5):
    return filters.threshold_multiotsu(data, classes)

def get_edge(mask):
    mm = mask
    mm = filters.gaussian(mm)
    mm = filters.sobel(mm)
    thr = filters.threshold_multiotsu(mm, 3)
    mm[mm<thr[1]] = False
    mm[mm>thr[1]] = True

    mm = morphology.dilation(morphology.erosion(mm))
    mm = morphology.dilation(mm)
    mm = morphology.dilation(mm)
    mm = morphology.dilation(mm)
    mm = morphology.dilation(mm)
    mm = morphology.dilation(mm)
    mm = morphology.dilation(mm)

    mm = morphology.remove_small_objects(mm.astype(bool), 2048)
    mm = filters.sobel(mm)
    mm = morphology.remove_small_objects(mm.astype(bool), 512)

    return mm

def make_mask(data, threshold, preprocess=False):

    mask = np.zeros(data.shape, dtype=int)
    mask[data > threshold] = 1

    if preprocess==True:
        mask = morphology.erosion(mask)
        mask = morphology.dilation(mask)

        mask = morphology.remove_small_holes(mask)
        # Small objects - fragments 들을 제거할지 말지?
        mask = morphology.remove_small_objects(mask, min_size=128)

    # for masking process
    mask = mask.astype(int)
    return mask

def make_array_3ch(data):
    return np.stack((data,)*3, axis=-1)

def segment_anything(img):

    sam = sam_model_registry["default"](checkpoint=r"C:\rkka_Projects\Cell\sam_checkpoint\sam_vit_h_4b8939.pth")
    mask_generator = SamAutomaticMaskGenerator(sam)
    sam.to(device='cuda')

    masks = mask_generator.generate(img)

    return masks
