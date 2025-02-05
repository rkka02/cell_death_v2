import numpy as np

def image_normalization(image, min=None, max=None):
    # Min-max normalization (float64 -> uint8) : like png file
    if min is None and max is None:
        min_val, max_val = image.min(), image.max()
        if max_val - min_val < 1e-10:
            image = np.zeros_like(image, dtype=np.uint8)
        else:
            image = ((image - min_val) / (max_val - min_val) * 255).astype(np.uint8)
    else:
        min_val, max_val = min, max
        if max_val - min_val < 1e-10:
            image = np.zeros_like(image, dtype=np.uint8)
        else:
            image = ((image - min_val) / (max_val - min_val) * 255).astype(np.uint8)
        
    return image
