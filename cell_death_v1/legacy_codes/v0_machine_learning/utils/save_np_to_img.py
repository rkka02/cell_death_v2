import numpy as np
from PIL import Image

def save_np_to_img(array, destination_path):
    array = (array - np.min(array))/(np.max(array) - np.min(array)) * 255
    image = Image.fromarray(array.astype(int)).convert('RGB')
    image.save(destination_path)
    return 0