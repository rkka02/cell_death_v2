import torchio
import torchvision.transforms as T
import numpy as np

augmentation_transform_3d = torchio.transforms.Compose([
    torchio.transforms.RandomFlip(axes=(1, 2), flip_probability=0.5),
    torchio.transforms.RandomAffine(degrees=(90, 0, 0), default_pad_value='mean')
])

augmentation_transform_2d = T.Compose([
    T.RandomHorizontalFlip(p=0.5),
    T.RandomVerticalFlip(p=0.5),
    T.RandomRotation(degrees=90),
])