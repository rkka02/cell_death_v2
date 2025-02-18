def crop_patch(slice_2d, crop_size, overlap=False, overlap_proportion=None):
    crop_h, crop_w = crop_size
    height, width = slice_2d.shape

    patches = []
    
    if overlap==True:
        # If stride is not provided, default to 50% overlap
        if overlap_proportion is None:
            overlap_proportion = 0.5
            
        stride_h, stride_w = int(crop_h * (1-overlap_proportion)), int(crop_w * (1-overlap_proportion))
        
        # Loop over the image using the stride
        for top in range(0, height - crop_h + 1, stride_h):
            for left in range(0, width - crop_w + 1, stride_w):
                patches.append((top, left))
    else:
        # Number of vertical/horizontal steps (discard remainder if it doesn't fit exactly)
        num_steps_y = height // crop_h
        num_steps_x = width // crop_w
        # For each tile
        for iy in range(num_steps_y):
            for ix in range(num_steps_x):
                top = iy * crop_h
                left = ix * crop_w
                patches.append((top, left))
                
    return patches