def crop_patch(slice_2d, crop_size):
    crop_h, crop_w = crop_size
    height, width = slice_2d.shape
    # Number of vertical/horizontal steps (discard remainder if it doesn't fit exactly)
    num_steps_y = height // crop_h
    num_steps_x = width // crop_w
    patches = []
    # For each tile
    for iy in range(num_steps_y):
        for ix in range(num_steps_x):
            top = iy * crop_h
            left = ix * crop_w
            patches.append((top, left))
    return patches