import torch
"""
def calculate_edge_mean_2d(tensor):
    return torch.mean((tensor[0,0,:] + tensor[0,1,:] + tensor[0,:,0] + tensor[0,:,1])/4)
"""

def calculate_edge_mean_2d(tensor):
    # tensor: shape (1, H, W)
    _, H, W = tensor.shape
    top = tensor[0, 0, :]       # shape: (W,)
    bottom = tensor[0, H-1, :]   # shape: (W,)
    left = tensor[0, :, 0]       # shape: (H,)
    right = tensor[0, :, W-1]    # shape: (H,)
    # Concatenate all border pixels and average
    edge_pixels = torch.cat([top, bottom, left, right])
    return edge_pixels.mean()
