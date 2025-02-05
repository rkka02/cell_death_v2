import torch

def calculate_edge_mean_2d(tensor):
    return torch.mean((tensor[0,0,:] + tensor[0,1,:] + tensor[0,:,0] + tensor[0,:,1])/4)
