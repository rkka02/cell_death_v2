from torch.utils.data import Dataset

class AugmentedDatasetWrapper(Dataset):
    def __init__(self, dataset, num_repeats=3):
        self.dataset = dataset
        self.num_repeats = num_repeats

    def __len__(self):
        # Length is the original dataset length multiplied by the number of repeats
        return len(self.dataset) * self.num_repeats

    def __getitem__(self, idx):
        # Map the index to a sample in the original dataset
        original_idx = idx % len(self.dataset)
        sample, label = self.dataset[original_idx]
        
        # Each time we fetch the sample, it will have an augmentation applied
        return sample, label