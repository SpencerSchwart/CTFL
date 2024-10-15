import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

class FlowFiedDataset(Dataset):
    def __init__ (self, file_path):
        data = np.loadtxt(file_path, skiprows=1)

        self.x = data[:, 0]
        self.y = data[:, 1]
        self.u_x = data[:, 2]
        self.u_y = data[:, 3]
        self.p = data[:, 4]

        grid_size = int(np.sqrt(len(self.x)))
        self.u_x = self.u_x.reshape((grid_size, grid_size))
        self.u_y = self.u_y.reshape((grid_size, grid_size))
        self.p = self.p.reshape((grid_size, grid_size))

        self.velocity = np.stack([self.u_x, self.u_y], axis=0)

    def __len__(self):
        return 1
    
    def __getitem__(self, idx):
        return torch.tensor(self.velocity, dtype=torch.float32), torch.tensor(self.p, dtype=torch.float32)
    
dataset = FlowFiedDataset('data.txt')

dataloader = DataLoader(dataset, batch_size=1, shuffle=True)
