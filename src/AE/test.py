import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import time
import torch.nn.functional as F
from torchvision import datasets
from torchvision import transforms
from torch.utils.data import DataLoader
from torch.utils.data import SubsetRandomSampler
from torch.utils.data import sampler
import numpy as np
import matplotlib.colors as mcolors

# Hyperparameters
RANDOM_SEED = 49
LEARNING_RATE = 0.0005
BATCH_SIZE = 256
NUM_EPOCHS = 30
NUM_CLASSES = 10

def get_dataloaders_mnist(batch_size, num_workers=0,
                          train_transforms=None, test_transforms=None):
    if train_transforms is None:
        train_transforms = transforms.ToTensor()    
    if test_transforms is None:
        test_transforms = transforms.ToTensor()
    train_dataset = datasets.MNIST(root='data',
                                  train=True,
             transform=train_transforms,download=True)
    valid_dataset = datasets.MNIST(root='data',
                                   train=True,
                                   transform=test_transforms)
    test_dataset = datasets.MNIST(root='data',
                                  train=False,
                                  transform=test_transforms)
    train_loader = DataLoader(dataset=train_dataset,
                              batch_size=batch_size,
                              num_workers=num_workers,
                              shuffle=True)
    test_loader = DataLoader(dataset=test_dataset,
                             batch_size=batch_size,
                             num_workers=num_workers,
                             shuffle=False)
    return train_loader, test_loader

train_loader, valid_loader, test_loader = get_dataloaders_mnist(batch_size=BATCH_SIZE, num_workers=2)

# Checking the dataset
print('Training Set:\n')
for images, labels in train_loader:
    print('Image batch dimensions:', images.size())
    print('Image label dimensions:', labels.size())
    print(labels[:10])