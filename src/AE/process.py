import numpy as np
import torch
from torch.utils.data import Dataset
import matplotlib.pyplot as plt


def plot_tensor(input_tensor, x, y):
    # Convert PyTorch tensor to NumPy array
    input_numpy = input_tensor.cpu().numpy()  # Shape: [3, height, width]

    # Extract the three channels: u_x, u_y, p
    u_x = input_numpy[0].flatten()  # First channel, flatten to 1D array
    u_y = input_numpy[1].flatten()  # Second channel
    p = input_numpy[2].flatten()    # Third channel

    # Plot u_x using scatter
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    sc1 = plt.scatter(x, y, c=u_x, cmap='jet', marker='o')
    plt.colorbar(sc1)
    plt.title('u_x')

    # Plot u_y using scatter
    plt.subplot(1, 3, 2)
    sc2 = plt.scatter(x, y, c=u_y, cmap='jet', marker='o')
    plt.colorbar(sc2)
    plt.title('u_y')

    # Plot p using scatter
    plt.subplot(1, 3, 3)
    sc3 = plt.scatter(x, y, c=p, cmap='jet', marker='o')
    plt.colorbar(sc3)
    plt.title('Pressure p')

    plt.show()



def load_and_process_data(file_path):
    # Step 1: Read the text file
    data = np.loadtxt(file_path, skiprows=1)  # Skip the header row

    # Step 2: Extract the relevant columns
    x = data[:, 0]    # x-coordinates (you may not need this for the network)
    y = data[:, 1]    # y-coordinates (you may not need this for the network)
    u_x = data[:, 2]  # u_x velocity component
    u_x /= np.max(np.abs(u_x))
    u_y = data[:, 3]  # u_y velocity component
    u_y /= np.max(np.abs(u_y))
    p = data[:, 4]    # pressure
    p /= np.max(np.abs(p))

    # Step 3: Reshape the data
    # Assuming the data is on a uniform grid and you know the grid dimensions
    # For example, reshape into (height, width) grid
    # grid_size = int(np.sqrt(len(u_x)))  # Assuming square grid
    # num_points = len(u_x)
    # grid_size = int(np.sqrt(num_points))
    # print(grid_size)
    #u_x = u_x.reshape(grid_size, grid_size)  # Trim excess data if needed
    #u_y = u_y.reshape(grid_size, grid_size)
    #p = p.reshape(grid_size, grid_size)

    # Step 3: Calculate the actual grid dimensions
    x_unique = len(np.unique(x))
    y_unique = len(np.unique(y))

    # Step 4: Reshape the data into (y_unique, x_unique)
    u_x = u_x.reshape(y_unique, x_unique)
    u_y = u_y.reshape(y_unique, x_unique)
    p = p.reshape(y_unique, x_unique)

    # print(u_x.shape, u_y.shape, p.shape)
    # plt.scatter(x, y, c=u_x, cmap='jet')
    # plt.imshow(u_x, cmap='jet')
    # plt.pcolormesh(x, y, u_x, cmap='jet', shading='auto')
    # plt.colorbar()
    # plt.show()

    # Step 4: Stack the channels (u_x, u_y, p) together as needed for NN input
    # Stack into a 3-channel input: [u_x, u_y, p]
    input_data = np.stack([u_x, u_y, p], axis=0)  # Shape: [3, height, width]
    # print(input_data.shape)
    # Step 5: Convert to a PyTorch tensor
    input_tensor = torch.tensor(input_data, dtype=torch.float32)
    # plot_tensor(input_tensor, x, y)

    return input_tensor

# Example usage:
# file_path = 'data/data-re40-1.txt'  # Replace with the actual file path

class FlowDataset(Dataset):
    def __init__(self, file_paths):
        self.file_paths = file_paths
    
    def __len__(self):
        return len(self.file_paths)
    
    def __getitem__(self, idx):
        # Load and process the data
        file_path = self.file_paths[idx]
        data_tensor = load_and_process_data(file_path)  # Function from earlier

        return data_tensor, data_tensor  # Input and target are the same for autoencoders

# Example usage
# file_paths = ['your_file1.txt', 'your_file2.txt']  # Replace with your actual file paths
# dataset = FlowDataset(file_paths)

def plot_raw_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)
    x = data[:, 0]
    y = data[:, 1]
    u_x = data[:, 2]

    plt.scatter(x, y, c=u_x, cmap='jet')
    plt.colorbar()
#    plt.show()

# plot_raw_data('data/data.txt')
input_tensor = load_and_process_data('data/data.txt')

def plot_tensor(input_tensor):
    # Convert PyTorch tensor to NumPy array
    input_numpy = input_tensor.cpu().numpy()  # Shape: [3, height, width]

    # Extract the three channels: u_x, u_y, p
    u_x = input_numpy[0]  # First channel
    u_y = input_numpy[1]  # Second channel
    p = input_numpy[2]    # Third channel

    # Plot u_x
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.imshow(u_x, cmap='jet')
    plt.colorbar()
    plt.title('u_x')

    # Plot u_y
    plt.subplot(1, 3, 2)
    plt.imshow(u_y, cmap='jet')
    plt.colorbar()
    plt.title('u_y')

    # Plot p
    plt.subplot(1, 3, 3)
    plt.imshow(p, cmap='jet')
    plt.colorbar()
    plt.title('Pressure p')

    plt.show()

# Example usage:
# plot_tensor(input_tensor)