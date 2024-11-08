import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
from process import FlowDataset
from model import Autoencoder
import matplotlib.pyplot as plt
import numpy as np

file_paths = ['data/data-6000000-5.txt']  # Replace with your actual file paths
dataset = FlowDataset(file_paths)

# Create a DataLoader for the dataset
dataloader = DataLoader(dataset, batch_size=4, shuffle=True)  # Adjust batch size as needed

def train_autoencoder(model, dataloader, num_epochs=100):
    # Set up the loss function and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-4)

    # Set model to training mode
    model.train()

    # Loop over epochs
    for epoch in range(num_epochs):
        running_loss = 0.0
        for i, (inputs, _) in enumerate(dataloader):
            # Move data to the appropriate device (CPU or GPU)
            inputs = inputs.to(device)

            # Zero the parameter gradients
            optimizer.zero_grad()

            # Forward pass
            outputs = model(inputs)

            # Compute loss
            loss = criterion(outputs, inputs)

            # Backward pass and optimization
            loss.backward()
            optimizer.step()

            # Accumulate loss for reporting
            running_loss += loss.item()

        # Print loss at the end of each epoch
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {running_loss / len(dataloader):.4f}')

# Initialize the autoencoder model
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = Autoencoder().to(device)

# Call the training function
train_autoencoder(model, dataloader, num_epochs=2000)

def visualize_reconstruction(model, dataloader, num_images=5):
    model.eval()  # Set model to evaluation mode
    with torch.no_grad():  # Disable gradient computation for evaluation
        for inputs, _ in dataloader: 
            inputs = inputs.to(device)
            outputs = model(inputs)
            
            # Convert tensors to numpy arrays
            inputs_np = inputs.cpu().numpy()
            outputs_np = outputs.cpu().numpy()

            # We assume you have x and y coordinates from your dataset
            x = np.linspace(0, inputs_np.shape[2], inputs_np.shape[2])  # Example for x
            y = np.linspace(0, inputs_np.shape[3], inputs_np.shape[3])  # Example for y
            X, Y = np.meshgrid(x, y)
            X = X.flatten()
            Y = Y.flatten()

            for i in range(num_images):
                fig, ax = plt.subplots(3, 2, figsize=(10, 12))  # 3 rows for u_x, u_y, and p; 2 columns for original and reconstructed

                # Calculate min and max for shared colorbar across original and reconstructed
                u_x_min, u_x_max = min(inputs_np[i][0].min(), outputs_np[i][0].min()), max(inputs_np[i][0].max(), outputs_np[i][0].max())
                u_y_min, u_y_max = min(inputs_np[i][1].min(), outputs_np[i][1].min()), max(inputs_np[i][1].max(), outputs_np[i][1].max())
                p_min, p_max = min(inputs_np[i][2].min(), outputs_np[i][2].min()), max(inputs_np[i][2].max(), outputs_np[i][2].max())

                # Plot original u_x
                c1 = ax[0, 0].scatter(Y, X, c=inputs_np[i][0].flatten(), cmap='jet', s=5, vmin=u_x_min, vmax=u_x_max)
                ax[0, 0].set_title('Original u_x', fontsize=10)
                #ax[0, 0].set_xlabel('x')
                ax[0, 0].set_ylabel('y')

                # Plot reconstructed u_x
                c2 = ax[0, 1].scatter(Y, X, c=outputs_np[i][0].flatten(), cmap='jet', s=5, vmin=u_x_min, vmax=u_x_max)
                ax[0, 1].set_title('Reconstructed u_x', fontsize=10)
                # ax[0, 1].set_xlabel('x')
                # ax[0, 1].set_ylabel('y')

                # Plot original u_y
                c3 = ax[1, 0].scatter(Y, X, c=inputs_np[i][1].flatten(), cmap='jet', s=5, vmin=u_y_min, vmax=u_y_max)
                ax[1, 0].set_title('Original u_y', fontsize=10)
                # ax[1, 0].set_xlabel('x')
                ax[1, 0].set_ylabel('y')

                # Plot reconstructed u_y
                c4 = ax[1, 1].scatter(Y, X, c=outputs_np[i][1].flatten(), cmap='jet', s=5, vmin=u_y_min, vmax=u_y_max)
                ax[1, 1].set_title('Reconstructed u_y', fontsize=10)
                # ax[1, 1].set_xlabel('x')
                # ax[1, 1].set_ylabel('y')

                # Plot original pressure
                c5 = ax[2, 0].scatter(Y, X, c=inputs_np[i][2].flatten(), cmap='jet', s=5, vmin=p_min, vmax=p_max)
                ax[2, 0].set_title('Original Pressure p', fontsize=10)
                ax[2, 0].set_xlabel('x')
                ax[2, 0].set_ylabel('y')

                # Plot reconstructed pressure
                c6 = ax[2, 1].scatter(Y, X, c=outputs_np[i][2].flatten(), cmap='jet', s=5, vmin=p_min, vmax=p_max)
                ax[2, 1].set_title('Reconstructed Pressure p', fontsize=10)
                ax[2, 1].set_xlabel('x')
                # ax[2, 1].set_ylabel('y')

                # Add shared colorbars
                fig.colorbar(c1, ax=ax[0, :], orientation='vertical')
                fig.colorbar(c3, ax=ax[1, :], orientation='vertical')
                fig.colorbar(c5, ax=ax[2, :], orientation='vertical')

                # plt.tight_layout()
                plt.show()

            break  # Only do this for the first batch

# Visualize reconstruction for a few images
visualize_reconstruction(model, dataloader)
