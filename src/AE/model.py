import torch.nn as nn
import torch.nn.functional as F

num_kernel = 5
num_stride = 2
num_padding = 4

class Autoencoder(nn.Module):
    def __init__(self):
        super(Autoencoder, self).__init__()

        # Encoder
        self.encoder = nn.Sequential(
            nn.Conv2d(3, 16, kernel_size=num_kernel, stride=num_stride, padding=num_padding),  # From 3 channels -> 16 channels
            nn.ReLU(),
            nn.Conv2d(16, 32, kernel_size=num_kernel, stride=num_stride, padding=num_padding),  # From 16 channels -> 32 channels
            nn.ReLU(),
            nn.Conv2d(32, 64, kernel_size=num_kernel, stride=num_stride, padding=num_padding),  # From 32 channels -> 64 channels
            nn.ReLU(),
            nn.Conv2d(64, 128, kernel_size=num_kernel, stride=num_stride, padding=num_padding),  # From 32 channels -> 64 channels
            nn.ReLU(),
            nn.Conv2d(128, 256, kernel_size=num_kernel, stride=num_stride, padding=num_padding),  # From 32 channels -> 64 channels
            nn.ReLU()
        )

        # Decoder
        self.decoder = nn.Sequential(
            nn.ConvTranspose2d(256, 128, kernel_size=num_kernel, stride=num_stride, padding=num_padding, output_padding=1),  # Upsample to 64 channels
            nn.ReLU(),
            nn.ConvTranspose2d(128, 64, kernel_size=num_kernel, stride=num_stride, padding=num_padding, output_padding=0),  # Upsample to 64 channels
            nn.ReLU(),
            nn.ConvTranspose2d(64, 32, kernel_size=num_kernel, stride=num_stride, padding=num_padding, output_padding=0),   # Upsample to 32 channels
            nn.ReLU(),
            nn.ConvTranspose2d(32, 16, kernel_size=num_kernel, stride=num_stride, padding=num_padding, output_padding=1),   # Upsample to 16 channels
            nn.ReLU(),
            nn.ConvTranspose2d(16, 3, kernel_size=num_kernel, stride=num_stride, padding=num_padding, output_padding=0),    # Back to 3 channels
            # nn.Sigmoid()  # Output in range [0, 1]
            nn.Tanh()
        )


    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        # x = F.interpolate(x, size=(103, 154), mode='bilinear', align_corners=False)
        return x
