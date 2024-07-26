import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def peskin_delta(x, h):
    abs_x = np.abs(x)
    result = np.zeros_like(x)
    mask = abs_x <= 2 * h
    result[mask] = (1 + np.cos(np.pi * abs_x[mask] / (2 * h))) / (4 * h)
    return result

# Define grid
h = 1.0
x = np.linspace(-3*h, 3*h, 100)
y = np.linspace(-3*h, 3*h, 100)
X, Y = np.meshgrid(x, y)

# Compute Peskin delta function values
Z = peskin_delta(X, h) * peskin_delta(Y, h)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_xlabel('x/h')
ax.set_ylabel('y/h')
ax.set_zlabel('δ')
# ax.set_title('Peskin Delta Function in 3D')

plt.show()
