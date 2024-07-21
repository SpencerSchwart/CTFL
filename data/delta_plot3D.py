import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def peskin_delta(r):
    abs_r = np.abs(r)
    if abs_r <= 1:
        return (1/8) * (3 - 2*abs_r + np.sqrt(1 + 4*abs_r - 4*abs_r**2))
    elif 1 < abs_r <= 2:
        return (1/8) * (5 - 2*abs_r - np.sqrt(-7 + 12*abs_r - 4*abs_r**2))
    else:
        return 0

# Generate 2D grid of points
x = np.linspace(-2.5, 2.5, 100)
y = np.linspace(-2.5, 2.5, 100)
X, Y = np.meshgrid(x, y)

# Calculate r and delta(r) for each point on the grid
Z = np.zeros_like(X)
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        r = np.sqrt((X[i, j]-0.5)**2 + (Y[i, j]-0.5)**2)
        Z[i, j] = peskin_delta(r)

# Plot the 3D surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_xlabel('x (# of cells)')
ax.set_ylabel('y (# of cells)')
ax.set_zlabel('δ_d(r)')
ax.set_title('3D Plot of Peskin Delta Function')

plt.show()
