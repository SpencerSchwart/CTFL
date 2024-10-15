import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch, Rectangle

# Load the data from the text file (replace 'log' with your file path)
data = np.loadtxt('log')

# Extract the 2nd and 4th columns, with the correct sign for y
x = data[:, 1]  # 2nd column
y = -data[:, 3]  # 4th column with a negative sign

xmin = 15301
xmax = 16585

# Create the main plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(x, y, 'b-', linewidth=1)

ax.set_xlabel("Time", fontsize=18, family='sans-serif')
ax.set_ylabel(r"$C_D$", fontsize=18, family='sans-serif')

ax.set_xlim(min(x), 150)  # Automatically set x range
ax.set_ylim(0, 2)  # Automatically set y range

ax.set_yticks([0.5, 1, 1.5, 2])

ax.tick_params(axis='both', which='both', direction='in', labelsize=12, length=6)
ax.xaxis.set_ticks_position('both')  # Add ticks to both top and bottom
ax.yaxis.set_ticks_position('both')  # Add ticks to both left and right

# Add an inset plot, zooming in on points between xmin and xmax
ax_inset = fig.add_axes([0.3, 0.2, 0.35, 0.35])  # [left, bottom, width, height]
ax_inset.plot(x, y, 'b-', linewidth=1)
ax_inset.set_xlim(x[xmin], x[xmax])  # Zoom in on x values between xmin and xmax
ax_inset.set_ylim(min(y[xmin:xmax]), max(y[xmin:xmax]))  # Zoom in on y values between xmin and xmax
ax_inset.set_yticks([1.4, 1.45])


# Add a rectangle in the main plot to show where the zoomed-in area is
rect = Rectangle((x[xmin], min(y[xmin:xmax])), x[xmax] - x[xmin], max(y[xmin:xmax]) - min(y[xmin:xmax]), 
                 edgecolor='gray', facecolor='none', linewidth=2)
ax.add_patch(rect)

# Add connecting lines between the inset and the zoomed-in area in the main plot
con = ConnectionPatch(xyA=(x[xmin], max(y[xmin:xmax])), xyB=(x[xmin], max(y[xmin:xmax])), 
                      coordsA="data", coordsB="data", axesA=ax_inset, axesB=ax, color="gray", lw=1)
fig.add_artist(con)

con2 = ConnectionPatch(xyA=(x[xmax], min(y[xmin:xmax])), xyB=(x[xmax], min(y[xmin:xmax])), 
                       coordsA="data", coordsB="data", axesA=ax_inset, axesB=ax, color="gray", lw=1)
fig.add_artist(con2)

# Save and display the plot
plt.savefig('cdvt2.pdf', format='pdf', dpi=300)
plt.show()
