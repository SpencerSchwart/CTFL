import numpy as np
import matplotlib.pyplot as plt

def peskin_delta(r):
    abs_r = np.abs(r)
    if abs_r <= 1:
        return (1/8) * (3 - 2*abs_r + np.sqrt(1 + 4*abs_r - 4*abs_r**2))
    elif 1 < abs_r <= 2:
        return (1/8) * (5 - 2*abs_r - np.sqrt(-7 + 12*abs_r - 4*abs_r**2))
    else:
        return 0

r_values = np.linspace(-3,3,500)
delta_values = np.array([peskin_delta(r) for r in r_values])

plt.plot(r_values, delta_values, label='Peskin Delta Function')
plt.xlabel('Delta Δ')
plt.ylabel('δ_d(r)')
plt.title('Peskin Delta Function')
plt.grid(True)
plt.legend()
plt.show()
