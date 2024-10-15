import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file_path = '/home/spencer/basilisk/CTFL/src/IBM-test/cylindernew6/freq-10-300.dat'
data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['count', 'x', 'y', 'time', 'u_x', 'u_y', 'p'])

time = data['time'].values
u_y = data['u_y'].values

Fs = 1 / (time[1] - time[0])
L = len(time)
T = 1 / Fs

X = np.fft.fft(u_y)

f = np.fft.fftfreq(L, T)

plt.plot(f[:L//2], 2.0/L * np.abs(X[:L//2]))

plt.xlim(0, 2)
plt.ylim(0, 1)

plt.xlabel('Vortex Shredding Frequency (Hz)')
plt.ylabel('Magnitude')
plt.title('x = 10*D')
plt.show()
