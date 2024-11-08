import numpy as np
import matplotlib.pyplot as plt

# Create a range of values from -10 to 10
x = np.linspace(-10, 10, 100)

# Apply the ReLU function
relu = np.maximum(0, x)
sigmoid = 1/(1+np.exp(-x))
tanh = 2/(1+np.exp(-2*x))-1

# Plot
plt.plot(x, relu, label='ReLU')
plt.plot(x, sigmoid, label='Sigmoid')
plt.plot(x, tanh, label='Tanh')
plt.xlabel('Input')
plt.ylabel('Output')
plt.title('ReLU Activation Function')
plt.grid(True)
plt.legend()
plt.show()
