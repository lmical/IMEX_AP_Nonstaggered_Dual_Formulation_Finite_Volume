import numpy as np
import matplotlib.pyplot as plt

# Load the data
filename = "den250"  # Change this to the actual filename
data = np.loadtxt(filename)

# Extract columns
x, y, rho = data[:, 0], data[:, 1], data[:, 2]

# Create a structured grid (assuming a structured Fortran output)
nx = len(np.unique(x))
ny = len(np.unique(y))
X = x.reshape((ny, nx))
Y = y.reshape((ny, nx))
Rho = rho.reshape((ny, nx))

# Plot using a jet colormap
plt.figure(figsize=(8, 6))
plt.pcolormesh(X, Y, Rho, shading='auto', cmap='jet')
plt.colorbar(label='Density')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Density Contour from Fortran Output')
plt.show()

