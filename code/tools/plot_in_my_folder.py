import numpy as np
import matplotlib.pyplot as plt

# Step 1: Load the coordinate data
x_coords = np.loadtxt("MESH_CoordinateX.dat", skiprows=1)  # Skip the first line (header)
y_coords = np.loadtxt("MESH_CoordinateY.dat", skiprows=1)  # Skip the first line (header)

# Step 2: Load the solution (density or other variables) data
data = np.loadtxt("SOLUTION_OCT_0000002.5000000.dat", skiprows=1)  # Skip the first line (header)

# Step 3: Extract only the Density column (first column in the file)
density = data[:, 0]  # First column corresponds to Density

# Step 4: Assuming nElemsX and nElemsY are the grid dimensions
nElemsX = len(x_coords)
nElemsY = len(y_coords)

# Step 5: Reshape the density data to match the grid dimensions
# Make sure the total number of elements matches nElemsX * nElemsY
if len(density) == nElemsX * nElemsY:
    density = density.reshape((nElemsY, nElemsX))  # Shape should match (nElemsY, nElemsX)
else:
    raise ValueError(f"Data size does not match the grid dimensions. Expected {nElemsX * nElemsY} values, but got {len(density)}")

# Step 6: Create the meshgrid for plotting
X, Y = np.meshgrid(x_coords, y_coords)

# Step 7: Plot using pcolormesh
plt.figure(figsize=(8, 6))
plt.pcolormesh(X, Y, density, shading='auto', cmap='jet')  # Use 'jet' colormap
plt.colorbar(label='Density')  # Add color bar with label
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Density Plot using pcolormesh')
plt.show()
