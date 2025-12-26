# vorticity_plot.py

import numpy as np
import matplotlib.pyplot as plt

# === USER INPUTS ===
rho=np.pi/15.0


# Define velocity components
def u_func(x, y):
    if y<=np.pi:
        return np.tanh((y - np.pi/2.0 )/rho)
    else:
        return np.tanh((3.0/2.0*np.pi - y )/rho)

def v_func(x, y):
    return 0.05*np.sin(x)

# Domain and grid
a, b = 0, 2*np.pi   # x-domain
c, d = 0, 2*np.pi   # y-domain
nx, ny = 100, 100      # number of grid points

plot_3d = False        # True = 3D, False = 2D
# ====================

# Create grid
x = np.linspace(a, b, nx)
y = np.linspace(c, d, ny)
X, Y = np.meshgrid(x, y)
hx = (b - a) / (nx - 1)
hy = (d - c) / (ny - 1)

# Initialize u, v arrays
u = np.zeros((ny, nx))
v = np.zeros((ny, nx))

# Fill u and v using loops (no vectorization)
for j in range(ny):
    for i in range(nx):
        u[j, i] = u_func(x[i], y[j])
        v[j, i] = v_func(x[i], y[j])

# Compute vorticity ω = v_x - u_y using central differences
omega = np.zeros_like(u)

for j in range(1, ny - 1):
    for i in range(1, nx - 1):
        dv_dx = (v[j, i+1] - v[j, i-1]) / (2 * hx)
        du_dy = (u[j+1, i] - u[j-1, i]) / (2 * hy)
        omega[j, i] = dv_dx - du_dy

# Plotting
if plot_3d:
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, omega, cmap='coolwarm', edgecolor='none')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('vorticity')
    ax.set_title('3D Vorticity Field')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
else:
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(omega, extent=[a, b, c, d], origin='lower', cmap='coolwarm', aspect='auto')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('2D Vorticity Field')
    fig.colorbar(im, ax=ax)

plt.tight_layout()
plt.show()
