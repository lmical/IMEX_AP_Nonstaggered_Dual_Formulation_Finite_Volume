import numpy as np
import matplotlib.pyplot as plt

def switching_function(x):
    x0 = 0.133
    x1 = 0.46
    res = np.zeros_like(x)

    # Region 1
    mask1 = x <= x0
    res[mask1] = 1.0 - x[mask1]**14

    # Region 2
    mask2 = (x > x0) & (x < x1)
    res[mask2] = (
        np.exp(1.0 - 1.0 / (1.0 - ((x[mask2] - x0) / (x1 - x0))**2)) *
        ((1.0 - x0**14) - (1.0 - x1)**14) + (1.0 - x1)**14
    )

    # Region 3
    mask3 = x >= x1
    res[mask3] = (1.0 - x[mask3])**14

    return res

def sw4(x):
    x0 = 0.45
    alpha = 2e-2
    x = np.asarray(x)
    result = np.zeros_like(x, dtype=float)
    mask = x > x0
    result[mask] = 1.0 - np.exp(-((x[mask] - x0)**2) / alpha)
    return result    



# Generate and plot the function
x_vals = np.linspace(0, 1., 1000)
y_vals = switching_function(x_vals)
z_vals = sw4(x_vals)

plt.figure(figsize=(8, 4))
plt.plot(x_vals, y_vals, label="Switching_Function(x)", color="red")
plt.plot(x_vals, z_vals, label="Switching_Function_work(x)", color="blue")
plt.xlabel("x")
plt.ylabel("Switching_Function(x)")
plt.title("Plot of Switching_Function with x0=0.133, x1=0.45")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
