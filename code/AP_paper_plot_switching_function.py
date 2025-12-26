import numpy as np
import matplotlib.pyplot as plt

# Define the switching function s(ε)
# Desired behavior:
#   s(1) = 0
#   s(0) = 1

def switching_function(eps):
    eps = np.asarray(eps)
    
    x0 = 0.15
    x1 = 0.4
    r = x1 - x0

    result = np.zeros_like(eps, dtype=float)

    # Region 1: ε ≤ x₀ → s(ε) = 1 - ε¹⁴
    mask1 = eps <= x0
    result[mask1] = 1.0 - eps[mask1]**14

    # Region 2: x₀ < ε < x₁ → smooth transition
    mask2 = (eps > x0) & (eps < x1)
    smooth_term = np.exp(1.0 - 1.0 / (1.0 - ((eps[mask2] - x0) / r)**2))
    a = 1.0 - x0**14
    b = (1.0 - x1)**14
    result[mask2] = smooth_term * (a - b) + b

    # Region 3: ε ≥ x₁ → s(ε) = (1 - ε)¹⁴
    mask3 = eps >= x1
    result[mask3] = (1.0 - eps[mask3])**14

    return result


# Define ε = 1/x over a range of x
x = np.linspace(1.0, 10.0, 100000)
eps = 1.0 / x

# Plot the switching function s(ε) vs 1/ε
plt.plot(x, switching_function(eps), label=r'$s(\varepsilon)$')
plt.xlabel(r"$\frac{1}{\varepsilon}$", fontsize=12)
# plt.ylabel(r"$s(\varepsilon)$", fontsize=12)
# plt.title("Smooth Switching Function", fontsize=14)
plt.xlim([1, 10])
# plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("switching_function.pdf")
plt.show()
