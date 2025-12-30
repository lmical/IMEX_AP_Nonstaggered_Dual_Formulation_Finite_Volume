import numpy as np


# !*I want 0 for eps=1
# !*I want 1 for eps=0


def final(eps):
    eps = np.asarray(eps)
    x0 = 0.15
    x1 = 0.4

    r=x1-x0

    result = np.zeros_like(eps, dtype=float)
    
    # Region 1: eps <= x0 → result = 1
    mask1 = eps <= x0
    result[mask1] = 1.0-eps[mask1]**14

    # Region 2: x0 < eps < x0 + r → result = 1 / (1 - ((eps - x0)/r)^2)
    mask2 = (eps > x0) & (eps < x1)
    result[mask2] =np.exp( 1.0 - 1.0 / (1.0 - ( ( (eps[mask2] - x0) / (x1-x0) ) ** 2) ) )*( (1.0-x0**14)-(1.0-x1)**14 )+(1.0-x1)**14


    # Region 3: eps >= x0 + r → result = 0 (already set by initialization)
    mask3 = eps >= x1
    result[mask3] = (1-eps[mask3])**14


    return result






x=np.linspace(1.0,10.0,100000)
eps=1/x

import matplotlib.pyplot as plt

plt.plot(x, final(eps), label=r'$s(\varepsilon)$')
plt.xlim([1,10])
plt.xlabel(r"$\frac{1}{\varepsilon}$",fontsize=30)
plt.tick_params(axis='both', labelsize=20)
plt.legend(fontsize=20,loc="lower right")
# plt.grid(True)
plt.tight_layout()
plt.savefig("switching_function.pdf")
# plt.show()


