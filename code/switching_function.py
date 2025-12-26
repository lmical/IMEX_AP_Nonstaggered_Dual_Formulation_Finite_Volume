import numpy as np


def sw1(eps):
    return np.exp(-2000.0*eps**6)

def sw2(eps):
    return np.exp(-3000*eps**8)

def sw3(eps):
    return np.exp(-50000.0*eps**10)

def sw4(eps):
    eps = np.asarray(eps)
    result = np.zeros_like(eps)
    x0=0.45
    alpha=2e-2
    mask = eps > x0
    result[mask] = 1.0-np.exp(-((eps[mask] - x0)) ** 2/ alpha)
    return result



def sw5(eps):
    eps = np.asarray(eps)
    x0 = 0.15
    r = 3e-1

    result = np.zeros_like(eps, dtype=float)
    
    # Region 1: eps <= x0 → result = 1
    mask1 = eps <= x0
    result[mask1] = 1.0

    # Region 2: x0 < eps < x0 + r → result = 1 / (1 - ((eps - x0)/r)^2)
    mask2 = (eps > x0) & (eps < x0 + r)
    result[mask2] = np.exp(1.0 - 1.0 / (1.0 - ((np.abs(eps[mask2] - x0) / r) ** 2)))

    # Region 3: eps >= x0 + r → result = 0 (already set by initialization)

    return result
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




def aaa(eps):
    eps = np.asarray(eps)
    x0 = 0.15
    x1 = 0.4

    r=x1-x0

    result = np.zeros_like(eps, dtype=float)
    
    # Region 1: eps <= x0 → result = 1
    mask1 = eps <= x0
    result[mask1] = 1.0-eps[mask1]**1

    # Region 2: x0 < eps < x0 + r → result = 1 / (1 - ((eps - x0)/r)^2)
    mask2 = (eps > x0) & (eps < x1)
    result[mask2] =np.exp( 1.0 - 1.0 / (1.0 - ( ( (eps[mask2] - x0) / (x1-x0) ) ** 2) ) )*( (1.0-x0**1)-(1.0-x1)**1 )+(1.0-x1)**1


    # Region 3: eps >= x0 + r → result = 0 (already set by initialization)
    mask3 = eps >= x1
    result[mask3] = (1-eps[mask3])**1


    return result




#    return (1.0-eps**4)

# def sw4(eps):
#     return np.exp(-3000*eps**8)*(1.0-eps**8)**2


# # Option 1: Power function
# def f_power(x, alpha=1.0):
#     return 1 - x**alpha

# # Option 2: Logarithmic decay
# def f_log(x):
#     return -np.log(x) / (-np.log(x) + 1)

# # Option 3: Exponential decay (normalized)
# def f_exp(x):
#     return np.exp(-0.1 / (1 - x))



x=np.linspace(1.0,10.0,100000)
eps=1/x

import matplotlib.pyplot as plt

# plt.plot(x, sw1(eps), label='Nan')
# plt.plot(x, sw2(eps), label='current')
# plt.plot(x, sw3(eps), label='ok for smooth')
# plt.plot(x, sw4(eps), label='ok for shocks')
# plt.plot(x, sw5(eps), label='question')
plt.plot(x, final(eps), label='final')
# plt.plot(x, f_power(eps), label='power')
# plt.plot(x, f_log(eps), label='log')
# plt.plot(x, f_exp(eps), label='exp')
# plt.plot(x, aaa(eps), label='aaa')

plt.xlabel("x")
# plt.xlim([9.9999999,10])
# plt.ylim([1-1e-5,1+1e-6])
plt.ylabel("Function values")
plt.legend()
plt.grid(True)
plt.show()

