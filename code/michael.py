import numpy as np

eps=np.array([1e-1,1e-2,1e-3,1e-4,1e-5,1e-6])
eps2=eps**2
err=np.array([3.2407994254371057E-006,3.6324976040378850E-008,5.0213086966665547E-010,5.6987299807461506E-012,5.7270818133465379E-014,6.5389588275316379E-016])

import matplotlib.pyplot as plt


plt.loglog(eps2,err)
plt.xlabel("eps**2")
plt.ylabel("error")
plt.show()