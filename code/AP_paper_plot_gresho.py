import os
import numpy as np
import matplotlib.pyplot as plt

# Parameters
eps = 1e-6
CFL = 0.475 #0.2
gmm = 1.4
n_levels = 13  # Number of contour levels
time_scheme=-82

name_time_scheme={-52:"IMEX_DeC2_Prim_Crank_Nicolson_Cons",-82:"IMEX_DeC2_Prim_Explicit_Cons"}

# File paths
folder = f"/home/lorenzo/Dropbox/Lorenzo/Git_ActiveFlux/Codes/code_conservative_Euler_2D/IMEX_HYPERBOLIC_TRICK_RHO_P_STAR_STAGE_DEPENDENT/gresho_longer_time_gamma1.4_128/eps{eps}/space_reconstruction27/time_scheme_{time_scheme}/K0.0/CFL{CFL}/"


print(folder)
# filename = "WC_SOLUTION_OCT_0000000.0000000.dat"
filename = "WC_SOLUTION_OCT_0000001.0000000.dat"

# Check folder
if not os.path.isdir(folder):
    raise ValueError(f"Folder not found: {folder}")

# Load mesh
mesh_files = sorted([f for f in os.listdir(folder) if f.startswith("MESH_")])
x = np.loadtxt(os.path.join(folder, mesh_files[0]), skiprows=1)
y = np.loadtxt(os.path.join(folder, mesh_files[1]), skiprows=1)
X, Y = np.meshgrid(x, y)

# Load solution
data = np.loadtxt(os.path.join(folder, filename), skiprows=1)
ro, u, v, p = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Reshape for 2D plotting
RO = ro.reshape(X.shape)
U = u.reshape(X.shape)
V = v.reshape(X.shape)
P = p.reshape(X.shape)

# Compute Mach number
p0  = 0.5
ro0 = 1.0
MACH = np.sqrt(U**2 + V**2) / np.sqrt(gmm * p0 / ro0) #np.sqrt(gmm * P / RO)

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
cf = ax.contourf(X, Y, MACH, levels=n_levels, cmap='coolwarm')
# cl = ax.contour(X, Y, MACH, levels=n_levels, colors='black', linewidths=0.5)
# ax.clabel(cl, inline=True, fontsize=8, fmt='%.2f')

# Colorbar and labels
cbar = fig.colorbar(cf, ax=ax)
# cbar.set_label(r'$Ma$')
# ax.set_title("Mach Number")
# ax.set_xlabel("x")
# ax.set_ylabel("y")
ls=20
cbar.ax.tick_params(labelsize=ls)  # Choose your desired fontsize
ax.tick_params(axis='both', labelsize=ls)

# Save figure
plt.tight_layout()

if filename.endswith("0000000.0000000.dat"):
    plt.savefig(f"gresho_initial.pdf", bbox_inches="tight")
else:
    plt.savefig("gresho_T1_time_scheme_"+name_time_scheme[time_scheme]+f"_eps{eps}_CFL{CFL}.pdf", bbox_inches="tight")


# plt.show()
