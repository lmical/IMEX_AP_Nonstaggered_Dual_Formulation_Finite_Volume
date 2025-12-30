import os
import numpy as np
import matplotlib.pyplot as plt

# Parameters
CFL = 0.475
gmm = 1.4

time_scheme=-82
name_time_scheme={-52:"IMEX_DeC2_Prim_Crank_Nicolson_Cons",-82:"IMEX_DeC2_Prim_Explicit_Cons"}


# File paths
folder = f"/home/lorenzo/Desktop/2_WORK/Research/2024_2026_NCSU/Alina_Alex_IMEX_AP_Nonstaggered_Dual_Formulation_Finite_Volume/code/IMEX_HYPERBOLIC_TRICK_RHO_P_STAR_STAGE_DEPENDENT/baroclinic_vorticity_generation/eps0.05/space_reconstruction27/time_scheme_"+str(time_scheme)+"/K0.0/CFL"+str(CFL)+"/"
filename = "WC_SOLUTION_OCT_0000000.0000000.dat"; T=0
# filename = "WC_SOLUTION_OCT_0000010.0000000.dat"; T=10
# filename = "WC_SOLUTION_OCT_0000020.0000000.dat"; T=20

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
MACH = np.sqrt(U**2 + V**2) / np.sqrt(gmm * P / RO)

# Plot
fig, ax_ro = plt.subplots(figsize=(20, 8))
minRO = np.min(RO)
maxRO = np.max(RO)
if minRO == maxRO:
    maxRO += 1e-8  # or choose a small range manually


ax_ro.set_aspect('equal')

print(minRO,maxRO)

# Plotting
artificialmin=0.0
artificialmax=2.0
levels = np.linspace(artificialmin, artificialmax, 50)
# Contour lines (black lines on top)
cont_lines = ax_ro.contour(X, Y, RO, levels=levels, colors='black', linewidths=0.5)






# Optional: label the contour lines
# ax_ro.clabel(cont_lines, inline=True, fontsize=8, fmt="%.2f")



cont_ro = ax_ro.contourf(X, Y, RO, cmap='coolwarm', levels=levels, vmin=artificialmin, vmax=artificialmax)
cbar_ro = fig.colorbar(cont_ro, ax=ax_ro, shrink=0.4)
# Set font size for colorbar tick labels
ls=14
cbar_ro.ax.tick_params(labelsize=ls)  # Choose your desired fontsize
ax_ro.tick_params(axis='both', labelsize=ls)
# Set colorbar tick labels
cbar_ro.set_ticks(np.arange(artificialmin, artificialmax+0.01, 0.2))  # ticks at 0.0, 0.5, ..., 2.0



ax_ro.set_title(r"t="+f"{T}",fontsize=25)
ax_ro.tick_params(axis='both', which='both', labelsize=20)
cbar_ro.ax.tick_params(labelsize=20)
cbar_ro.set_ticks([0.0, 0.5, 1.0, 1.5, 2.0])



# cbar_ro.set_label(r"$\rho$")
# ax_ro.set_title(r"$\rho$")
# ax_ro.set_xlabel('X')
# ax_ro.set_ylabel('Y')

# Save figure
plt.tight_layout()
# plt.savefig("baroclinic_vorticity_generation_at_time_"+str(T)+"_time_scheme_"+name_time_scheme[time_scheme]+f"_eps0.05_CFL{CFL}.pdf", bbox_inches="tight")
# plt.savefig("baroclinic_vorticity_generation_at_time_"+str(T)+f"_eps0.05_CFL{CFL}.pdf", bbox_inches="tight")
plt.savefig("baroclinic_"+filename+f"_eps0.05_CFL{CFL}.pdf", bbox_inches="tight")
# plt.savefig(f"gresho_initial.pdf", bbox_inches="tight")
# plt.show()
