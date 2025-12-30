import os
import numpy as np
import matplotlib.pyplot as plt

# Parameters
eps = 1e-4
CFL = 0.475 #0.1 #0.475 
gmm = 1.4

time_scheme=-82
name_time_scheme={-52:"IMEX_DeC2_Prim_Crank_Nicolson_Cons",-82:"IMEX_DeC2_Prim_Explicit_Cons"}


# File paths
folder = f"/home/lorenzo/Desktop/2_WORK/Research/2024_2026_NCSU/Alina_Alex_IMEX_AP_Nonstaggered_Dual_Formulation_Finite_Volume/code/IMEX_HYPERBOLIC_TRICK_RHO_P_STAR_STAGE_DEPENDENT/double_shear_layer_256/eps"+str(eps)+"/space_reconstruction27/time_scheme_"+str(time_scheme)+"/K0.0/CFL"+str(CFL)+"/"
# filename = "WC_SOLUTION_OCT_0000000.0000000.dat"; T=0
filename = "WC_SOLUTION_OCT_0000006.0000000.dat"; T=6
# filename = "WC_SOLUTION_OCT_0000010.0000000.dat"; T=10

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

VORT = np.zeros_like(RO)
NX = X.shape[0]
NY = Y.shape[1]

dx = X[1, 1] - X[1, 0]
dy = Y[1, 1] - Y[0, 1]

# print(dx, dy)


for indi in range(1, NX - 1):
    for indj in range(1, NY - 1):
        VORT[indi, indj] = (V[indi, indj+1] - V[indi, indj-1]) / (2 * dx) - (U[indi+1, indj] - U[indi-1, indj]) / (2 * dy)

# Plot
fig, axv = plt.subplots(figsize=(8, 6))

# Clear the axis and update the plot with new data
axv.clear()
# Plotting
artificialmin = -5
artificialmax = 5
if artificialmin == artificialmax:
    artificialmax += 1e-8  # or choose a small range manually
levels = np.linspace(artificialmin, artificialmax, 100)
# Contour lines (black lines on top)
# cont_lines = axv.contour(X, Y, VORT, levels=levels, colors='black', linewidths=0.5)




cont_vort = axv.contourf(X, Y, VORT, levels=levels, cmap='jet')
cbar_vort = fig.colorbar(cont_vort, ax=axv)
# axv.set_title(r"$vorticity$")
# axv.set_xlabel('X')
# axv.set_ylabel('Y')

cbar_vort.set_ticks(np.arange(artificialmin, artificialmax+0.01, 1))  # ticks at 0.0, 0.5, ..., 2.0
ls=20
cbar_vort.ax.tick_params(labelsize=ls)  # Choose your desired fontsize
axv.tick_params(axis='both', labelsize=ls)
cbar_vort.set_ticks([-4, -2, 0, 2, 4])



# Save figure
plt.tight_layout()
if filename.endswith("0000000.0000000.dat"):
    plt.title("Initial condition",fontsize=25)
    plt.savefig("double_shear_layer_initial_vorticity.pdf", bbox_inches="tight")
else:
    plt.title(r"$\varepsilon = 10^{%d}$" % int(np.log10(eps)), fontsize=25)
    # plt.savefig("double_shear_layer_at_time_"+str(T)+"_time_scheme_"+name_time_scheme[time_scheme]+f"_eps"+str(eps)+"_CFL"+str(CFL)+".pdf", bbox_inches="tight")
    plt.savefig("double_shear_layer_"+filename+f"_eps"+str(eps)+"_CFL"+str(CFL)+".pdf", bbox_inches="tight")
# plt.show()
