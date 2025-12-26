import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

variable_name = "W_X"

# Set up command-line options
parser = OptionParser()
parser.add_option("-d", "--dir", action="store", type="string", dest="d", default="tsunami", help="Directory name")
parser.add_option("-v", "--value", action="store", type="float", dest="value", default=0.01, help="Some value")
parser.add_option("-g", "--gmm", action="store", type="float", dest="gmm", default=1.4, help="GMM value")

options, args = parser.parse_args()

folder = options.d
eps = options.value
gmm = options.gmm

print("Eps:", eps)
print("Gmm:", gmm)

# Validate folder
if not os.path.isdir(folder):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

# Global index for plotting steps
i_plot = 0

def plot_mach(folder_name):
    global i_plot, cbar_Ma  # Declare cbar_Ma as global here

    folder_name = folder_name + "/"
    delimiter_in = ' '
    headerlines_in = 1

    # Import MESH data
    import_mesh = [file for file in os.listdir(folder_name) if file.startswith(variable_name + '_MESH_')]
    import_mesh.sort()
    
    x = np.loadtxt(folder_name + import_mesh[0], delimiter=delimiter_in, skiprows=headerlines_in)
    y = np.loadtxt(folder_name + import_mesh[1], delimiter=delimiter_in, skiprows=headerlines_in)
    
    X, Y = np.meshgrid(x, y)

    # Import SOLUTION data
    import_solution = [file for file in os.listdir(folder_name) if file.startswith(variable_name + '_SOLUTION_OCT')]
    import_solution.sort()
    nfiles_solution = len(import_solution)

    # Setup the plot
    fig, ax_Ma = plt.subplots(figsize=(8, 6))
    cbar_Ma = None  # Initialize colorbar for Mach

    def plot_one_step(i):
        global cbar_Ma  # Declare as global to modify within the function
        filename = import_solution[i]
        fig.suptitle(filename + " right:next, left:prev, up:-10, down:+10")
        
        # Load solution data
        mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)
        u = mydata_solution[:, 1]
        v = mydata_solution[:, 2]
        p = mydata_solution[:, 3]
        
        U = np.reshape(u, (X.shape[0], Y.shape[1]))
        V = np.reshape(v, (X.shape[0], Y.shape[1]))
        P = np.reshape(p, (X.shape[0], Y.shape[1]))

        MACH = np.sqrt(U ** 2 + V ** 2) / np.sqrt(gmm * P / 1.0)

        # Clear previous plot
        ax_Ma.clear()
        
        # Plot Mach number
        cont_Ma = ax_Ma.contourf(X, Y, MACH, cmap='coolwarm')
        if cbar_Ma:
            cbar_Ma.remove()
        cbar_Ma = fig.colorbar(cont_Ma, ax=ax_Ma)
        cbar_Ma.set_label(r'$Ma$')
        ax_Ma.set_title(r'Mach Number ($Ma$)')
        ax_Ma.set_xlabel('X')
        ax_Ma.set_ylabel('Y')

    plot_one_step(i_plot)

    def change_file(event):
        global i_plot
        sys.stdout.flush()
        if event.key == 'right':
            i_plot = min(i_plot + 1, nfiles_solution - 1)
        elif event.key == 'left':
            if i_plot > 0:
                i_plot -= 1
        elif event.key == 'down':
            i_plot = min(i_plot + 10, nfiles_solution - 1)
        elif event.key == 'up':
            i_plot = max(i_plot - 10, 0)

        plot_one_step(i_plot)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('key_press_event', change_file)
    plt.show()

if __name__ == '__main__':
    plot_mach(folder)
