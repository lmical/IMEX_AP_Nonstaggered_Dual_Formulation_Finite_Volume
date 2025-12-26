import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from optparse import OptionParser

variable_name="W_X"

parser = OptionParser()
parser.add_option("-d", "--dir", action="store", type="string", dest="d", default="tsunami")
options, args = parser.parse_args()
folder = options.d

aspect_ratio = "quad"

if aspect_ratio == "quad":
    aspect_ratio_dimensions = [4, 4]
elif aspect_ratio == "rect_X":
    aspect_ratio_dimensions = [14, 4]
else: 
    print("Aspect ratio undefined")
    sys.exit(1)

gmm = 1.4

print(folder)
if not os.path.isdir(folder):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

i_plot = 0

cbar_vort = None

def plotting(folder_name):
    global i_plot, cbar_vort
    folder_name = folder_name + "/"
    
    delimiter_in = ' '
    headerlines_in = 1

    # Importing the MESH DATA
    import_mesh = [file for file in os.listdir(folder_name) if file.startswith(variable_name + '_MESH_')]
    import_mesh.sort()
    nfiles_mesh = len(import_mesh)
    for i in range(nfiles_mesh):
        filename = import_mesh[i]
        mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
        if i == 0:
            x = mydata_mesh
            x_ini = min(x)
            x_end = max(x)
        elif i == 1:
            y = mydata_mesh
            y_ini = min(y)
            y_end = max(y)

    X, Y = np.meshgrid(x, y)

    # Importing the SOLUTION DATA
    import_solution = [file for file in os.listdir(folder_name) if file.startswith(variable_name + '_SOLUTION_')]
    import_solution.sort()
    nfiles_solution = len(import_solution)

    fig, axv = plt.subplots(figsize=(8, 6))

    # Initialize colorbars with empty plots
    cont_vort = axv.contourf(X, Y, np.zeros_like(X), cmap='viridis')
    cbar_vort = fig.colorbar(cont_vort, ax=axv)
    cbar_vort.set_label('vort')


    def plot_one_step(i):
        global cbar_vort
        filename = import_solution[i]
        fig.suptitle(filename + " right:next, left:prev, up:-10, down:+10")
        mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)
        ro = mydata_solution[:, 0]
        u = mydata_solution[:, 1]
        v = mydata_solution[:, 2]
        p = mydata_solution[:, 3]
        phi = mydata_solution[:, 4]

        RO = np.reshape(ro, (X.shape[0], Y.shape[1]))
        U = np.reshape(u, (X.shape[0], Y.shape[1]))
        V = np.reshape(v, (X.shape[0], Y.shape[1]))
        P = np.reshape(p, (X.shape[0], Y.shape[1]))
        PHI = np.reshape(phi, (X.shape[0], Y.shape[1]))

        ROU = RO * U
        ROV = RO * V

        ENERGY = P / (gmm - 1) + 0.5 * RO * (U ** 2 + V ** 2)


        VORT = np.zeros_like(RO)
        NX = X.shape[0]
        NY = Y.shape[1]

        dx = X[1, 1] - X[1, 0]
        dy = Y[1, 1] - Y[0, 1]

        # print(dx, dy)


        for indi in range(1, NX - 1):
            for indj in range(1, NY - 1):
                VORT[indi, indj] = (V[indi + 1, indj] - V[indi - 1, indj]) / (2 * dy) - (U[indi, indj + 1] - U[indi, indj - 1]) / (2 * dx)

        # Clear the axis and update the plot with new data
        axv.clear()
        cont_vort = axv.contourf(X, Y, VORT, cmap='viridis', vmin=-5, vmax=5)
        if cbar_vort:
            cbar_vort.remove()
        cbar_vort = fig.colorbar(cont_vort, ax=axv)
        axv.set_title(r"$vorticity$")
        axv.set_xlabel('X')
        axv.set_ylabel('Y')
        axv.set_aspect(aspect_ratio_dimensions[0] / aspect_ratio_dimensions[1])  # Set aspect ratio


        print(np.max(VORT),np.min(VORT))

    plot_one_step(i_plot)

    def change_file(event):
        global i_plot
        sys.stdout.flush()
        if event.key == 'right':
            i_plot += 1
        elif event.key == 'left':
            if i_plot > 0:
                i_plot -= 1
            else:
                return
        elif event.key == 'down':
            i_plot = min(i_plot + 10, nfiles_solution - 1)
        elif event.key == 'up':
            i_plot = max(i_plot - 10, 0)

        plot_one_step(i_plot)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('key_press_event', change_file)
    plt.show()

if __name__ == '__main__':
    plotting(folder)
