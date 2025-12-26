import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser



parser = OptionParser()
parser.add_option("-d", "--dir", action="store", type="string", dest="d", default="tsunami", help="Directory name")
parser.add_option("-v", "--value", action="store", type="float", dest="value", default=1.0, help="Some value")
parser.add_option("-g", "--gmm", action="store", type="float", dest="gmm", default=1.4, help="GMM value")  # New line

options, args = parser.parse_args()

folder = options.d
eps = options.value
gmm = options.gmm  # Update gmm with the value from the command line

print("Eps",eps)
print("Gmm",gmm)

print(folder)
if not os.path.isdir(folder):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

i_plot = 0

# Global colorbar variables
cbar_ro = cbar_u = cbar_v = cbar_p = cbar_qx = cbar_qy = cbar_E = cbar_Ma = None

def plotting(folder_name, what_plot='all'):
    global i_plot, cbar_ro, cbar_u, cbar_v, cbar_p, cbar_qx, cbar_qy, cbar_E, cbar_Ma
    folder_name = folder_name + "/"

    delimiter_in = ' '
    headerlines_in = 1

    # Importing the MESH DATA
    import_mesh = [file for file in os.listdir(folder_name) if file.startswith('MESH_')]
    import_mesh.sort()
    nfiles_mesh = len(import_mesh)
    for i in range(nfiles_mesh):
        if i == 0:
            filename = import_mesh[i]
            mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
            x = mydata_mesh
        if i == 1:
            filename = import_mesh[i]
            mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
            y = mydata_mesh

    X, Y = np.meshgrid(x, y)

    # Importing the SOLUTION DATA
    import_solution = [file for file in os.listdir(folder_name) if file.startswith('WC_Err_SOLUTION_OCT')]
    import_solution.sort()
    nfiles_solution = len(import_solution)

    if what_plot == 'density':
        fig, axs = plt.subplots(1, 1, figsize=(30, 5))
    elif what_plot == 'primitive':
        fig, axs = plt.subplots(1, 4, figsize=(30, 5))
    elif what_plot == 'conservative':
        fig, axs = plt.subplots(1, 4, figsize=(30, 5))
    else:
        fig, axs = plt.subplots(2, 4, figsize=(30, 10))


    if what_plot=="density":
        ax_ro = axs
    elif what_plot=="primitive":
        ax_ro = axs[0]
        ax_u =  axs[1]
        ax_v =  axs[2]
        ax_p =  axs[3]
    elif what_plot=="conservative":
        ax_ro = axs[0]
        ax_qx =  axs[1]
        ax_qy =  axs[2]
        ax_E  =  axs[3]
    elif what_plot=="all":
        ax_ro = axs[0, 0]
        ax_u = axs[0, 1] 
        ax_v = axs[0, 2] 
        ax_p = axs[0, 3] 
        ax_qx = axs[1, 0] 
        ax_qy = axs[1, 1] 
        ax_E  = axs[1, 2] 
        ax_Ma = axs[1, 3] 

    # Initialize colorbars with empty plots
    cont_ro = ax_ro.contourf(X, Y, np.zeros_like(X), cmap='Blues')
    cbar_ro = fig.colorbar(cont_ro, ax=ax_ro)
    cbar_ro.set_label(r"$\rho$")

    if what_plot == 'primitive' or what_plot == 'all':
        cont_u = ax_u.contourf(X, Y, np.zeros_like(X), cmap='viridis')
        cbar_u = fig.colorbar(cont_u, ax=ax_u)
        cbar_u.set_label('u')

        cont_v = ax_v.contourf(X, Y, np.zeros_like(X), cmap='viridis')
        cbar_v = fig.colorbar(cont_v, ax=ax_v)
        cbar_v.set_label('v')

        cont_p = ax_p.contourf(X, Y, np.zeros_like(X), cmap='viridis')
        cbar_p = fig.colorbar(cont_p, ax=ax_p)
        cbar_p.set_label('p')

    if what_plot == 'conservative' or what_plot == 'all':
        if ax_qx:
            cont_qx = ax_qx.contourf(X, Y, np.zeros_like(X), cmap='viridis')
            cbar_qx = fig.colorbar(cont_qx, ax=ax_qx)
            cbar_qx.set_label(r'$q_x$')

        if ax_qy:
            cont_qy = ax_qy.contourf(X, Y, np.zeros_like(X), cmap='viridis')
            cbar_qy = fig.colorbar(cont_qy, ax=ax_qy)
            cbar_qy.set_label(r'$q_y$')

        if ax_E:
            cont_E = ax_E.contourf(X, Y, np.zeros_like(X), cmap='viridis')
            cbar_E = fig.colorbar(cont_E, ax=ax_E)
            cbar_E.set_label(r'$E$')

    if what_plot == 'all':
        if ax_Ma:
            cont_Ma = ax_Ma.contourf(X, Y, np.zeros_like(X), cmap='viridis')
            cbar_Ma = fig.colorbar(cont_Ma, ax=ax_Ma)
            cbar_Ma.set_label(r'$Ma$')


    def plot_one_step(i):
        global cbar_ro, cbar_u, cbar_v, cbar_p, cbar_qx, cbar_qy, cbar_E, cbar_Ma
        filename = import_solution[i]
        fig.suptitle(filename + " right:next, left:prev, up:-10, down:+10")
        mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)
        ro     = mydata_solution[:, 0]
        rou    = mydata_solution[:, 1]
        rov    = mydata_solution[:, 2]
        energy = mydata_solution[:, 3]
        phi    = mydata_solution[:, 4]

        RO     = np.reshape(ro,  (X.shape[0], Y.shape[1]))
        ROU    = np.reshape(rou,   (X.shape[0], Y.shape[1]))
        ROV    = np.reshape(rov,   (X.shape[0], Y.shape[1]))
        ENERGY = np.reshape(energy,   (X.shape[0], Y.shape[1]))
        PHI    = np.reshape(phi, (X.shape[0], Y.shape[1]))

        U=ROU/RO
        V=ROV/RO

        P=(gmm-1.0)*(ENERGY-0.5*RO*(U**2+V**2)*eps**2)

        MACH = np.sqrt(U ** 2 + V ** 2)/np.sqrt(gmm*P/1.0)


        # Clear previous plots and colorbars
        if what_plot!="density":
            for ax in axs.flatten():
                ax.clear()
        else:
            axs.clear()

        # Plotting
        RO=ROU+ROV
        C=25
        avg=np.average(RO)
        print(avg)

        cont_ro = ax_ro.contourf(X, Y, RO, cmap='Blues')
        if cbar_ro:
            cbar_ro.remove()
        cbar_ro = fig.colorbar(cont_ro, ax=ax_ro)
        cbar_ro.set_label(r"$\rho$")
        ax_ro.set_title(r"$\rho$")
        ax_ro.set_xlabel('X')
        ax_ro.set_ylabel('Y')

        if what_plot == 'primitive' or what_plot == 'all':
            cont_u = ax_u.contourf(X, Y, U, cmap='viridis')
            if cbar_u:
                cbar_u.remove()
            cbar_u = fig.colorbar(cont_u, ax=ax_u)
            cbar_u.set_label('u')
            ax_u.set_title('u')
            ax_u.set_xlabel('X')
            ax_u.set_ylabel('Y')

            cont_v = ax_v.contourf(X, Y, V, cmap='viridis')
            if cbar_v:
                cbar_v.remove()
            cbar_v = fig.colorbar(cont_v, ax=ax_v)
            cbar_v.set_label('v')
            ax_v.set_title('v')
            ax_v.set_xlabel('X')
            ax_v.set_ylabel('Y')

            cont_p = ax_p.contourf(X, Y, P, cmap='viridis')
            if cbar_p:
                cbar_p.remove()
            cbar_p = fig.colorbar(cont_p, ax=ax_p)
            cbar_p.set_label('p')
            ax_p.set_title('p')
            ax_p.set_xlabel('X')
            ax_p.set_ylabel('Y')

        if what_plot == 'conservative' or what_plot == 'all':
            if ax_qx:
                cont_qx = ax_qx.contourf(X, Y, ROU, cmap='viridis')
                if cbar_qx:
                    cbar_qx.remove()
                cbar_qx = fig.colorbar(cont_qx, ax=ax_qx)
                cbar_qx.set_label(r'$q_x$')
                ax_qx.set_title(r'$q_x$')
                ax_qx.set_xlabel('X')
                ax_qx.set_ylabel('Y')

            if ax_qy:
                cont_qy = ax_qy.contourf(X, Y, ROV, cmap='viridis')
                if cbar_qy:
                    cbar_qy.remove()
                cbar_qy = fig.colorbar(cont_qy, ax=ax_qy)
                cbar_qy.set_label(r'$q_y$')
                ax_qy.set_title(r'$q_y$')
                ax_qy.set_xlabel('X')
                ax_qy.set_ylabel('Y')

            if ax_E:
                cont_E = ax_E.contourf(X, Y, ENERGY, cmap='viridis')
                if cbar_E:
                    cbar_E.remove()
                cbar_E = fig.colorbar(cont_E, ax=ax_E)
                cbar_E.set_label(r'$E$')
                ax_E.set_title(r'$E$')
                ax_E.set_xlabel('X')
                ax_E.set_ylabel('Y')

        if what_plot == 'all':
            if ax_Ma:
                cont_Ma = ax_Ma.contourf(X, Y, MACH, cmap='viridis')
                if cbar_Ma:
                    cbar_Ma.remove()
                cbar_Ma = fig.colorbar(cont_Ma, ax=ax_Ma)
                cbar_Ma.set_label(r'$Ma$')
                ax_Ma.set_title(r'$Ma$')
                ax_Ma.set_xlabel('X')
                ax_Ma.set_ylabel('Y')


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
    plotting(folder,'density')
    # plotting(folder,'primitive')
    # plotting(folder,'conservative')
    # plotting(folder,'all')
