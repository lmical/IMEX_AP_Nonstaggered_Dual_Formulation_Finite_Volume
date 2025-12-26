import os, sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser


variable="X" #"X" #"Y"

folder = os.getcwd()


nametest="sod"
gmm=1.4


print(folder)
if not os.path.isdir(folder):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

i_plot = 0
def plotting(folder_name, what_plot='all'):
    global i_plot
    folder_name = folder_name + "/"
    
    delimiter_in = ' '
    headerlines_in = 1

    files_mesh = folder_name + 'MESH_*.dat'
    files_solution = folder_name + 'SOLUTION_*.dat'

    # Importing the MESH DATA
    import_mesh = [file for file in os.listdir(folder_name) if file.startswith('MESH_')]
    import_mesh.sort()
    nfiles_mesh = len(import_mesh)
    for i in range(nfiles_mesh):
        if i == 0:
            filename = import_mesh[i]
            mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
            x = mydata_mesh
            x_ini = min(x)
            x_end = max(x)
            Nx    = len(x)
        if i == 1:
            filename = import_mesh[i]
            mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
            y = mydata_mesh
            y_ini = min(y)
            y_end = max(y)
            Ny    = len(y)

    X, Y = np.meshgrid(x, y)

    print(Nx,Ny)

    # Importing the SOLUTION DATA
    import_solution = [file for file in os.listdir(folder_name) if file.startswith('SOLUTION_OCT')]
    import_solution.sort()
    nfiles_solution = len(import_solution)

    # fig = plt.figure(1, figsize=(20,8))
    fig = plt.figure(1, figsize=(20,4))


    if what_plot=='density':
        ax_ro = fig.add_subplot(111)
    elif what_plot=='primitive':
        ax_ro = fig.add_subplot(131)
        if variable=="X":
            ax_u  = fig.add_subplot(132)
        elif variable=="Y":
            ax_v  = fig.add_subplot(132)
        else:
            print("Variable not implemented")
            print(variable)
            quit()
        ax_p  = fig.add_subplot(133)
    elif what_plot=='conservative':
        ax_ro = fig.add_subplot(131)
        if variable=="X":
            ax_qx = fig.add_subplot(132)
        elif variable=="Y":
            ax_qy = fig.add_subplot(132)
        else:
            print("Variable not implemented")
            print(variable)
            quit()
        ax_E  = fig.add_subplot(133)
    elif what_plot=='all':
        ax_ro = fig.add_subplot(231)
        if variable=="X":
            ax_u  = fig.add_subplot(232)
        elif variable=="Y":
            ax_v  = fig.add_subplot(232)
        else:
            print("Variable not implemented")
            print(variable)
            quit()
        ax_p  = fig.add_subplot(233)
        if variable=="X":
            ax_qx = fig.add_subplot(234)
        elif variable=="Y":
            ax_qy = fig.add_subplot(234)
        else:
            print("Variable not implemented")
            print(variable)
            quit()

        ax_E  = fig.add_subplot(235)

    else:
        print("Plotting option not available")
        print("It was")
        print(what_plot)
        quit()



    filename = import_solution[-1]
    fig.suptitle(filename + " right:next, left:prev, up:-10, down:+10")
    mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)

    ro       = mydata_solution[:, 0]
    rou      = mydata_solution[:, 1]
    rov      = mydata_solution[:, 2]
    energy   = mydata_solution[:, 3]
    phi      = mydata_solution[:, 4]

    RO       = np.reshape(ro , (X.shape[0], Y.shape[1]))
    ROU      = np.reshape(rou  , (X.shape[0], Y.shape[1]))
    ROV      = np.reshape(rov  , (X.shape[0], Y.shape[1]))
    ENERGY   = np.reshape(energy  , (X.shape[0], Y.shape[1]))
    PHI      = np.reshape(phi, (X.shape[0], Y.shape[1]))

    U=ROU/RO
    V=ROV/RO



    P=(gmm-1.0)*(ENERGY-0.5*RO*(U**2+V**2))


    namevariable=variable
    if variable=="X":
        z=x
        RO       = RO[0,:] 
        U        = U[0,:]  
        V        = V[0,:]  
        P        = P[0,:]  
        ROU      = ROU[0,:]  
        ROV      = ROV[0,:]  
        ENERGY   = ENERGY[0,:]  
        PHI = PHI[0,:]            
    elif variable=="Y":
        z=y
        RO  = RO[:,0] 
        U   = U[:,0]  
        V   = V[:,0]  
        P   = P[:,0]  
        ROU      = ROU[:,0]  
        ROV      = ROV[:,0]  
        ENERGY   = ENERGY[:,0]  
    else:
        print("Variable not implemented",variable)
        quit()

    # Surface plots
    ax_ro.clear()
    ax_ro.plot(z, RO, " .", markersize=2.0, color="magenta")
    ax_ro.set_title(r"$\rho$")
    ax_ro.set_xlabel(namevariable)
    # ax_ro.set_ylabel(r'$\rho$')


    if what_plot=='primitive' or what_plot=='all':

        if variable=="X":
            ax_u.clear()
            ax_u.plot(z, U, " .", markersize=2.0, color="magenta")
            ax_u.set_title('u')
            ax_u.set_xlabel(namevariable)
            # ax_u.set_ylabel('u')
        elif variable=="Y":
            ax_v.clear()
            ax_v.plot(z, V, " .", markersize=2.0, color="magenta")
            ax_v.set_title('v')
            ax_v.set_xlabel(namevariable)
            # ax_v.set_ylabel('v')
        else:
            print("Variable not implemented")
            print(variable)
            quit()

        ax_p.clear()
        ax_p.plot(z, P, " .", markersize=2.0, color="magenta")
        ax_p.set_title('p')
        ax_p.set_xlabel(namevariable)
        # ax_p.set_ylabel('p')

    if what_plot=='conservative' or what_plot=='all':
        if variable=="X":
            ax_qx.clear()
            ax_qx.plot(z, ROU, " .", markersize=2.0, color="magenta")
            ax_qx.set_title(r'$\rho u$')
            ax_qx.set_xlabel(namevariable)
            # ax_qx.set_ylabel(r'$\rho u$')
        elif variable=="Y":
            ax_qy.clear()
            ax_qy.plot(z, ROV, " .", markersize=2.0, color="magenta")
            ax_qy.set_title(r'$\rho v$')
            ax_qy.set_xlabel(namevariable)
            # ax_qy.set_ylabel(r'$\rho v$')
        else:
            print("Variable not implemented")
            print(variable)
            quit()

        ax_E.clear()
        ax_E.plot(z, ENERGY, " .", markersize=2.0, color="magenta")
        ax_E.set_title(r'$E$')
        ax_E.set_xlabel(namevariable)
        # ax_E.set_ylabel(r'$E$')


    if what_plot=='all':

        pass
        # ax_E.set_ylabel(r'$E$')
    
    # plt.show()
    plt.savefig(nametest+"_U.png", format="png", bbox_inches="tight",dpi=300)







if __name__=='__main__':
    
    # plotting(folder,'density')
    plotting(folder,'primitive')
    # plotting(folder,'conservative')
    # plotting(folder,'all')