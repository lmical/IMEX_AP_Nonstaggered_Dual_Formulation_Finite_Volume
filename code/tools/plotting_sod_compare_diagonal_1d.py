teststocompare=[] 

coordinate="X" #"X" #"Y" #NB:Actually y is not implemented

teststocompare.append(["AF/PP1/M2SSPRK2", "SAFFVO MUSCL"   ,"-","#1f77b4"]) 
teststocompare.append(["Rusanov/M2SSPRK2", "FV MUSCL"   ,"-","#ff7f0e"]) 


colors = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#7b94c7", "#fdae6b", "#b3e2cd", "#ff9896", "#a2c8ec",
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99"
]

import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl

#test="../AnotherFolder";
test="./"; #Actual folder, not used but in case you want to run it from a specific folder you can use it to implement it


firsttime=True
gmm=1.4
col=0
for test in teststocompare:

    folder     = test[0] #Folder to visit if possible
    namescheme = test[1] #Name for the legend
    ls         = test[2] #Linestyle
    cl         = test[3] #color

    delimiter_in = ' '
    headerlines_in = 1

    if os.path.isdir(folder):  #CONDITION: Is it a folder? If yes go on
        # Importing the MESH DATA
        import_mesh = [file for file in os.listdir(folder) if file.startswith('MESH_')]
        import_mesh.sort()
        nfiles_mesh = len(import_mesh)
        for i in range(nfiles_mesh):
            if i == 0:
                filename = import_mesh[i]
                mydata_mesh = np.loadtxt(folder+"/"+filename, delimiter=delimiter_in, skiprows=headerlines_in)
                x = mydata_mesh
                x_ini = min(x)
                x_end = max(x)
            if i == 1:
                filename = import_mesh[i]
                mydata_mesh = np.loadtxt(folder+"/"+filename, delimiter=delimiter_in, skiprows=headerlines_in)
                y = mydata_mesh
                y_ini = min(y)
                y_end = max(y)

        X, Y = np.meshgrid(x, y)        


        # Importing the SOLUTION DATA
        import_solution = [file for file in os.listdir(folder) if file.startswith('SOLUTION_OCT')]
        import_solution.sort()
        filename=import_solution[-1]


        mydata_solution = np.loadtxt(folder+"/"+ filename, skiprows=headerlines_in)

        ro     = mydata_solution[:, 0]
        rou    = mydata_solution[:, 1]
        rov    = mydata_solution[:, 2]
        energy = mydata_solution[:, 3]
        phi    = mydata_solution[:, 4]

        RO     = np.reshape(ro , (X.shape[0], Y.shape[1]))
        ROU    = np.reshape(rou  , (X.shape[0], Y.shape[1]))
        ROV    = np.reshape(rov  , (X.shape[0], Y.shape[1]))
        ENERGY = np.reshape(energy  , (X.shape[0], Y.shape[1]))
        PHI    = np.reshape(phi, (X.shape[0], Y.shape[1]))


        RO     = np.diag(RO)
        ROU    = np.diag(ROU)
        ROV    = np.diag(ROV)
        ENERGY = np.diag(ENERGY)

        U=ROU/RO
        V=ROV/RO
        modV=np.sqrt(U**2+V**2)
        z=np.diag(X)

        RO   = RO[z > 0]
        ROU    = ROU[z > 0]
        ROV    = ROV[z > 0]
        ENERGY = ENERGY[z > 0]
        modV = modV[z > 0]
        z    = z[z>0]

        U=ROU/RO
        V=ROV/RO

        P=(gmm-1.0)*(ENERGY-0.5*RO*(U**2+V**2))

        #passing to the diagonal
        z    = np.sqrt(2)*z

        if firsttime:
            fig, axs = pl.subplots(1,3, figsize=(12, 4)) #Array of subplots
            # fig.suptitle("")
            firsttime=False

        lw=1.5

        axs[0].plot(z,RO,label=(namescheme),linewidth=lw,color=cl)
        axs[0].set_title(r"$\rho$")
        axs[0].set_xlabel("r")
        axs[0].grid()
        # axs[0].set_ylabel("y")
        # axs[0].legend("H")

        #u
        axs[1].plot(z,modV,label=(namescheme),linewidth=lw,color=cl)
        axs[1].set_title(r"$\sqrt{u^2+v^2}$")
        axs[1].set_xlabel("r")
        axs[1].grid()


        #q
        axs[2].plot(z,P,label=(namescheme),linewidth=lw,color=cl)
        axs[2].set_title("p")
        axs[2].set_xlabel("r")
        axs[2].grid()

 
fig.tight_layout()
for row in axs: 
    row.grid(which='major', color='#666666', linestyle='-',alpha=0.2)
    row.minorticks_on()  
    row.grid(which='minor', color='#666666', linestyle='-',alpha=0.2) 

# Move the legend outside the plot
pl.legend(loc='lower left', fontsize='7', bbox_to_anchor=(1.05, 0.45))

params = {'mathtext.default': 'regular' }   
pl.savefig("plotting_compare_diag.pdf", format="pdf", bbox_inches="tight")
# # pl.savefig("plotting_compare.png",dpi=600)
pl.show()
