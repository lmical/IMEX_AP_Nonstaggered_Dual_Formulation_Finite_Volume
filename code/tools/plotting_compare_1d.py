teststocompare=[] 

coordinate="X" #"X" #"Y" #NB:Actually y is not implemented

#PP
teststocompare.append(["AF/PP1/M2SSPRK2", "SAFFVO MUSCL"   ,"-","#1f77b4"]) 
# teststocompare.append(["AF/PP1/M21SSPRK2","SAFFVO k-MUSCL" ,"-","#ff7f0e"])
# teststocompare.append(["AF/PP1/M23SSPRK2","SAFFVO VL-MUSCL","-","#2ca02c"])
# teststocompare.append(["AF/PP1/M24SSPRK2","SAFFVO M-MUSCL" ,"-","#d62728"]) 


#NOPP
teststocompare.append(["AF/NOPP/M2SSPRK2", "NO-PP SAFFVO MUSCL"   ,"-","#ff7f0e"]) 
# teststocompare.append(["AF/NOPP/M21SSPRK2","NO-PP SAFFVO k-MUSCL" ,"-","#ff7f0e"])
# teststocompare.append(["AF/NOPP/M23SSPRK2","NO-PP SAFFVO VL-MUSCL","-","#2ca02c"])
# teststocompare.append(["AF/NOPP/M24SSPRK2","NO-PP SAFFVO M-MUSCL" ,"-","#d62728"]) 


#FV
# teststocompare.append(["Rusanov/M2SSPRK2", "FV MUSCL"   ,"-","#ff7f0e"]) 
# teststocompare.append(["rusanov/M21SSPRK2","FV k-MUSCL" ,"--","#ff7f0e"])
# teststocompare.append(["rusanov/M23SSPRK2","FV VL-MUSCL","--","#2ca02c"])
# teststocompare.append(["rusanov/M24SSPRK2","FV M-MUSCL" ,"--","#d62728"]) 

#PP
# teststocompare.append(["AF/PP1/M2SSPRK2", "SAFFVO MUSCL"   ,"-","#1f77b4"]) 
# teststocompare.append(["AF/PP1/M21SSPRK2","SAFFVO k-MUSCL" ,"-","#ff7f0e"])
# teststocompare.append(["AF/PP1/M23SSPRK2","SAFFVO VL-MUSCL","-","#2ca02c"])
# teststocompare.append(["AF/PP1/M24SSPRK2","SAFFVO M-MUSCL" ,"-","#d62728"]) 


#NOPP
# teststocompare.append(["AF/NOPP/M2SSPRK2", "NO-PP SAFFVO MUSCL"   ,"--","#1f77b4"]) 
# teststocompare.append(["AF/NOPP/M21SSPRK2","NO-PP SAFFVO k-MUSCL" ,"--","#ff7f0e"])
# teststocompare.append(["AF/NOPP/M23SSPRK2","NO-PP SAFFVO VL-MUSCL","--","#2ca02c"])
# teststocompare.append(["AF/NOPP/M24SSPRK2","NO-PP SAFFVO M-MUSCL" ,"--","#d62728"]) 


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

    if os.path.isdir(folder):  #CONDITION: Is it a folder? If yes go on
        filename = "MESH_Coordinate"+coordinate+".dat"
        delimiter_in = ' '
        headerlines_in = 1
        mydata_mesh = np.loadtxt(folder+"/"+ filename, delimiter=delimiter_in, skiprows=headerlines_in)
        z = mydata_mesh
        z_ini = min(z)
        z_end = max(z)
        Nz    = len(z)
        
        # Importing the SOLUTION DATA
        import_solution = [file for file in os.listdir(folder) if file.startswith('SOLUTION_OCT')]
        import_solution.sort()
        filename=import_solution[-1]


        mydata_solution = np.loadtxt(folder+"/"+ filename, skiprows=headerlines_in)

        ro   = mydata_solution[:, 0]
        rou    = mydata_solution[:, 1]
        rov    = mydata_solution[:, 2]
        energy    = mydata_solution[:, 3]
        phi  = mydata_solution[:, 4]

        RO     = ro[:Nz]
        ROU    = rou[:Nz]
        ROV    = rov[:Nz]
        ENERGY = energy[:Nz]
        PHI    = phi[:Nz]

        U=ROU/RO
        V=ROV/RO

        P=(gmm-1.0)*(ENERGY-0.5*RO*(U**2+V**2))

        if firsttime:
            fig, axs = pl.subplots(1,3, figsize=(12, 4)) #Array of subplots
            # fig.suptitle("")
            firsttime=False

        lw=1.5

        axs[0].plot(z,RO,label=(namescheme),linewidth=lw,color=cl)
        axs[0].set_title(r"$\rho$")
        axs[0].set_xlabel("x")
        axs[0].grid()
        # axs[0].set_ylabel("y")
        # axs[0].legend("H")

        #u
        axs[1].plot(z,ROU/RO,label=(namescheme),linewidth=lw,color=cl)
        axs[1].set_title("u")
        axs[1].set_xlabel("x")
        axs[1].grid()


        #q
        axs[2].plot(z,P,label=(namescheme),linewidth=lw,color=cl)
        axs[2].set_title("p")
        axs[2].set_xlabel("x")
        axs[2].grid()

 
fig.tight_layout()
for row in axs: 
    row.grid(which='major', color='#666666', linestyle='-',alpha=0.2)
    row.minorticks_on()  
    row.grid(which='minor', color='#666666', linestyle='-',alpha=0.2) 

# Move the legend outside the plot
pl.legend(loc='lower left', fontsize='7', bbox_to_anchor=(1.05, 0.45))

params = {'mathtext.default': 'regular' }   
pl.savefig("plotting_compare.pdf", format="pdf", bbox_inches="tight")
# # pl.savefig("plotting_compare.png",dpi=600)
pl.show()
