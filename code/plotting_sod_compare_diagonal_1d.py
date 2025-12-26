teststocompare=[] 

nametest="2D_Sod"
solution_file="SOLUTION_OCT_0000000.2500000.dat"
epsilon_test=1.0 #1.0,0.1,0.01
scheme_test="AF" #NB: Only for AF
post_processing_test=1 
# setting_test="EX" #"IMEX","EX"
space_reconstruction_test=23 #[20,21,22,23,24]
CFL_test=0.45


coordinate="X" #"X" #"Y" #NB:Actually y is not implemented


gmm=1.4


teststocompare=[] 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, space_reconstruction_test, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, space_reconstruction_test, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, space_reconstruction_test, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, space_reconstruction_test, 5, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, space_reconstruction_test, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, space_reconstruction_test, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, space_reconstruction_test, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, space_reconstruction_test, 5, CFL_test]) 

# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 20, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 20, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 20, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 20, 5, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 20, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 20, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 20, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 20, 5, CFL_test]) 



# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 21, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 21, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 21, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 21, 5, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 21, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 21, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 21, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 21, 5, CFL_test]) 

# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 22, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 22, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 22, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 22, 5, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 22, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 22, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 22, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 22, 5, CFL_test]) 

# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 23, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 23, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 23, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 23, 5, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 23, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 23, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 23, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 23, 5, CFL_test]) 

# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 24, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 24, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 24, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 24, 5, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 24, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 24, 3, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 24, 4, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 24, 5, CFL_test]) 


# teststocompare.append([epsilon_test, "EX",   scheme_test,  post_processing_test, 23, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 23, 2, CFL_test]) 
# teststocompare.append([epsilon_test, "EX",          "FV",  post_processing_test, 23, 2, CFL_test]) 

teststocompare.append([epsilon_test, "IMEX", scheme_test,  post_processing_test, 23, 2, CFL_test]) 


import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl


####################################################
colors_setting = {"IMEX":"#1f77b4","EX":"#ff7f0e"}
linestyles_order_time        = {2:"-",3:"--",4:"-.",5:":"}
markers_space_reconstruction = {20:"o",21:"^",22:"*",23:"v",24:"s"}
lw=0.1
ms=2
####################################################

import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl

#test="../AnotherFolder";
test="./"; #Actual folder, not used but in case you want to run it from a specific folder you can use it to implement it

fig, axs = pl.subplots(1,3, figsize=(20, 8)) #Array of subplots
# fig.suptitle("")
firsttime=False

col=0
for indt, test in enumerate(teststocompare): #Loop on the schemes
    epsilon              = test[0]
    setting              = test[1]
    scheme               = test[2]
    post_processing      = test[3]
    space_reconstruction = test[4]
    order_time           = test[5]
    CFL                  = test[6]


    foldName="safety_tests_all/"+nametest+"/eps"+str(epsilon)+"/"+scheme+"/PP"+str(post_processing)+"/"+setting+"/space_reconstruction_"+str(space_reconstruction)+"/DeC"+str(order_time)+"/CFL"+str(CFL)


    delimiter_in = ' '
    headerlines_in = 1

    if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
        # Importing the MESH DATA
        import_mesh = [file for file in os.listdir(foldName) if file.startswith('MESH_')]
        import_mesh.sort()
        nfiles_mesh = len(import_mesh)
        for i in range(nfiles_mesh):
            if i == 0:
                filename = import_mesh[i]
                mydata_mesh = np.loadtxt(foldName+"/"+filename, delimiter=delimiter_in, skiprows=headerlines_in)
                x = mydata_mesh
                x_ini = min(x)
                x_end = max(x)
            if i == 1:
                filename = import_mesh[i]
                mydata_mesh = np.loadtxt(foldName+"/"+filename, delimiter=delimiter_in, skiprows=headerlines_in)
                y = mydata_mesh
                y_ini = min(y)
                y_end = max(y)

        X, Y = np.meshgrid(x, y)        


        # Importing the SOLUTION DATA
        filename=solution_file
        if os.path.isfile(foldName+"/"+ filename):

            mydata_solution = np.loadtxt(foldName+"/"+ filename, skiprows=headerlines_in)

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


            RO   = RO[z > 0]
            U    = U[z > 0]
            V    = V[z > 0]
            P    = P[z > 0]
            modV = modV[z > 0]
            z    = z[z>0]

            #passing to the diagonal
            z    = np.sqrt(2)*z


            #Plot
            if setting=="EX":
                labeltime="DeC"
            elif setting=="IMEX":
                labeltime="IMEXDeC"
            else:
                print("Error in label time")
                quit()

            #ro
            if scheme=="FV":
                axs[0].plot(z,RO,markersize=ms,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_order_time[order_time], linewidth=1.5,label=scheme+" MUSCL"+str(space_reconstruction)+" "+str(labeltime)+str(order_time),color="k")            
            else:
                axs[0].plot(z,RO,markersize=ms,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_order_time[order_time], linewidth=1.5,label=scheme+" MUSCL"+str(space_reconstruction)+" "+str(labeltime)+str(order_time),color=colors_setting[setting])
            axs[0].set_title(r"$\rho$")
            axs[0].set_xlabel("r")
            axs[0].grid()


            #u
            if scheme=="FV":
                axs[1].plot(z,modV,markersize=ms,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_order_time[order_time], linewidth=1.5,label=scheme+" MUSCL"+str(space_reconstruction)+" "+str(labeltime)+str(order_time),color="k")
            else:
                axs[1].plot(z,modV,markersize=ms,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_order_time[order_time], linewidth=1.5,label=scheme+" MUSCL"+str(space_reconstruction)+" "+str(labeltime)+str(order_time),color=colors_setting[setting])
            axs[1].set_title(r"$\sqrt{u^2+v^2}$")
            axs[1].set_xlabel("r")
            axs[1].grid()


            #q
            if scheme=="FV":
                axs[2].plot(z,P,markersize=ms,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_order_time[order_time], linewidth=1.5,label=scheme+" MUSCL"+str(space_reconstruction)+" "+str(labeltime)+str(order_time),color="k")
            else:
                axs[2].plot(z,P,markersize=ms,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_order_time[order_time], linewidth=1.5,label=scheme+" MUSCL"+str(space_reconstruction)+" "+str(labeltime)+str(order_time),color=colors_setting[setting])
            axs[2].set_title("p")
            axs[2].set_xlabel("r")
            axs[2].grid()

 
fig.tight_layout()
for row in axs: 
    row.grid(which='major', color='#666666', linestyle='-',alpha=0.2)
    row.minorticks_on()  
    row.grid(which='minor', color='#666666', linestyle='-',alpha=0.2) 

# Move the legend outside the plot
pl.legend(loc='lower left', fontsize='7', bbox_to_anchor=(1.05, 0.05))

params = {'mathtext.default': 'regular' }   
pl.savefig(nametest+"_compare_diag.pdf", format="pdf", bbox_inches="tight")
# # pl.savefig(nametest+"_compare_diag.png",dpi=600)
pl.show()
