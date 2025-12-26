name="WC"


# nametest="gresho_256"
# nametest="smooth_gresho_256"
# nametest="low_mach_vortex_256"
nametest="low_mach_vortex_512"
# nametest="low_mach_vortex_64_512"
# nametest="unsteady_vortex_256"

# nametest="gresho_400"
# nametest="smooth_gresho_400"
# nametest="low_mach_vortex_400"
# nametest="unsteady_vortex_400"

K=0.0 #1.0


# name_base_folder="IMEX_non_replacing_primitive_variables"
# name_base_folder="IMEX_standard"
# name_base_folder="IMEX_standard_hyperbolicity_trick_modified"
# name_base_folder="IMEX_new_K_function_of_epsilon"
# name_base_folder    ="IMEX_new_K4" 
# name_base_folder    ="IMEX_replace_K4" 
# name_base_folder    ="IMEX_try_new_replace_K4" 
# name_base_folder    ="IMEX_last_K4"
# name_base_folder    ="IMEX_256"
# name_base_folder    ="IMEX_FINAL_SIMULATIONS"
name_base_folder= "IMEX_DEFINITIVE"
CFL                 =0.1
ord=2
if ord==1:
    space_reconstruction=1
    order_time          =-1
else:
    space_reconstruction=27
    time_scheme         =-22


teststocompare=[] 
teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1.0)  , CFL]) 
teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1e-1)  , CFL]) 
teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1e-2)  , CFL]) 
teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1e-3)  , CFL]) 
# teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1e-4)  , CFL]) 
# teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1e-5)  , CFL]) 
# teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1e-6)  , CFL]) 
# # #NO #teststocompare.append([name_base_folder,nametest,space_reconstruction,time_scheme,K,str(1e-7)  , CFL]) 


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
import matplotlib.pyplot as plt

#test="../AnotherFolder";
test="./"; #Actual folder, not used but in case you want to run it from a specific folder you can use it to implement it


firsttime=True

fig = plt.figure(1, figsize=(20,6))
ax_ro = fig.add_subplot(141)
ax_u  = fig.add_subplot(142)
ax_v  = fig.add_subplot(143)
ax_p  = fig.add_subplot(144)

ax_ro.clear()
ax_ro.set_title(r"$\rho$")
ax_ro.set_xlabel('N')
ax_ro.set_ylabel('Error')


ax_u.clear()
ax_u.set_title(r"$u$")
ax_u.set_xlabel('N')
ax_u.set_ylabel('Error')

ax_v.clear()
ax_v.set_title(r"$v$")
ax_v.set_xlabel('N')
ax_v.set_ylabel('Error')

ax_p.clear()
ax_p.set_title(r"$p$")
ax_p.set_xlabel('N')
ax_p.set_ylabel('Error')



for test in teststocompare:

    iname_base_folder     = test[0] 
    inametest             = test[1] 
    ispace_reconstruction = test[2]      
    itime_scheme          = test[3] 
    iK_coefficient        = test[4]
    iepsilon              = test[5]
    iCFL                  = test[6]

    foldName=iname_base_folder+"/"+inametest+"/eps"+str(iepsilon)+"/space_reconstruction"+str(ispace_reconstruction)+"/time_scheme_"+str(itime_scheme)+"/K"+str(iK_coefficient)+"/CFL"+str(iCFL)

    print(foldName)

    if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
        count=0
        errorfiles=[]
        for file in os.listdir(foldName): #CONDITION: Is there more than 1 error files?
            if file.startswith(name+"_ErrorL1"):
                count=count+1
                errorfiles.append(file)
        if count>1: #If yes go on, MAKE THE CONVERGENCE ANALYSIS
            errorfiles.sort() #Order the files
            #Check where you are
            print("You are in the folder "+foldName+" where you have "+str(count)+" L1 errors")
            #Defining the errors and mesh parameters vectors
            numberofelements = np.zeros(count)
            errors = np.zeros((count,4))

            #Opening the file to write the rates of convergence
            fid = open(foldName+"/convergence.tex",'w')

            print("   N       L1 error rho      Order rho      L1 Error u     Order u      L1 Error v     Order v         L1 Error p        Order p")
            fid.write("   N   & L1 error rho  &  Order rho & L1 Error u  &  Order u & L1 Error v  &  Order v & L1 Error p  &  Order p\\ \n")  


            for indi in range(count): #Process the files
                try:
                    fileID = open(foldName+"/"+errorfiles[indi],'r')
                    for indl, line in enumerate(fileID): # iterate over each line
                        #print(indl,line)
                        if indl==0:

                            nx, ny, rho, u, v, p = line.split() # split it by whitespace
                            numberofelements[indi]=nx
                            errors[indi,0] = float(rho) 
                            errors[indi,1] = float(u)
                            errors[indi,2] = float(v) 
                            errors[indi,3] = float(p) 


                    fileID.close()    
                except:
                    errors[indi,:]=float("NaN")
                if indi>0:
                    order = np.array([-np.log(errors[indi,j]/errors[indi-1,j])/(np.log(numberofelements[indi]/numberofelements[indi-1])) for j in range(4)])
                else:
                    order = np.zeros(4)
                print('{:4d}'.format(int(numberofelements[indi])),"        ",format(errors[indi,0], '.3e'),"        ", format(order[0], '.3f'),"        ",format(errors[indi,1], '.3e'),"        ", format(order[1], '.3f'),"        ",format(errors[indi,2], '.3e'),"        ", format(order[2], '.3f'),"        ",format(errors[indi,3], '.3e'),"        ", format(order[3], '.3f'))
                fid.write(" "+str(int(numberofelements[indi]))+"  &   "+format(errors[indi,0], '.3e')+"  &  "+format(order[0], '.3f')+"  &  "+format(errors[indi,1], '.3e')+" & "+format(order[1], '.3f')+"  &  "+format(errors[indi,2], '.3e')+" & "+format(order[2], '.3f')+"  &  "+format(errors[indi,3], '.3e')+" & "+format(order[3], '.3f')+" \\\\ \n")  
            fid.close()


            ax_ro.loglog(numberofelements,errors[:,0],label="space"+str(ispace_reconstruction)+" time"+str(itime_scheme)+" eps"+str(iepsilon)+" K"+str(iK_coefficient))
            ax_u.loglog(numberofelements,errors[:,1] ,label="space"+str(ispace_reconstruction)+" time"+str(itime_scheme)+" eps"+str(iepsilon)+" K"+str(iK_coefficient))
            ax_v.loglog(numberofelements,errors[:,2] ,label="space"+str(ispace_reconstruction)+" time"+str(itime_scheme)+" eps"+str(iepsilon)+" K"+str(iK_coefficient))
            ax_p.loglog(numberofelements,errors[:,3] ,label="space"+str(ispace_reconstruction)+" time"+str(itime_scheme)+" eps"+str(iepsilon)+" K"+str(iK_coefficient))



ax_ro.grid(True)
ax_u.grid(True)
ax_v.grid(True)
ax_p.grid(True)

ax_ro.loglog(numberofelements,numberofelements**(-1),"--",linewidth=1,label="order 1")
ax_ro.loglog(numberofelements,1e1*numberofelements**(-2),"--",linewidth=1,label="order 2")

ax_u.loglog(numberofelements,numberofelements**(-1),"--",linewidth=1,label="order 1")
ax_u.loglog(numberofelements,1e1*numberofelements**(-2),"--",linewidth=1,label="order 2")

ax_v.loglog(numberofelements,numberofelements**(-1),"--",linewidth=1,label="order 1")
ax_v.loglog(numberofelements,1e1*numberofelements**(-2),"--",linewidth=1,label="order 2")

ax_p.loglog(numberofelements,numberofelements**(-1),"--",linewidth=1,label="order 1")
ax_p.loglog(numberofelements,1e1*numberofelements**(-2),"--",linewidth=1,label="order 2")

lighter_color = (0.7, 0.7, 0.7)  # RGB values for a lighter gray color
# pl.loglog(numberofelements,200*numberofelements[0]**(-2)/numberofelements[0]**(-1)*numberofelements**(-1),":",linewidth=1.8,label="order 1",color=lighter_color)
# pl.loglog(numberofelements,200*numberofelements**(-2),":",linewidth=1.8,label="order 2",color="k")

ax_ro.legend(loc='lower left',fontsize='7')
params = {'mathtext.default': 'regular' }   
plt.savefig("convergence_same_K_"+name_base_folder+"_"+nametest+"_K"+str(K)+"_order"+str(ord)+".pdf", format="pdf", bbox_inches="tight")
# pl.savefig("convergenceCOMPARE.png",dpi=600)
# plt.show()
