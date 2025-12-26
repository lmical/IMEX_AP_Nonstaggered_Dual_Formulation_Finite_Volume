teststocompare=[] 


num=10
whichvariable="primitives"
name="W_X"


# nametest="gresho_vortex"
# nametest="smooth_gresho_vortex"
# nametest="smooth_vortex_eps"
nametest="smooth_vortex_eps_modified"


scheme              ="AF"
setting="IMEX" #"IMEX","EX"
post_processing     =0
CFL                 =0.25
ord=1
if ord==1:
    space_reconstruction=1
    order_time          =1
else:
    space_reconstruction=23
    order_time          =2

teststocompare.append([num,whichvariable,nametest,str(1.0)  ,scheme,post_processing,setting,space_reconstruction,order_time,CFL]) 
teststocompare.append([num,whichvariable,nametest,str(0.1)  ,scheme,post_processing,setting,space_reconstruction,order_time,CFL]) 
teststocompare.append([num,whichvariable,nametest,str(0.01) ,scheme,post_processing,setting,space_reconstruction,order_time,CFL]) 
teststocompare.append([num,whichvariable,nametest,str(0.001),scheme,post_processing,setting,space_reconstruction,order_time,CFL]) 
teststocompare.append([num,whichvariable,nametest,str(0.0001),scheme,post_processing,setting,space_reconstruction,order_time,CFL]) 


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



col=0
for test in teststocompare:

    inum           = test[0] 
    iwhichvariable = test[1] 
    inametest      = test[2]      
    ieps           = test[3] 
    ischeme        = test[4]
    ipp            = test[5]
    isetting       = test[6]
    ispace         = test[7]
    itime          = test[8]
    iCFL           = test[9]


    foldName="FINAL_RESULTS_"+str(inum)+"_"+iwhichvariable+"/"+inametest+"/eps"+str(ieps)+"/"+ischeme+"/PP"+str(ipp)+"/"+isetting+"/space_reconstruction_"+str(ispace)+"/DeC"+str(itime)+"/CFL"+str(iCFL)
            #/FINAL_RESULTS_10_primitives/smooth_gresho_vortex/eps0.1/AF/PP0/IMEX/space_reconstruction_1/DeC1/CFL0.25

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

                            nx, ny, rho, rhou, rhov, E = line.split() # split it by whitespace
                            numberofelements[indi]=nx
                            errors[indi,0] = float(rho) 
                            errors[indi,1] = float(rhou)
                            errors[indi,2] = float(rhov) 
                            errors[indi,3] = float(E) 


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


            ax_ro.loglog(numberofelements,errors[:,0],label="$\epsilon$="+str(ieps))
            ax_u.loglog(numberofelements,errors[:,1] ,label="$\epsilon$="+str(ieps))
            ax_v.loglog(numberofelements,errors[:,2] ,label="$\epsilon$="+str(ieps))
            ax_p.loglog(numberofelements,errors[:,3] ,label="$\epsilon$="+str(ieps))




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
plt.savefig("convergence_"+name+"_"+nametest+"_K"+str(num)+"_order"+str(order_time)+".pdf", format="pdf", bbox_inches="tight")
# pl.savefig("convergenceCOMPARE.png",dpi=600)
plt.show()
