# We need a couple of packages in this chapter
import numpy as np  
# This is the basic package in python with all the numerical functions

import matplotlib.pyplot as plt 
# This package allows to  plot

from timeit import timeit
import csv
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl
from soupsieve import select


# variables="primitives" #"all" #"primitives"
# variable_for_convergence="W_X" #"U", "W_X", "W_Y"
# nametest="unsteady_vortex"
# epsilon_vector=[1.0] #[1.0,0.1,0.01]
# scheme_vector=["AF"] #["FV","AF"] #NB: Only for AF
# post_processing_vector=[0] 
# setting_vector=["EX"] #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24] #[23]
# order_time_vector=[2] #[2,3,4,5] #[2,3]
# CFL=0.45

variables="primitives" #"all" #"primitives"
variable_for_convergence="W_Y" #"U", "W_X", "W_Y"
nametest="gresho_vortex"
epsilon_vector=[1.0] #[1.0,0.1,0.01]
scheme_vector=["AF"] #NB: Only for AF
post_processing_vector=[0] 
setting_vector=["IMEX","EX"] #"IMEX","EX"
space_reconstruction_vector=[23] #[20,21,22,23,24] #[23]
order_time_vector=[2] #[2,3,4,5] #[2,3]
CFL=0.45

refabscissae=np.array([25.0,50.0,100.0,200.0,400.0,800.0,1600.0])

####################################################
name_post_processing = {0:"NOPP",1:"PP"}
colors_setting = {"IMEX":"#1f77b4","EX":"#ff7f0e"}
linestyles_order_time        = {2:"-",3:"--",4:"-.",5:":"}
markers_space_reconstruction = {20:"o",21:"^",22:"*",23:"v",24:"s"}
ms=6.0
####################################################

nVar=4


# colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
#    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


#Addressing the proper files
error_file_start=""
if variable_for_convergence=="U":
    error_file_start=""
elif variable_for_convergence=="W_X":
    error_file_start="W_X_"
elif variable_for_convergence=="W_Y":
    error_file_start="W_Y_"
else:
    print("Variable for convergence not available")
    print(variable_for_convergence)
    quit()




#test="../AnotherFolder";
test="./"; #Actual folder, not used but in case you want to run it from a specific folder you can use it to implement it


fig = pl.figure(figsize=(20, 20)) #Array of subplots


for epsilon in epsilon_vector:
    for setting in setting_vector:
        for scheme in scheme_vector:
            for post_processing in post_processing_vector:
                for space_reconstruction in space_reconstruction_vector:
                    for order_time in order_time_vector:
                        foldName="tests_"+variables+"/"+nametest+"/eps"+str(epsilon)+"/"+scheme+"/PP"+str(post_processing)+"/"+setting+"/space_reconstruction_"+str(space_reconstruction)+"/DeC"+str(order_time)+"/CFL"+str(CFL)
                        if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
                            #print(bffolder)
                            count=0
                            errorfiles=[]

                            for file in os.listdir(foldName): #CONDITION: Is there more than 1 error files?
                                if file.startswith(error_file_start+"ErrorL1"):
                                    count=count+1
                                    errorfiles.append(file)
                            if count>1: #If yes go on, MAKE THE CONVERGENCE ANALYSIS
                                errorfiles.sort() #Order the files
                                #Check where you are
                                print("You are in the folder"+foldName+" where you have "+str(count)+" L1 errors")
                                #Defining the errors and mesh parameters vectors
                                numberofcells_X = np.zeros(count)
                                numberofcells_Y = np.zeros(count)
                                errors = np.zeros((count,nVar))

                                #Opening the file to write the rates of convergence
                                fid = open(foldName+"/"+error_file_start+"convergence.tex",'w')
                                towrite="  N"
                                for indj in range(nVar):
                                    towrite+="      L1 error      Order"
                                towrite+=" \\ \n"
                                print(towrite)
                                fid.write(towrite)  
                                for indi in range(count): #Process the files
                                    try:
                                        fileID = open(foldName+"/"+errorfiles[indi],'r')
                                        for indl, line in enumerate(fileID): # iterate over each line
                                            # print(indl,line)
                                            if indl==0:
                                                # Split the line by whitespace
                                                values = line.split()
                                                # Convert the values to float and create a NumPy array
                                                vec = np.array([float(value) for value in values])
                                                numberofcells_X[indi] = int(vec[0]) 
                                                numberofcells_Y[indi] = int(vec[1]) 
                                                errors[indi,:] = vec[2:]
                                        fileID.close()    
                                    except:
                                        errors[indi,:]=float("NaN")
                                    if indi>0:
                                        experimentalorder = np.array([-np.log(errors[indi,j]/errors[indi-1,j])/(np.log(numberofcells_X[indi]/numberofcells_X[indi-1])) for j in range(nVar)])
                                    else:
                                        experimentalorder = np.zeros(nVar)
                                    towrite=""
                                    towrite+=" "+str(int(numberofcells_X[indi]))
                                    for indj in range(nVar):
                                        towrite+="  &   "+format(errors[indi,indj], '.3e')+"  &  "+format(experimentalorder[indj], '.3f')
                                    towrite+=" \\ \n"
                                    print(towrite)
                                    fid.write(towrite)  
                                fid.close()
                                #Plot
                                if setting=="EX":
                                    labeltime="DeC"
                                elif setting=="IMEX":
                                    labeltime="IMEXDeC"
                                else:
                                    print("Error in label time")
                                    quit()
                                pl.loglog(numberofcells_X,errors[:,0],markersize=ms,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_order_time[order_time], linewidth=1.5,label=scheme+" MUSCL"+str(space_reconstruction)+" "+str(labeltime)+str(order_time)+" "+name_post_processing[post_processing],color=colors_setting[setting])


pl.loglog(refabscissae,1e1*refabscissae**(-1),":",linewidth=2,label="order 2",color="0.70")
pl.loglog(refabscissae,2e2*refabscissae**(-2),":",linewidth=2,label="order 2",color="0.50")
pl.loglog(refabscissae,1e3*refabscissae**(-3),":",linewidth=2,label="order 3",color="0.30")
# pl.loglog(refabscissae,1e8*refabscissae**(-5),":",linewidth=2,label="order 5",color="0.50")
# pl.loglog(refabscissae,1e11*refabscissae**(-7),":",linewidth=2,label="order 7",color="0.40")
# pl.loglog(refabscissae[:-1],1e11*refabscissae[:-1]**(-9),":",linewidth=2,label="order 9",color="0.30")
# pl.loglog(refabscissae[:-1],5e13*refabscissae[:-1]**(-11),":",linewidth=2,label="order 11",color="0.20")
# pl.loglog(refabscissae[:-2],5e15*refabscissae[:-2]**(-13),":",linewidth=2,label="order 13",color="0.10")

pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
params = {'mathtext.default': 'regular' }   
pl.xlabel("$N_{x}$")
pl.ylabel("$L^1$ error density")
#fig.suptitle("DeC error convergence")
fig.set_tight_layout(True)
# pl.grid(True,which="both")
pl.grid()
pl.savefig(nametest+"_"+variables+"_"+error_file_start+"convergence_Euler_2D.pdf", format="pdf", bbox_inches="tight")
# pl.savefig(nametest+"convergence_Euler_2D.png",dpi=600)
# pl.show()





