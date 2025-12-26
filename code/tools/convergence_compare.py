teststocompare=[] 

#PP
# teststocompare.append(["AF/PP1/M2SSPRK2", "SAFFVO MUSCL"   ,"-","#1f77b4"]) 
# teststocompare.append(["AF/PP1/M21SSPRK2","SAFFVO k-MUSCL" ,"-","#ff7f0e"])
# teststocompare.append(["AF/PP1/M23SSPRK2","SAFFVO VL-MUSCL","-","#2ca02c"])
# teststocompare.append(["AF/PP1/M24SSPRK2","SAFFVO M-MUSCL" ,"-","#d62728"]) 


#NOPP
# teststocompare.append(["AF/NOPP/M2SSPRK2", "NO-PP SAFFVO MUSCL"   ,"-","#1f77b4"]) 
# teststocompare.append(["AF/NOPP/M21SSPRK2","NO-PP SAFFVO k-MUSCL" ,"-","#ff7f0e"])
# teststocompare.append(["AF/NOPP/M23SSPRK2","NO-PP SAFFVO VL-MUSCL","-","#2ca02c"])
# teststocompare.append(["AF/NOPP/M24SSPRK2","NO-PP SAFFVO M-MUSCL" ,"-","#d62728"]) 


#FV
# teststocompare.append(["rusanov/M2SSPRK2", "FV MUSCL"   ,"--","#1f77b4"]) 
# teststocompare.append(["rusanov/M21SSPRK2","FV k-MUSCL" ,"--","#ff7f0e"])
# teststocompare.append(["rusanov/M23SSPRK2","FV VL-MUSCL","--","#2ca02c"])
# teststocompare.append(["rusanov/M24SSPRK2","FV M-MUSCL" ,"--","#d62728"]) 

#PP
teststocompare.append(["AF/PP1/M2SSPRK2", "SAFFVO MUSCL"   ,"-","#1f77b4"]) 
teststocompare.append(["AF/PP1/M21SSPRK2","SAFFVO k-MUSCL" ,"-","#ff7f0e"])
teststocompare.append(["AF/PP1/M23SSPRK2","SAFFVO VL-MUSCL","-","#2ca02c"])
teststocompare.append(["AF/PP1/M24SSPRK2","SAFFVO M-MUSCL" ,"-","#d62728"]) 


#NOPP
teststocompare.append(["AF/NOPP/M2SSPRK2", "NO-PP SAFFVO MUSCL"   ,"--","#1f77b4"]) 
teststocompare.append(["AF/NOPP/M21SSPRK2","NO-PP SAFFVO k-MUSCL" ,"--","#ff7f0e"])
teststocompare.append(["AF/NOPP/M23SSPRK2","NO-PP SAFFVO VL-MUSCL","--","#2ca02c"])
teststocompare.append(["AF/NOPP/M24SSPRK2","NO-PP SAFFVO M-MUSCL" ,"--","#d62728"]) 


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

col=0
for test in teststocompare:

    folder     = test[0] #Folder to visit if possible
    namescheme = test[1] #Name for the legend
    ls         = test[2] #Linestyle
    cl         = test[3] #color

    if os.path.isdir(folder):  #CONDITION: Is it a folder? If yes go on
        count=0
        errorfiles=[]
        for file in os.listdir(folder): #CONDITION: Is there more than 1 error files?
            if file.startswith("ErrorL1"):
                count=count+1
                errorfiles.append(file)
        if count>1: #If yes go on, MAKE THE CONVERGENCE ANALYSIS
            errorfiles.sort() #Order the files
            #Check where you are
            print("You are in the folder "+folder+" where you have "+str(count)+" L1 errors")
            #Defining the errors and mesh parameters vectors
            numberofelements = np.zeros(count)
            errors = np.zeros((count,4))
            #Opening the file to write the rates of convergence
            fid = open(folder+"/convergence_density.tex",'w')

            print("  N     L1 error ro     Order ro     ")
            fid.write("  N   & L1 error ro  &  Order ro \\ \n")  

            for indi in range(count): #Process the files
                try:
                    fileID = open(folder+"/"+errorfiles[indi],'r')
                    numberofelements[indi] = int(errorfiles[indi][-8:-4])
                    for indl, line in enumerate(fileID): # iterate over each line
                        #print(indl,line)
                        if indl==0:
                            ro, rou, rov, E = line.split() # split it by whitespace
                            errors[indi,0] = float(ro) 
                            errors[indi,1] = float(rou) 
                            errors[indi,2] = float(rov) 
                            errors[indi,3] = float(E) 
                    fileID.close()    
                except:
                    errors[indi,:]=float("NaN")
                if indi>0:
                    order = np.array([-np.log(errors[indi,j]/errors[indi-1,j])/(np.log(numberofelements[indi]/numberofelements[indi-1])) for j in range(4)])
                else:
                    order = np.zeros(2)
                print(int(numberofelements[indi]),"     ",format(errors[indi,0], '.3e'),"     ", format(order[0], '.3f'))
                fid.write(" "+str(int(numberofelements[indi]))+"  &   "+format(errors[indi,0], '.3e')+"  &  "+format(order[0], '.3f')+" \\ \n")  
            fid.close()
            #Plot
            if firsttime:
                fig=pl.figure()
                firsttime=False
            pl.loglog(numberofelements,errors[:,0],ls,linewidth=2,label=(namescheme),color=cl)


pl.grid()
# pl.loglog(numberofelements,errors[0,0]/numberofelements[0]**(-1)*numberofelements**(-1),"--",linewidth=1,label="order 1")
# pl.loglog(numberofelements,errors[0,0]/numberofelements[0]**(-2)*numberofelements**(-2),"--",linewidth=1,label="order 2")
lighter_color = (0.7, 0.7, 0.7)  # RGB values for a lighter gray color
pl.loglog(numberofelements,200*numberofelements[0]**(-2)/numberofelements[0]**(-1)*numberofelements**(-1),":",linewidth=1.8,label="order 1",color=lighter_color)
pl.loglog(numberofelements,200*numberofelements**(-2),":",linewidth=1.8,label="order 2",color="k")

pl.legend(loc='lower left',fontsize='7')
params = {'mathtext.default': 'regular' }   
pl.xlabel("$N$")
pl.ylabel("$L^1$ error")
pl.savefig("convergenceCOMPARE.pdf", format="pdf", bbox_inches="tight")
# pl.savefig("convergenceCOMPARE.png",dpi=600)
pl.show()
