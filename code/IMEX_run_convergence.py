import numpy as np
import shutil, os
import re
from joblib import Parallel, delayed
import multiprocessing
import sys



# numbertest=200
# nametest="low_mach_vortex_64_512"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes          =[-82]
# epsilons             =[1e-0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6] #[1e-0,1e-1,1e-2,1e-3,1e-4] #[1.0]
# # K_coefficients       =[0.0,1.0,2.0,5.0,10.0]
# K_coefficients       =[0.0]
# N_refinements = {1:4,2:4,20:4,21:4,22:4,23:4,24:4,25:4,26:4,27:4,-1:4,-2:4,-20:4,-21:4,-22:4,-23:4,-24:4,-25:4,-26:4,-27:4}
# starting_elements_X=64
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.2


numbertest=-99
starting_elements_X=128
nametest="gresho_longer_time_gamma1.4_"+str(starting_elements_X)
space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
time_schemes          =[-82]
epsilons             =[1e-6]#,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]#,1e-3,1e-4,1e-5,1e-6] #[1.0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6] #[1.0]
K_coefficients       =[0.0]
N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
starting_elements_Y=starting_elements_X
starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
speed_estimate=0
CFL=0.475 #0.2


# numbertest=301
# nametest="baroclinic_vorticity_generation"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-52]
# epsilons             =[0.05] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=800
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=160
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475


# numbertest=300
# starting_elements_X=256 #200
# nametest="double_shear_layer_"+str(starting_elements_X)
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes          =[-52]
# epsilons             = [1e-0]#[1.0,1e-1,1e-2,1e-3,1e-4] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.475 #0.1




# numbertest=8
# nametest="sod1.0"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-52]
# epsilons             =[1.0] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.475 #0.2

# numbertest=88
# nametest="sod1.0_with_timestep_relaxation"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-52]
# epsilons             =[1.0] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.475 #0.2

# numbertest=-8 
# eps=0.3 #0.3 (works <=0.05) #0.6 (works <=0.25) #0.9
# nametest="sod"+str(eps)
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-82]
# epsilons             =[eps] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.2 #0.475 #0.2


# numbertest=-80
# eps=0.3 #0.3 #0.6 #0.9
# nametest="sod"+str(eps)+"_with_timestep_relaxation"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-52]
# epsilons             =[eps] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.475 #0.2






# numbertest=7
# starting_elements_X=400 #200
# nametest="sod2d1.0"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes          =[-52]
# epsilons             = [1.0]#[1.0,1e-1,1e-2,1e-3,1e-4] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.475 #0.2

# numbertest=77
# nametest="sod2d1.0_with_timestep_relaxation"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-52]
# epsilons             =[1.0] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.2



# numbertest=-7
# eps=0.3 #0.3 (works <=0.05) #0.6 (works <=0.2) #0.9
# nametest="sod2d"+str(eps)
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes          =[-82]
# epsilons             = [eps]#[1.0,1e-1,1e-2,1e-3,1e-4] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=400 #200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.1 #0.475 #0.2


# numbertest=-70
# eps=0.3 #0.3 #0.6 #0.9
# nametest="sod2d"+str(eps)+"_with_timestep_relaxation"
# space_reconstructions =[27] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes          =[-52,-82]
# epsilons             = [eps]#[1.0,1e-1,1e-2,1e-3,1e-4] #[1.0]
# K_coefficients       =[0.0]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=400 #200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.475 #0.2



###########################################3
###########################################3
###########################################3






####################################################
name_reconstructed_variable = {0:"cons",    1:"char"}
name_riemann_solver         = {-2:"central", -1:"LF", 0:"rusanov", 1:"exact", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:"FORCEalpha"}
name_speed_estimate         = {0:"standard",1:"riemann"}
####################################################


name_base_folder    ="IMEX_HYPERBOLIC_TRICK_RHO_P_STAR_STAGE_DEPENDENT" 


# #Clean, Compile
# instruction="make clean; make"
# os.system(instruction)

#Create simulations folders and copy main there
for space_reconstruction in space_reconstructions:
    for time_scheme in time_schemes:
        for K_coefficient in K_coefficients:
            for epsilon in epsilons:
                foldName=name_base_folder+"/"+nametest+"/eps"+str(epsilon)+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/K"+str(K_coefficient)+"/CFL"+str(CFL)
                instruction =""
                instruction+="mkdir -p "+foldName+" \n" #creation of the folder
                os.system(instruction)
                instruction =""
                instruction+="cp bin/main "+foldName+"/main \n" #copying DATA/don1d
                os.system(instruction)

#Run simulations
for space_reconstruction in space_reconstructions:
    for time_scheme in time_schemes:
        for K_coefficient in K_coefficients:
            for epsilon in epsilons:
                foldName=name_base_folder+"/"+nametest+"/eps"+str(epsilon)+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/K"+str(K_coefficient)+"/CFL"+str(CFL)

                elements_X=np.zeros(N_refinements[space_reconstruction])
                elements_X[0]=starting_elements_X_dict[space_reconstruction]
                for ind in range(1,N_refinements[space_reconstruction]):
                    elements_X[ind]=elements_X[ind-1]*2

                elements_Y=np.zeros(N_refinements[space_reconstruction])
                elements_Y[0]=starting_elements_Y_dict[space_reconstruction]
                for ind in range(1,N_refinements[space_reconstruction]):
                    elements_Y[ind]=elements_Y[ind-1]*2

                instruction =""
                instruction+="cd "+foldName+" \n" #move in the folder
                for indi,element_X in enumerate(elements_X):        
                    if test_type=="2D":
                        instruction+="./main "+str(numbertest)+" "+str(int(element_X))+" "+str(int(elements_Y[indi]))+" "+str(space_reconstruction)+" "+str(space_reconstruction)+" "+str(time_scheme)+" "+str(CFL)+" "+str(epsilon)+" "+str(K_coefficient)+" "+" \n" 
                    elif test_type=="1D":
                        instruction+="./main "+str(numbertest)+" "+str(int(element_X))+" "+str(int(elements_Y[0]))+" "+str(space_reconstruction)+" "+str(space_reconstruction)+" "+str(time_scheme)+" "+str(CFL)+" "+str(epsilon)+" "+str(K_coefficient)+" "+" \n" 
                    else:
                        print("Test type not recognized")
                        quit()


                print(instruction)
                os.system(instruction)
