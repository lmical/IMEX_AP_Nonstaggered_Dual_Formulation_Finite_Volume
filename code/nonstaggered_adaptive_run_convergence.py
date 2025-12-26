import numpy as np
import shutil, os
import re
from joblib import Parallel, delayed
import multiprocessing
import sys


# numbertest=3
# # nametest="unsteady_vortex_velocity"
# nametest="unsteady_vortex_momentum_imex"
# space_reconstructions =[-26] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-102]
# post_processings     =[5]
# epsilons             =[1.0]
# N_refinements = {1:6,2:6,20:6,21:6,22:6,23:6,24:6,25:6,26:6,-1:6,-2:6,-20:6,-21:6,-22:6,-23:6,-24:6,-25:6,-26:6}
# starting_elements_X=50
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3

# numbertest=101
# nametest="smooth_gresho"
# space_reconstructions =[-26] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[-102]
# post_processings     =[5]
# epsilons             =[0.1]
# N_refinements = {1:6,2:6,20:6,21:6,22:6,23:6,24:6,25:6,26:6,-1:6,-2:6,-20:6,-21:6,-22:6,-23:6,-24:6,-25:6,-26:6}
# starting_elements_X=50
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3


# numbertest=5
# # nametest="unsteady_vortex_momentum"
# nametest="advection_density_momentum"
# space_reconstructions =[-1] #-np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[1]
# K_coefficients       =[0.025]
# N_refinements = {1:6,2:6,20:6,21:6,22:6,23:6,24:6,25:6,26:6,-1:6,-2:6,-20:6,-21:6,-22:6,-23:6,-24:6,-25:6,-26:6}
# starting_elements_X=50
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3


# numbertest=150
# nametest="advection_sin4"
# orders=[3,5,7,9,11,13]
# meshes = {3:9,5:7,7:6,9:5,11:4,13:4}
# starting_elements=20
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =0
# speed_estimate=0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0

# numbertest=151
# nametest="advection_sin4_longer_domain"
# orders=[3,5,7,9,11,13]
# meshes = {3:6,5:6,7:5,9:4,11:3,13:3}
# starting_elements_dict = {3:160,5:80,7:80,9:40,11:40,13:40}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =0
# speed_estimate=0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0

# numbertest=250
# nametest="composite_wave_200"
# orders=[3,5,7]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=200
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=1
# riemann_solver        =-1
# speed_estimate        =0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


# numbertest=251
# starting_elements=200
# nametest="composite_wave_longer_time_"+str(starting_elements)
# orders=[3,5,7]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=1
# riemann_solver        =-2
# speed_estimate        =0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0

# numbertest=151
# nametest="advection_sin4_longer_domain_HO"
# orders=[3,5,7,9,11,13,15,17]#[3,5,7,9,11,13,15,17,19,21,23,25,27,29,31] #[3,5,7,9,11,13]
# meshes = {3:8,5:8,7:6,9:4,11:3,13:3,15:3,17:3,19:3,21:3,23:3,25:3,27:3,29:3,31:3}
# starting_elements=40
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =0
# speed_estimate=0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0

# numbertest=151
# nametest="advection_sin4_longer_domain_communication"
# orders=[3,5,7] #[3,5,7,9,11,13]
# meshes = {3:6,5:6,7:5,9:4,11:3,13:3}
# starting_elements_dict = {3:160,5:80,7:80,9:40,11:40,13:40}  # use elements of 'orders' as keys
# reconstructed_variable=1
# riemann_solver        =6
# speed_estimate=0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


#####################
#CFL 0.8
#100 elements
#####################
#Let us assume standard speed estimate.
#The only simulation crashing is cons with exact riemann solver for order 13.
#All the others manage to arrive to the final time.
#It is possible to make also this setting work if using the exact riemann solver speed estimate.
#In order to run with standard speed estimate, for order 13, with such setting, the CFL must be decreased to 0.8.
#We perform the comparison for CFL=0.8
#The results with reconstruction of conserved variables present some slight oscillations in the central plateau of velocity and pressure.
#Results with rusanov and exact rieman solver are comparable. 
#Rusanov is a bit more diffusive and prevents over- and undershoots when reconstructing characteristic variables 
#Slight over and undershoots are present for the exact riemann solver even when characteristic variables are reconstructed
#Overall the best setting seems to be characteristic variables with Rusanov
#####################
# CFL 0.95
#####################
# Comparison numerical fluxes for reconstruction of characteristic variables
# CFL 0.95 works up to order 7 for all fluxes
#####################
# numbertest=8
# nametest="sod"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95


#####################
#CFL 0.95
#200 elements
#####################
#All simulations run without any problems, for all orders with all settings.
#And differences cannot be appreciated at the macroscopic level.
#The substantial differences are in the density and in the momentum, especially in the first one.
#The zoom show that reconstructing characteristic variables is better.
#For very high orders, the reconstruction of conserved variables affect a bit the scheme and peaks are not well captured. WENO9 outperforms higher order.
#The reconstruction of characteristic variables produces instead much cleaner results.
#The exact Riemann solver produces better results than rusanov and this is particularly true when conserved variables are reconstructed.
#For reconstruction of characteristic variables only slight improvements can be appreciated when the exact Riemann solver is adopted.
#Overall the best setting seems to be characteristic variables with exact Riemann solver.
#The results obtained with high order are always better besides the mentioned case: reconstruction of conserved variables.
#####################
# numbertest=-300
# nametest="shock_turbulence_interaction_shu_osher"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=400
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95


#####################
#CFL 0.95
#1000 elements
#####################
# numbertest=-301
# nametest="shock_turbulence_interaction_titarev_toro"
# space_reconstructions = -np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[12]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1}
# starting_elements_X=1000
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.95

#######################
#NOT WORKING
#CFL=0.1 is not enough
#NB: I think it is a problem of reconstruction getting negative and not of CFL, the simulation does not crash at the beginning.
#It does not even work for reconstruction of char with rusanov and standard speed estimate
#######################
# numbertest=350
# starting_elements_X=800
# nametest="woodward_colella_"+str(starting_elements_X)
# space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29 #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[0.2] #[0.15,0.2,0.25,0.3] #[0.01,0.05,0.1,0.5,1.0]
# WhichVariableForLimiting_PCSD=2
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95

#WRONG BCs
# numbertest=-302
# starting_elements_X=32000
# nametest="shock_turbulence_modified_example_1_"+str(starting_elements_X)
# space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29  #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[2.0] #[0.02,0.025,0.015] #[0] #[0.01,0.05,0.1,0.5,1.0] #1
# WhichVariableForLimiting_PCSD =2 #1,2,3
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95

#WRONG BCs
# numbertest=-303
# starting_elements_X=1200 #To run
# nametest="shock_turbulence_modified_example_2_"+str(starting_elements_X)
# space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29 #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[0.015] #[0.015,0.02,0.025] #[0.0] #[0.01,0.05,0.1,0.5,1.0]
# WhichVariableForLimiting_PCSD =2 #1,2,3
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95

#CORRECT BCs
# numbertest=-304
# starting_elements_X=800 #12800
# nametest="shock_turbulence_modified_example_1_free_BC_"+str(starting_elements_X)
# space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29  #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[1.0] #[1.4] #[1.4] #,1.3,1.4,1.5,1.6,1.7,1.8,1.9] #[2.0] #[2.0] #[2.0,3.0,4.0,5.0] #[0.02,0.025,0.015] #[0] #[0.01,0.05,0.1,0.5,1.0] #1
# WhichVariableForLimiting_PCSD =2 #1,2,3
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95


# numbertest=151
# starting_elements_X=1600
# nametest="sin4_1D"
# space_reconstructions     = [29] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29  #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[0.0] #[2.0] #[2.0,3.0,4.0,5.0] #[0.02,0.025,0.015] #[0] #[0.01,0.05,0.1,0.5,1.0] #1
# WhichVariableForLimiting_PCSD =2 #1,2,3
# N_refinements = {28:5,29:5}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95

#CORRECT BCs
# numbertest=-305
# starting_elements_X=1200 #To run
# nametest="shock_turbulence_modified_example_2_free_BC_"+str(starting_elements_X)
# space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29 #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[0.015] #[0.015,0.02,0.025] #[0.0] #[0.01,0.05,0.1,0.5,1.0]
# WhichVariableForLimiting_PCSD =2 #1,2,3
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95

# numbertest=10
# starting_elements_X=250
# nametest="implosion_problem_"+str(starting_elements_X)
# space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29 #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[25] #20,30,40 #[0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0] #[0.0] #[0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0] #K=25 BEST
# WhichVariableForLimiting_PCSD =2 #1,2,3
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95



# numbertest=11
# starting_elements_X=1000
# nametest="2d_riemann_problem_"+str(starting_elements_X)
# space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# space_reconstructions_fix =  29 #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# K_coefficients       =[30] #20,30,40 #[0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0] #[0.0] #[0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0] #K=25 BEST
# WhichVariableForLimiting_PCSD =2 #1,2,3
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=starting_elements_X
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.4 #0.95

numbertest=12
starting_elements_X=256
nametest="RT_instability_"+str(starting_elements_X)
space_reconstructions     = [28] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
space_reconstructions_fix =  29 #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
time_schemes         =[3]
K_coefficients       =[50.0] #20,30,40 #[0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0] #[0.0] #[0.01,0.05,0.1,0.5,1.0,5.0,10.0,50.0] #K=25 BEST
WhichVariableForLimiting_PCSD =2 #1,2,3
N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1,-26:1,-28:1,-29:1}
starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
starting_elements_Y=1024
starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
speed_estimate=0
CFL=0.4 #0.95


#####################
#CFL 0.95
#100 elements
#####################
#We did not encounter any problem in running this test.
#Results with reconstruction of characteristic variables are much cleaner.
#When reconstructing conserved variables, one gets little spurious oscillations in the plateaux as for Sod.
#Results with exact Riemann solver and Rusanov are similar.
#Very high order methods, from order 9 on, produce more pronounced, yet small, over and undershoots.
#The higher diffusion of Rusanov, compared to the exact Riemann solver, is able to kill them.
#Overall the nicest results are the one with reconstruction of characteristic variables and Rusanov numerical flux.
#Higher order methods seem to better capture the solution details.
#####################
# numbertest=401
# nametest="RP1"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95



#######################
#NOT WORKING
#CFL=0.1 is not enough
#NB: I think it is a problem of reconstruction getting negative and not of CFL
#######################
#From CFL smaller or equal than 0.85, the test works for order 3 if characteristic variables are reconstructed and Rusanov numerical flux is adopted.
#For all other settings and orders, even for order 3, even for speed estimate with exact Riemann solver, even for smaller CFL up to 0.05, the simulations crash.
#Again, we remark that, even though this is out from the goal of this investigation, the approach may benefit from extra asaptive strategies such sa a posteriori limiting.
#This is why we focused on a relaxed version of the problem.
#######################
# numbertest=402
# nametest="RP2"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=2
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95


#####################
#CFL 0.4
#100 elements
#####################
#Let us focus on reconstruction of characteristic variables.
#For CFL=0.95 only order 3 works with both numerical fluxes even for standard speed estimate.
#There are no simulation crashes for both numerical fluxes for CFL smaller or equal to 0.7 even for standard speed estimate.
#Actually, for CFL=0.75, Rusanov numerical flux works for all orders, while, the exact Riemann solver only up to order 11.
#With reconstruction of conserved variables variables, the CFL must be decreased.
#In such a case:
#with exact Riemann solver, from CFL=0.7 and CFL=0.6 only order 3 works even with Riemann speed estimate, from CFL=0.5 to CFL=0.05 only orders 3 and 5 work, 
#with Rusanov, from CFL=0.7 orders 3, 5 and 7 work, for CFL=0.6 orders 3,5,7 and 9 work, for CFL=0.5 all orders work beside order 13, for CFL=0.4 all orders work
#We comment the results for CFL=0.4.
#The reconstruction of characteristic variables outperforms the one one on conserved variables.
#This can be inferred from the plotted profiles for Rusanov numerical flux, while, for exact Riemann solver this is also clear from the fact that the simulation crases if conserved variables are reconstructed from order 7 on.
#Again, the diffusive character of Rusanov is able to clip some spurious over and undershoots.
#Overall, the best results are obtained for Rusanov with reconstruction of characteristic variables and Rusanov.
#Also in this case, the higher is the order, the closer are the results to the reference exact solution.
#####################
#CFL=0.8
#####################
# Comparison numerical fluxes for reconstruction of characteristic variables
#-1, 6
# CFL=0.95 fine for all
# 0,2,3,5
# CFL=0.95 fine for 3; CFL=0.85 fine up to 5; CFL=0.8 up to order 7
# 1
# CFL=0.95 fine for 3; CFL=0.8 fine up to 7.
# 4
# Does not work up to CFL=0.05 for all orders
#####################
# numbertest=412
# nametest="RP2_relaxed"
# space_reconstructions = -np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[12]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.95 !*I don't know


#####################
#CFL 0.95 up to order 9
#CFL 0.3 up to order 11
#100 elements
#####################
#Assuming standard speed estimate,
#for CFL=0.95 all the simulations work up to order 9, besides when characteristic variables are reconstructed and Rusanov is used. 
#In such a case all simulations crash. This is due to the underestimate of the local wave speed, when Riemann speed estimate is assumed, also this setting works up to order 9.
#Simulations work up to order 11 for CFL smaller or equal to 0.45 also with standard speed estimate for reconstruction of characteristic variables and exact Riemann solver.
#However, they only work up to order 9 for all the other configurations even for Riemann speed estimate. This also holds up to CFL=0.05, we did not test for smaller CFLs.
#We perform the comparison for CFL=0.45
#As one would expect from the previous description, the best setting is the reconstruction of characteristic variables in combination with exact Riemann solver.
#This is in fact also the only setting for which the simulation with order 11 does not crash.
#The reconstruction of conserved variables, as usual, causes little spurious oscillations.
#Focusing on the results obtained with reconstruction of characteristic variables and exact Riemann solver, we see that higher order methods better capture many solution details.
#######################
#CFL=0.7
#######################
# Comparison numerical fluxes for reconstruction of characteristic variables
#-1, 1, 2, 4, 5, 6
# CFL=0.95 fine for all
# 0
# CFL=0.8 up to order 5; 0.7 up to 7
# 3
# CFL=0.8 up to order 5; 0.75 up to 7
#######################
# numbertest=403
# nametest="RP3"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95

#####################
#CFL 0.95
#100 elements
#####################
#All settings experience no crashes for all orders for CFL=0.95, besides reconstruction of conserved variables with Rusanov which only works up to order 11 even with Riemann speed estimate.
#We could not manage to find a CFL for which order 13 was running for such a setting, we tested up to CFL=0.05 with the safe Riemann speed estimate.
#We perform the comparison for CFL=0.95.
#In this case, it is particularly evident how the reconstruction of characteristic variables is essential to reduce spurious oscillations.
#Reconstruction of conserved variables produces several ugly oscillations and terrible over- and undershoots which are much smaller if characteristic variables are reconstructed.
#Results with exact Riemann solver are slightly better than the ones obtained with Rusanov.
#Overall, the features of the solution are captured by all orders, however, let us notice that in this case, spurious oscillations are present in plateux and the higher is the order of the method, the higher these oscillations are.
#####################
# numbertest=404
# nametest="RP4"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95


#######################
#CFL=0.95 is not enough
#CFL=0.5 only works for order 13
#CFL=0.3 works from order 7 to 11
#CFL=0.1 works up to order 11
#100 elements
#######################
# The behavior of the methods on this test was quite irregular.
# Reconstruction of conserved variables with Rusanov always crashes for all orders and any value of $C_{CFL}:=0.95$ between 0.95 and 0.05.
# For $C_{CFL}:=0.95$, only reconstruction of characteristic variables with Rusanov for order 3 works. All other settings crash.
# For $C_{CFL}:=0.75$, reconstruction of characteristic variables with Rusanov works up to order 7. All other settings crash.
# For $C_{CFL}\leq 0.65$, reconstruction of characteristic variables with Rusanov works up to order 11. All other setting crash, besides some irregular exceptions involving exact Riemann solver for orders 9, 11 and 13 for some values of $C_{CFL}$, for both variable reconstructions.
# For $C_{CFL}:=0.3$, reconstruction of conserved variables with exact Riemann solver works for order higher than 5,
# while, reconstruction of characteristic variables with exact Riemann solver works from order 7 to 11.
# For $C_{CFL}:=0.2$, reconstruction of characteristic variables with exact Riemann solver works from order 5 to 11.
# For the same $C_{CFL}$, reconstruction of conserved variables works for all orders with exact Riemann solver.
# For $C_{CFL}:=0.1$, all settings work but reconstruction of characteristic variables with Rusanov for order 13 and reconstruction of conserved variables with Rusanov, which always crashes for all orders, as already anticipated.
# The results are analogous if Riemann speed estimate is adopted, i.e., problems in this case do not occur at the beginning of the simulations but rather later on due to the reconstructions producing negative values of density or pressure. As already stated, fixing this issue is not in the goals of this paper and is left for future works.
# Let us comment the results obtained for $C_{CFL}:=0.1$, reported in Figure~\ref{fig:Euler_1d_RP5_all_settings_primitive_variables}.
# As can be inferred from the previous discussion, reconstructing characteristic variables is mandatory here when increasing the order of accuracy.
# For reconstruction of conserved variables, no results are available with Rusanov. With exact Riemann solver, they are available for all orders but present many spurious oscillations, especially in the plateaux of density, velocity and pressure.
# Results with reconstruction of characteristic variables are much better.
# In such a case, as expected, results with Rusanov are less oscillatory, even though order 13 works only for exact Riemann solver.
# However, despite this, we can appreciate a much higher resolution when the exact Riemann solver is employed.
# The density peak is smeared by Rusanov, but is perfectly captured even by order 3 with reconstruction of conserved variables, when exact Riemann solver is employed.
# This can be well appreciated in Figure~\ref{fig:Euler_1d_RP5_zoom_density}, where we report some zooms on the density obtained with reconstruction of characteristic variables and both numerical fluxes.
# Results get better as the order increase, excluded order 13 which is characterized by ugly spurious overshoots not present for other orders.
#######################
#CFL=0.1
#######################
# Comparison numerical fluxes for reconstruction of characteristic variables
#-1
# CFL=0.35 up to order 3; CFL=0.3 up to order 5; CFL=0.25 up to order 7
# 0
# CFL=0.95 up to order 3; 0.85 up to 5; 0.75 up to 7
# 1
# CFL=0.15 up to order 7; 
# NB: 5 and 7 allow for higher CFLs
# order 5 runs for CFL=0.2; order 7 runs for CFL=0.3
# 2
# CFL=0.15 up to order 7; 
# NB: 5 and 7 allow for higher CFLs
# order 5 runs for CFL=0.2; order 7 runs for CFL=0.25
# 3
# CFL=0.4 up to order 7;
# NB: 5 and 7 allow for higher CFLs
# order 5 runs for CFL=0.6; order 7 runs for CFL=0.55
# 4
# CFL=0.45 up to order 7;
# NB: 5 and 7 allow for higher CFLs
# order 5 runs for CFL=0.5; order 7 runs for CFL=0.8
# 5
# CFL=0.1 up to order 7; 
# NB: 5 and 7 allow for higher CFLs
# order 5 runs for CFL=0.15; order 7 runs for CFL=0.15
# 6
# CFL=0.95 ok for all
#######################
# numbertest=405
# nametest="RP5"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=2
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95


#######################
#CFL=0.75
#######################
#For $C_{CFL}\leq 0.95$, exact Riemann solver and any choice of reconstructed variables, no simulation crashes occur.
#For $C_{CFL}:= 0.95$, Rusanov works up to order 5 for both reconstructed variables.
#Even decreasing $C_{CFL}$ up to 0.05, nothing changes for Rusanov with reconstruction of conserved variables.
#This is not the case for Rusanov with reconstruction of characteristic variables:
#for $C_{CFL}:= 0.85$, the related simulations work up to order 9;
#for $C_{CFL}:= 0.8$, they work up to order 11;
#for $C_{CFL}\leq 0.75$, they work for all orders.
#Nothing changes when Riemann speed estimate is adopted.
#We perform the comparison for CFL=0.75
#The reconstruction of characteristic variables results essential to reduce spurious oscillations.
#In fact, the results obtained with reconstruction of conserved variables and exact Rieman solver are quite oscillatory, while, the simulations crash from order 7 on if Rusanov is adopted.
#For reconstruction of characteristic variables all orders manage to arrive to the final time for both Riemann solver.
#Oscillations are not completely absent and they are mainly concentrated in the area of impact between the two streams.
#Overall, the best results are the ones obtained with Rusanov and exact Riemann solver.
#Again, we remark that a posteriori limiting or oder adaptive strategies would be very useful to decrease the oscillations.
#######################
#CFL=0.85
#######################
# Comparison numerical fluxes for reconstruction of characteristic variables
# 0
# CFL=0.95 up to order 5; CFL=0.85 up to order 7
# -1,1,2,3,4,5, 6
# CFL=0.95 ok for all
#######################
# numbertest=406
# nametest="RP6"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95


#######################
#CFL=0.95
#######################
# numbertest=420
# nametest="stationary_contact"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95


#######################
#CFL=0.95
#######################
# numbertest=421
# nametest="moving_contact"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=3
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95

#####################
#CFL 0.95
#100 elements
#####################
# numbertest=500
# nametest="lax"
# space_reconstructions = [-26] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.3 #0.95

#####################
#CFL 0.95
#100 elements
#####################
# numbertest=501
# nametest="lax_smaller_lambda"
# space_reconstructions = np.array([-1, -20, -21, -22, -23, -24]) #put sign -
# time_schemes         =[12]
# post_processings     =[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.95 



####################################################
name_reconstructed_variable = {0:"cons",    1:"char"}
name_riemann_solver         = {-2:"central", -1:"LF", 0:"rusanov", 1:"exact", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:"FORCEalpha"}
name_speed_estimate         = {0:"standard",1:"riemann"}
####################################################

if WhichVariableForLimiting_PCSD==1:
    name_base_folder="conserved_latest_adaptivity_based_on_rho"
elif WhichVariableForLimiting_PCSD==2:
    name_base_folder="conserved_latest_adaptivity_based_on_rhou"
elif WhichVariableForLimiting_PCSD==3:
    name_base_folder="conserved_latest_adaptivity_based_on_energy"
else:
    print("Wrong choice for WhichVariableForLimiting_PCSD")
    quit()


# #Clean, Compile
# instruction="make clean; make"
# os.system(instruction)

#Create simulations folders and copy main there
for space_reconstruction in space_reconstructions:
    for time_scheme in time_schemes:
        for K_coefficient in K_coefficients:
                foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/space_reconstruction_fix"+str(space_reconstructions_fix)+"/time_scheme_"+str(time_scheme)+"/K"+str(K_coefficient)+"/CFL"+str(CFL)
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
                foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/space_reconstruction_fix"+str(space_reconstructions_fix)+"/time_scheme_"+str(time_scheme)+"/K"+str(K_coefficient)+"/CFL"+str(CFL)

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
                        instruction+="./main "+str(numbertest)+" "+str(int(element_X))+" "+str(int(elements_Y[indi]))+" "+str(space_reconstruction)+" "+str(space_reconstructions_fix)+" "+str(time_scheme)+" "+str(CFL)+" "+str(K_coefficient)+" "+str(WhichVariableForLimiting_PCSD)+" "+" \n" 
                    elif test_type=="1D":
                        instruction+="./main "+str(numbertest)+" "+str(int(element_X))+" "+str(int(elements_Y[0]))+" "+str(space_reconstruction)+" "+str(space_reconstructions_fix)+" "+str(time_scheme)+" "+str(CFL)+" "+str(K_coefficient)+" "+str(WhichVariableForLimiting_PCSD)+" "+" \n" 
                    else:
                        print("Test type not recognized")
                        quit()


                print(instruction)
                os.system(instruction)
