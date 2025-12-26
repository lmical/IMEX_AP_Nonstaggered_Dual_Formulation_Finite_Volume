import numpy as np
import shutil, os
import re
# from joblib import Parallel, delayed
import multiprocessing
import sys

NX=0
NY=0

variables="all" #"all" #"primitives"
# folder="FINAL_RESULTS_10_CONSERVED_IMPLICIT_CRANC_NICOLSON_1e-12"
folder="FINAL_RESULTS_time_step"


# numbertest=3
# nametest="unsteady_vortex"
# epsilon_vector=[1.0] #[1.0,0.1,0.01]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 7 #6
# starting_elements=25
# # number_of_refinements = 1
# # starting_elements=1600

# numbertest=7
# nametest="2D_Sod"
# epsilon_vector=[1.0] #[1.0,0.1,0.01]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24] #[20,21,22,23,24]
# order_time_vector=[-22]# [2,3,4,5]
# CFL=0.25
# number_of_refinements = 1
# starting_elements=200


# numbertest=8
# nametest="1D_Sod_X"
# epsilon_vector=[1.0] #[1.0,0.1,0.01]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 1
# starting_elements=200


# numbertest=100
# nametest="gresho_vortex"
# epsilon_vector=[1.0] #[1.0,0.1,0.01,0.001]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 6 #6
# starting_elements=25
# # number_of_refinements = 1
# # starting_elements=800

# numbertest=99
# nametest="gresho_vortex_longer_time"
# epsilon_vector=[1e-5]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 1
# starting_elements=100


# numbertest=101
# nametest="smooth_gresho_vortex"
# epsilon_vector=[1.0,0.1,0.01,0.001]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 5 #6
# starting_elements=25
# # number_of_refinements = 1
# # starting_elements=800


# numbertest=200
# nametest="smooth_vortex_eps"
# epsilon_vector=[1.0,0.1,0.01,0.001]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 4 #5
# starting_elements=100
# # number_of_refinements = 1
# # starting_elements=1600



# numbertest=201
# nametest="smooth_vortex_eps_modified"
# epsilon_vector=[1.0,0.1,0.01,0.001]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 5 #6
# starting_elements=100
# # number_of_refinements = 1
# # starting_elements=1600



# numbertest=301
# nametest="baroclinic_vorticity_generation"
# epsilon_vector=[0.05] #[1.0,0.1,0.01]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.1
# number_of_refinements = 1
# starting_elements=800
# NX=800
# NY=160

# numbertest=302
# nametest="kelvin_helmholtz_instability"
# epsilon_vector=[0.001] #[0.1,0.01,0.001]
# scheme_vector=["AF"] #NB: Only for AF
# post_processing_vector=[5] #[0] 
# setting="IMEX" #"IMEX","EX"
# space_reconstruction_vector=[23] #[20,21,22,23,24]
# order_time_vector=[-22] #[2,3,4,5]
# CFL=0.25
# number_of_refinements = 1
# starting_elements=1000
# NX=1000
# NY=500



#########################################
#########################################
#########################################
#CHECK FROM HERE
#########################################
#########################################
#########################################

#######################
#NOT WORKING
#CFL=0.1 is not enough
#NB: I think it is a problem of reconstruction getting negative and not of CFL, the simulation does not crash at the beginning.
#It does not even work for reconstruction of char with rusanov and standard speed estimate
#######################
# numbertest=350
# nametest="woodward_colella"
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=200
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=1
# riemann_solver=0
# speed_estimate=1
# CFL=0.1
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0
# NX=200
# NY=5


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
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =1
# speed_estimate=0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


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
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=1
# riemann_solver        =1
# speed_estimate        =0
# CFL=0.05
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


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
# numbertest=412
# nametest="RP2_relaxed"
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =1
# speed_estimate        =0
# CFL=0.4
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


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
# numbertest=403
# nametest="RP3"
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =0
# speed_estimate        =0
# CFL=0.45
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0

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
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =1
# speed_estimate        =0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


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
# numbertest=405
# nametest="RP5"
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=1
# riemann_solver        =0
# speed_estimate        =0
# CFL=0.1
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0

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
# numbertest=406
# nametest="RP6"
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =1
# speed_estimate        =0
# CFL=0.75
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0

#####################
#CFL 0.95
#100 elements
#####################
# numbertest=500
# nametest="lax"
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=1
# riemann_solver        =1
# speed_estimate=0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


# numbertest=501
# nametest="lax_smaller_lambda"
# orders=[3,5,7,9,11,13]
# meshes = {3:1,5:1,7:1,9:1,11:1,13:1}
# starting_elements=100
# starting_elements_dict = {key: starting_elements for key in orders}  # use elements of 'orders' as keys
# reconstructed_variable=0
# riemann_solver        =0
# speed_estimate=0
# CFL=0.95
# RelaxedCFL=CFL
# NRelaxedTimeSteps=0


####################################################
#NAME OF TIME SCHEME
if setting == "EX":
    time_scheme={1:1,2:12,3:13,4:14,5:15}
elif setting == "IMEX":
    time_scheme={1:-1,2:-2,3:-3,4:-4,5:-5,-12:-12,-22:-22,-32:-32}
else:
    print("No acceptable setting", setting)
    quit()
####################################################


#Clean, Compile
instruction="make clean; make"
os.system(instruction)

#Create simulations folders and copy main there
for epsilon in epsilon_vector:
    for scheme in scheme_vector:
        for post_processing in post_processing_vector:
            for space_reconstruction in space_reconstruction_vector:
                for order_time in order_time_vector:
                    foldName=folder+"_"+variables+"/"+nametest+"/eps"+str(epsilon)+"/"+scheme+"/PP"+str(post_processing)+"/"+setting+"/space_reconstruction_"+str(space_reconstruction)+"/DeC"+str(order_time)+"/CFL"+str(CFL)
                    instruction =""
                    instruction+="mkdir -p "+foldName+" \n" #creation of the folder
                    os.system(instruction)
                    instruction =""
                    instruction+="cp bin/main "+foldName+"/main \n" #copying DATA/don1d
                    os.system(instruction)
                    print("Created:", foldName)


#Run simulations
for epsilon in epsilon_vector:
    for scheme in scheme_vector:
        for post_processing in post_processing_vector:
            for space_reconstruction in space_reconstruction_vector:
                for order_time in order_time_vector:
                    foldName=folder+"_"+variables+"/"+nametest+"/eps"+str(epsilon)+"/"+scheme+"/PP"+str(post_processing)+"/"+setting+"/space_reconstruction_"+str(space_reconstruction)+"/DeC"+str(order_time)+"/CFL"+str(CFL)
                    # print("Number of meshes",number_of_refinements)
                    elements=np.zeros(number_of_refinements)
                    elements[0]=starting_elements
                    for ind in range(1,number_of_refinements):
                        elements[ind]=elements[ind-1]*2
                    instruction =""
                    instruction+="cd "+foldName+" \n" #move in the folder
                    for element in elements:      
                        if nametest in ["baroclinic_vorticity_generation","kelvin_helmholtz_instability"]:
                            instruction+="./main "+str(numbertest)+" "+str(NX)+" "+str(NY)+" "+str(int(space_reconstruction))+" "+str(int(time_scheme[order_time]))+" "+str(CFL)+" "+str(post_processing)+" "+str(epsilon)+" \n"
                        else:
                            instruction+="./main "+str(numbertest)+" "+str(int(element))+" "+str(int(element))+" "+str(int(space_reconstruction))+" "+str(int(time_scheme[order_time]))+" "+str(CFL)+" "+str(post_processing)+" "+str(epsilon)+" \n"
                        
                    print(instruction)
                    os.system(instruction)
