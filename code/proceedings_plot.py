coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
lw=2. #2. for shock interaction
ms=3.0 #1. 
fs=35
#overwritten later if needed


# nametest="composite_wave_200"
# file="SOLUTION_OCT_0000008.0000000.dat"
# CFL_test=0.95 #OK
# rec_var_test       =1
# riemann_solver_test=5
# speed_estimate_test=0
# solution_type="exact solution"


# nametest="composite_wave_longer_time_200"
# file="SOLUTION_OCT_0001000.0000000.dat"
# CFL_test=0.95 #OK
# rec_var_test       =1
# riemann_solver_test=1
# speed_estimate_test=0
# solution_type="exact solution"

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
# nametest="sod"
# file="SOLUTION_OCT_0000000.2000000.dat"
# CFL_test=0.95 #0.8 #original #0.95 #numerical fluxes
# solution_type="exact solution"



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
# nametest="shock_turbulence_interaction_shu_osher"
# file="SOLUTION_OCT_0000001.8000000.dat"
# CFL_test=0.95 #OK
# lw=2. #for shock interaction
# ms=2. 
# solution_type="reference solution"


#######################
#######################
#######################
# nametest="shock_turbulence_modified_example_1_6400"
# file="SOLUTION_OCT_0000005.0000000.dat"
# lw=2. #for shock interaction
# ms=2. 
# CFL_test=0.4
# WhichVariableForLimiting_PCSD =2
# solution_type="reference solution"
# teststocompare=[] 
# teststocompare.append([29,29,3, 0.0, CFL_test]) 
# teststocompare.append([28,29,3, 1.0, CFL_test]) 
# teststocompare.append([28,28,3, 0.0, CFL_test]) 


#######################
#######################
#######################
# nametest="shock_turbulence_modified_example_1_free_BC_3200"
# file="SOLUTION_OCT_0000005.0000000.dat"
# lw=2. #for shock interaction
# ms=2. 
# CFL_test=0.4
# WhichVariableForLimiting_PCSD =2
# solution_type="reference solution"
# teststocompare=[] 
# teststocompare.append([29,29,3, 0.0, CFL_test]) 
# teststocompare.append([28,29,3, 1.0, CFL_test]) 
# teststocompare.append([28,29,3, 2.0, CFL_test]) 
# teststocompare.append([28,29,3, 3.0, CFL_test]) 
# teststocompare.append([28,29,3, 4.0, CFL_test]) 
# teststocompare.append([28,28,3, 0.0, CFL_test]) 


#######################
#######################
#######################
# nametest="shock_turbulence_modified_example_2_3200"
# file="SOLUTION_OCT_0000005.0000000.dat"
# lw=2. #for shock interaction
# ms=2. 
# CFL_test=0.4
# WhichVariableForLimiting_PCSD=2
# solution_type="reference solution"
# teststocompare=[] 
# teststocompare.append([29,29,3, 0.0, CFL_test]) 
# teststocompare.append([28,29,3, 0.015, CFL_test]) 
# teststocompare.append([28,28,3, 0., CFL_test]) 


#######################
#######################
#######################
nametest="shock_turbulence_modified_example_2_free_BC_1200" #NOT 1600
file="SOLUTION_OCT_0000005.0000000.dat"
lw=2. #for shock interaction
ms=2. 
CFL_test=0.4
WhichVariableForLimiting_PCSD=2
solution_type="reference solution"
teststocompare=[] 
teststocompare.append([28,28,3, 0., CFL_test]) 
teststocompare.append([28,29,3, 0.015, CFL_test]) 
# teststocompare.append([29,29,3, 0.0, CFL_test]) 




#######################
#######################
#######################
# nametest="woodward_colella_400"
# file="SOLUTION_OCT_0000000.0380000.dat"
# CFL_test=0.4
# WhichVariableForLimiting_PCSD=2
# solution_type="reference solution"
# teststocompare=[] 
# teststocompare.append([28,28,3, 0., CFL_test]) 
# teststocompare.append([28,29,3, 0.2, CFL_test]) 
# # teststocompare.append([29,29,3, 0.0, CFL_test]) 




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
# nametest="RP1"
# file="SOLUTION_OCT_0000000.2000000.dat"
# CFL_test=0.95 #OK
# solution_type="exact solution"


#######################
#NOT WORKING
# NOT WORKING, actually it works for order 3 only
#CFL=0.1 is not enough
#NB: I think it is a problem of reconstruction getting negative and not of CFL
#######################
#From CFL smaller or equal than 0.85, the test works for order 3 if characteristic variables are reconstructed and Rusanov numerical flux is adopted.
#For all other settings and orders, even for order 3, even for speed estimate with exact Riemann solver, even for smaller CFL up to 0.05, the simulations crash.
#Again, we remark that, even though this is out from the goal of this investigation, the approach may benefit from extra asaptive strategies such sa a posteriori limiting.
#This is why we focused on a relaxed version of the problem.
#######################
# nametest="RP2"
# file="SOLUTION_OCT_0000000.1500000.dat"
# CFL_test=0.25
# solution_type="exact solution"


#####################
#CFL 0.4
#100 elements
#####################
#With reconstruction of characteristic variables, there are no simulation crashes for both numerical fluxes for CFL smaller or equal to 0.7
#With reconstruction of conserved variables variables, the CFL must be decreased.
#In such a case:
#with exact Riemann solver, from CFL=0.7 and CFL=0.6 only order 3 works (even with Riemann speed estimate), from CFL=0.5 to CFL=0.05 only orders 3 and 5 work, 
#with Rusanov, from CFL=0.7 orders 3, 5 and 7 work, for CFL=0.6 orders 3,5,7 and 9 work, for CFL=0.5 all orders work beside order 13, for CFL=0.4 all orders work
#We comment the results for CFL=0.4.
#The reconstruction of characteristic variables outperforms the one one on conserved variables.
#This can be inferred from the plotted profiles for Rusanov numerical flux, while, for exact Riemann solver this is also clear from the fact that the simulation crases if conserved variables are reconstructed from order 7 on.
#Again, the diffusive character of Rusanov is able to clip some spurious over and undershoots.
#Overall, the best results are obtained for Rusanov with reconstruction of characteristic variables and Rusanov.
#Also in this case, the higher is the order, the closer are the results to the reference exact solution.
#####################
# nametest="RP2_relaxed"
# file="SOLUTION_OCT_0000000.1500000.dat"
# CFL_test=0.8  #0.4 #original #0.8 #numerical fluxes
# rec_var_test       =0
# riemann_solver_test=0
# speed_estimate_test=0
# solution_type="exact solution"


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
# nametest="RP3"
# file="SOLUTION_OCT_0000000.0120000.dat"
# CFL_test=0.95 #0.45 #original #0.7 #numerical fluxes
# solution_type="exact solution"


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
# nametest="RP4"
# file="SOLUTION_OCT_0000000.0350000.dat"
# CFL_test=0.95 #OK
# solution_type="exact solution"


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
# nametest="RP5"
# file="SOLUTION_OCT_0000000.0120000.dat"
# CFL_test=0.1 #OK
# solution_type="exact solution"


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
# nametest="RP6"
# file="SOLUTION_OCT_0000000.8000000.dat"
# CFL_test=0.85 #0.75 #original #0.85 #numerical fluxes
# solution_type="exact solution"

#######################
#CFL=0.95
#######################
# nametest="stationary_contact"
# file="SOLUTION_OCT_0000002.0000000.dat"
# # file="SOLUTION_OCT_0000005.0000000.dat"
# # file="SOLUTION_OCT_0000050.0000000.dat"
# # file="SOLUTION_OCT_0000500.0000000.dat"
# # file="SOLUTION_OCT_0005000.0000000.dat"
# CFL_test=0.5 #OK
# solution_type="exact solution"

#######################
#CFL=0.95
#######################
# nametest="moving_contact"
# file="SOLUTION_OCT_0000002.0000000.dat"
# CFL_test=0.5 #OK
# solution_type="exact solution"


#####################
#CFL 0.95
#100 elements
#####################
#All settings work with no problem for CFL=0.95 for standard speed stimate.
#Results are much better with reconstruction of characteristic variables.
#Exact Riemann solver gives better results than Rusanov.
#Overall, the best results are the ones obtained with reconstruction of characteristic variables and exact Riemann solver.
#####################
# nametest="lax"
# file="SOLUTION_OCT_0000001.3000000.dat"
# CFL_test=0.95
# solution_type="exact solution"


# nametest="lax_smaller_lambda"
# file="SOLUTION_OCT_0000001.3000000.dat"
# CFL_test=0.95
# solution_type="exact solution"

if not(nametest.startswith("shock_turbulence_interaction")):
    lwex=lw+1.5
else:
    lwex=lw


gmm=1.4



#Reconstructed variable, riemann_solver, speed_estimate, order, CFL, RelaxedCFL, NRelaxedTimeSteps



#Reconstructed variable, riemann_solver, speed_estimate, order, CFL, RelaxedCFL, NRelaxedTimeSteps
# teststocompare.append([28,29,3, 0.0, CFL_test]) 
# teststocompare.append([28,29,3, 0.025, CFL_test]) 
# teststocompare.append([28,29,3, 0.05, CFL_test]) 
# teststocompare.append([28,29,3, 0.075, CFL_test]) 
# teststocompare.append([28,29,3, 0.1, CFL_test]) 
# teststocompare.append([28,29,3, 0.5, CFL_test]) 
# teststocompare.append([28,29,3, 1.0, CFL_test]) 

# teststocompare.append([28,28,3, 0.1, CFL_test]) 
# teststocompare.append([28,29,3, 0.1, CFL_test]) 
# teststocompare.append([29,29,3, 0.1, CFL_test]) 
# teststocompare.append([28,28,3, 0.25, CFL_test]) 
# teststocompare.append([28,29,3, 0.25, CFL_test]) 
# teststocompare.append([29,29,3, 0.25, CFL_test]) 
# teststocompare.append([28,28,3, 0.5, CFL_test]) 
# teststocompare.append([28,29,3, 0.5, CFL_test]) 
# teststocompare.append([29,29,3, 0.5, CFL_test]) 
# teststocompare.append([28,28,3, 0.75, CFL_test]) 
# teststocompare.append([28,29,3, 0.75, CFL_test]) 
# teststocompare.append([29,29,3, 0.75, CFL_test]) 

# teststocompare.append([28,28,3, 0.5, CFL_test]) 

# teststocompare.append([29,29,3, 0.0, CFL_test]) 

# teststocompare.append([28,29,3, 0.015, CFL_test]) 
# teststocompare.append([28,29,3, 0.02 , CFL_test]) 
# teststocompare.append([28,29,3, 0.025, CFL_test]) 
# teststocompare.append([28,29,3, 1.0  , CFL_test]) 

# teststocompare.append([28,28,3, 0., CFL_test]) 



import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl


linestyles_space_reconstruction = {1: "-",
                                   10:"-",
                                   2: "-",
                                   20:"-",
                                   21:"-",
                                   22:"-",
                                   23:"-",
                                   24:"-",
                                   25:"-",
                                   26:"-",
                                   27:"-",
                                   28:"-",
                                   29:"-",
                                  -1: "-",
                                  -10:"-",
                                  -2: "-",
                                  -20:"-",
                                  -21:"-",
                                  -22:"-",
                                  -23:"-",
                                  -24:"-",
                                  -25:"-"}

markers_space_reconstruction   = { 1: "x",
                                   10:"x",
                                   2: "x",
                                   20:"x",
                                   21:"x",
                                   22:"x",
                                   23:"x",
                                   24:"x",
                                   25:"x",
                                   26:"x",
                                   27:"x",
                                   28:"x",
                                   29:"x",
                                  -1: "o",
                                  -10:"o",
                                  -2: "o",
                                  -20:"o",
                                  -21:"o",
                                  -22:"o",
                                  -23:"o",
                                  -24:"o",
                                  -25:"o",
                                  -26:"o",
                                  -27:"o",
                                  -28:"o",
                                  -29:"o"}


#Updated col
colors = {
     1:"#17becf",  #PWC
     2:"#d62728",  #MINMOD
     20:"#d62728", #MINMOD
     21:"#1f77b4", #
     22:"#ff7f0e", #
     23:"#2ca02c", #
     24:"#9467bd",  #
     25:"#98df8a",   # New
     26:"#c5b0d5",   # New
     27:"#c49c94",   # New
     28: "red",   # New
     29: "blue",   # New


    -1:"#17becf",  #PWC
    -2:"#d62728",  #MINMOD
    -20:"#d62728", #MINMOD
    -21:"#1f77b4", #
    -22:"#ff7f0e", #
    -23:"#2ca02c", #
    -24:"#9467bd",  #



     93:"#e377c2", #

     6:"#8c564b",  #
    -2:"#bcbd22",  #

    -6:"#ff9896",  #
     7:"#7f7f7f"   #
}

# colors = {
#     3:  "#1f77b4",   # 
#     5:  "#ff7f0e",   # 
#     7:  "#2ca02c",   # 
#     9:  "#d62728",   # 
#     11: "#9467bd",   # 
#     13: "#8c564b",   # 
#     15: "#17becf",   # 
#     17: "#bcbd22",   # 
#     19: "#e377c2",   # 
#     21: "#7f7f7f",   # 
#     23: "#ff9896",   # 
#     25: "#98df8a",   # New
#     27: "#c5b0d5",   # New
#     29: "#c49c94",   # New
#     31: "#f7b6d2",   # New
#     33: "#aec7e8",   # New
#     35: "#ffbb78",   # New
#     37: "#9edae5",   # New
#     39: "#dbdb8d",   # New
#     41: "#c7c7c7",   # New
#     43: "#bc80bd",   # New
#     45: "#6b6ecf",   # New
#     47: "#b5cf6b",   # New
#     49: "#e7969c",   # New
#     51: "#d6616b",   # New
#     53: "#843c39",   # New
#     55: "#8c6d31",   # New
#     57: "#d9d9d9",   # New
#     59: "#7b4173",   # New
#     61: "#b15928",   # New
#     63: "#6baed6"    # New
# }

####################################################
name_reconstructed_variable = {0:"cons",    1:"char"}
name_riemann_solver         = {-2:"central", -1:"LF", 0:"rusanov", 1:"exact", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:"FORCEalpha"}
name_speed_estimate         = {0:"standard",1:"riemann"}
####################################################

####################################################
label_reconstructed_variable = {0:"cons.",    1:"char."}
label_riemann_solver         = {-2:"Central", -1:"LxF", 0:"Rus", 1:"Ex.RS", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:r"FORCE-$\alpha$"}
label_time_scheme            =   { 1: "Euler"          ,
                                   2: "SSPRK2"         ,
                                   3: "SSPRK3"         ,
                                   4: "SSPRK4"         ,
                                   5: "RK65"           ,
                                   12:"SSPRK2"         , #NB:DeC2 is SSPRK2
                                   13:"DeC3"           ,
                                   14:"DeC4"           ,
                                   15:"DeC5" }
####################################################

nVar=3


params = {'mathtext.default': 'regular' } #parameters plot




fig = pl.figure(1, figsize=(10,10))
pl.title( 'density',fontsize=fs)
pl.xlabel('x',fontsize=fs)

firsttime=True
index_plot=0




##############################
# READ AND PLOT EXACT SOLUTION
##############################
# FILLING x_ref AND ref_sol
##############################
if (nametest.startswith("composite_wave")):
    solfoldName="reference_solutions/composite_wave"
elif (nametest.startswith("woodward_colella")):
    solfoldName="reference_solutions/woodward_colella_200000"
elif (nametest.startswith("shock_turbulence_modified_new")):
    solfoldName="reference_solutions/shock_turbulence_modified"
elif (nametest.startswith("shock_turbulence_modified_example_1_free_BC")):
    solfoldName="reference_solutions/shock_turbulence_interaction_example_1_free_BCs_200000"
elif (nametest.startswith("shock_turbulence_modified_example_2_free_BC")):
    solfoldName="reference_solutions/shock_turbulence_interaction_example_2_free_BCs_200000"
elif (nametest.startswith("shock_turbulence_modified_example_1")):
    solfoldName="reference_solutions/shock_turbulence_interaction_example_1_inflow_outflow_200000"
elif (nametest.startswith("shock_turbulence_modified_example_2")):
    solfoldName="reference_solutions/shock_turbulence_interaction_example_2_inflow_outflow_200000"
else:
    solfoldName="reference_solutions/"+nametest
    
if (1==1): #IN CASE YOU DO NOT WANT THE REFERENCE
    if os.path.isdir(solfoldName):  #CONDITION: Is it a folder? If yes go on
        namefile=solfoldName+'/'+"reference.out"
        if os.path.isfile(namefile): #If the file exists, read it
            print("You are in the folder "+solfoldName+" dealing with "+namefile)
            lines=[]
            with open(namefile, 'r') as a:
                for idx, line in enumerate(a.readlines()):
                    lines.append(line)

            x_ref=np.array([])
            ro_ref=np.array([])
            q_ref=np.array([])
            E_ref=np.array([])
            v_ref=np.array([])
            p_ref=np.array([])

            #Loop over the lines that we have read
            for idx, line in enumerate(lines): #rmk: 0-based numeration
                data=line.split()
                x_ref=np.append(x_ref,float(data[0]))
                ro_ref=np.append(ro_ref,float(data[1]))
                v_ref=np.append(v_ref,float(data[2]))
                p_ref=np.append(p_ref,float(data[3]))
                # E_ref=np.append(E_ref,float(data[4]))

            q_ref=ro_ref*v_ref
            E_ref=p_ref/(gmm-1)+0.5*ro_ref*v_ref**2

            pl.plot(x_ref,ro_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            index_plot=index_plot+1


if WhichVariableForLimiting_PCSD==1:
    name_base_folder="latest_adaptivity_based_on_rho"
elif WhichVariableForLimiting_PCSD==2:
    name_base_folder="latest_adaptivity_based_on_rhou"
elif WhichVariableForLimiting_PCSD==3:
    name_base_folder="latest_adaptivity_based_on_energy"
else:
    print("Wrong choice for WhichVariableForLimiting_PCSD")
    quit()

for indt, test in enumerate(teststocompare): #Loop on the schemes
    space_reconstruction    =test[0]
    space_reconstruction_fix=test[1]
    time_scheme             =test[2]
    K_coefficient         =test[3]
    CFL                     =test[4]

    if K_coefficient==0:
        foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/space_reconstruction_fix"+str(space_reconstruction_fix)+"/time_scheme_"+str(time_scheme)+"/K0.0"+"/CFL"+str(CFL)
    else:
        foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/space_reconstruction_fix"+str(space_reconstruction_fix)+"/time_scheme_"+str(time_scheme)+"/K"+str(K_coefficient)+"/CFL"+str(CFL)


    print()
    print("trying to enter",foldName)

    if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
        print("entered",foldName)
        filename = "MESH_Coordinate"+coordinate+".dat"
        if os.path.isfile(foldName+"/"+ filename):
            delimiter_in = ' '
            headerlines_in = 1
            mydata_mesh = np.loadtxt(foldName+"/"+ filename, delimiter=delimiter_in, skiprows=headerlines_in)
            z = mydata_mesh
            z_ini = min(z)
            z_end = max(z)
            Nz    = len(z)

            # Importing the SOLUTION DATA
            if os.path.isfile(foldName+"/"+ file):

                mydata_solution = np.loadtxt(foldName+"/"+ file, skiprows=headerlines_in)

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

    
                # pl.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=str(space_reconstruction)+" "+str(space_reconstruction_fix)+" "+label_time_scheme[time_scheme],color=colors[space_reconstruction_fix],alpha=0.7,zorder=3)
                if space_reconstruction==28 and space_reconstruction_fix==28:
                    pl.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label="LDCU",alpha=0.7,zorder=3)
                elif space_reconstruction==29 and space_reconstruction_fix==29:
                    pl.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label="o-LDCU",alpha=0.7,zorder=3)
                elif space_reconstruction==28 and space_reconstruction_fix==29:
                    pl.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label="a-LDCU K="+str(K_coefficient),alpha=0.7,zorder=3)
                else:
                    print("NO,I DO NO WANT TO PLOT THIS SETTING", space_reconstruction, space_reconstruction_fix)
                    quit()

                if (nametest.startswith("woodward_colella")):
                    pl.xlim([0.55,0.85])
                elif (nametest.startswith("shock_turbulence_modified_example_1")):
                    pl.xlim([-2,0.0])
                elif (nametest.startswith("shock_turbulence_modified_example_2")):
                    pl.xlim([9,11])
                else:
                    pass

                index_plot=index_plot+1

            




pl.grid()





# Create the legend with reordered handles and labels
pl.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)

pl.savefig(nametest+"_compare_plot_only_density_CFL"+str(CFL)+".pdf", format="pdf", bbox_inches="tight")
# pl.savefig(nametest+"_CFL"+str(CFL)+"_compare_char_plot.pdf", format="pdf", bbox_inches="tight")
# pl.show()

