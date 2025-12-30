coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
lw=2. #2. for shock interaction
ms=3.0 #1. 
fs=10
#overwritten later if needed

name_time_scheme={-52:"IMEX_DeC2_Prim_Crank_Nicolson_Cons",-82:"IMEX_DeC2_Prim_Explicit_Cons"}

eps=1.0
nametest="sod"+str(eps)
file    ="WC_SOLUTION_OCT_0000000.2000000.dat"
ref_file="SOLUTION_OCT_0000000.2000000.dat"
CFL_test=0.475 #0.2
solution_type="reference solution"
xspan=[0,1]
teststocompare=[] 
#epsilon, space_reconstruction, time_scheme, K_coefficient, CFL
teststocompare.append([eps, 27, -82, 0.0, CFL_test]) 

# eps=1.0
# nametest="sod"+str(eps)+"_with_timestep_relaxation"
# file    ="WC_SOLUTION_OCT_0000000.2000000.dat"
# ref_file="SOLUTION_OCT_0000000.2000000.dat"
# CFL_test=0.475 #0.2
# solution_type="reference solution"
# xspan=[0,1]
# teststocompare=[] 
# #epsilon, space_reconstruction, time_scheme, K_coefficient, CFL
# teststocompare.append([eps, 27, -82, 0.0, CFL_test]) 


# eps=0.6
# nametest="sod"+str(eps)
# file    ="WC_SOLUTION_OCT_0000000.0500000.dat"
# ref_file="SOLUTION_OCT_0000000.0500000.dat"
# CFL_test=0.475 #0.2
# solution_type="reference solution"
# xspan=[0,1]
# teststocompare=[] 
# #epsilon, space_reconstruction, time_scheme, K_coefficient, CFL
# teststocompare.append([eps, 27, -82, 0.0, CFL_test]) 


# eps=0.3
# nametest="sod"+str(eps)+"_with_timestep_relaxation"
# file    ="WC_SOLUTION_OCT_0000000.0500000.dat"
# ref_file="SOLUTION_OCT_0000000.0500000.dat"
# CFL_test=0.475 #0.2
# solution_type="reference solution"
# xspan=[0,1]
# teststocompare=[] 
# #epsilon, space_reconstruction, time_scheme, K_coefficient, CFL
# teststocompare.append([eps, 27, -82, 0.0, CFL_test]) 








if not(nametest.startswith("shock_turbulence_interaction")):
    lwex=lw+1.5
else:
    lwex=lw


gmm=1.4




import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl

pl.rcParams['xtick.labelsize'] = 4
pl.rcParams['ytick.labelsize'] = 4


linestyles_space_reconstruction = {1: "--",
                                   10:"--",
                                   2: "--",
                                   20:"--",
                                   21:"--",
                                   22:"--",
                                   23:"--",
                                   24:"--",
                                   25:"--",
                                  -1: "-",
                                  -10:"-",
                                  -2: "-",
                                  -20:"-",
                                  -21:"-",
                                  -22:"-",
                                  -23:"-",
                                  -24:"-",
                                  -25:"-",
                                  -26:"-",
                                  -27:"-"}

markers_space_reconstruction   = { 1: "x",
                                   10:"x",
                                   2: "x",
                                   20:"x",
                                   21:"x",
                                   22:"x",
                                   23:"x",
                                   24:"x",
                                   25:"x",
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
                                  -27:"o"}


#Updated col
colors = {
     1:"#17becf",  #PWC
     2:"#d62728",  #MINMOD
     20:"#d62728", #MINMOD
     21:"#1f77b4", #
     22:"#ff7f0e", #
     23:"#2ca02c", #
     24:"#9467bd",  #

    -1:"#17becf",  #PWC
    -2:"#d62728",  #MINMOD
    -20:"#d62728", #MINMOD
    -21:"#1f77b4", #
    -22:"#ff7f0e", #
    -23:"#2ca02c", #
    -24:"#9467bd",  #
    -26:"red", #AAAAAAAAAAAAAAAAAAAAAA
    -27:"red", #AAAAAAAAAAAAAAAAAAAAAA

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
order_of_plotting           = {1:9,2:8,3:7,5:6,7:5,9:4,11:3,13:2}
####################################################


nVar=3


params = {'mathtext.default': 'regular' } #parameters plot




fig = pl.figure(1, figsize=(20,5))
ax_ro = fig.add_subplot(131)
ax_ro.set_title( r'$\rho$',fontsize=fs+5)
# ax_ro.set_xlabel('x',fontsize=fs)
# ax_ro.set_ylabel('density',fontsize=fs)
ax_v  = fig.add_subplot(132)
ax_v.set_title( 'u',fontsize=fs+5)
# ax_v.set_xlabel('x',fontsize=fs)
# ax_v.set_ylabel('velocity',fontsize=fs)
ax_p  = fig.add_subplot(133)
ax_p.set_title( 'p',fontsize=fs+5)
# ax_p.set_xlabel('x',fontsize=fs)
# ax_p.set_ylabel('pressure',fontsize=fs)
# ax_q  = fig.add_subplot(234)
# ax_q.set_title( 'momentum')
# ax_q.set_xlabel('x')
# ax_q.set_ylabel('momentum')
# ax_E  = fig.add_subplot(235)
# ax_E.set_title( 'total energy')
# ax_E.set_xlabel('x')
# ax_E.set_ylabel('total energy')

firsttime=True
index_plot=0

ax_ro.tick_params(axis='both', which='major', labelsize=15)
ax_v.tick_params(axis='both', which='major',  labelsize=15)
ax_p.tick_params(axis='both', which='major',  labelsize=15)



##############################
# READ AND PLOT EXACT SOLUTION
##############################
# FILLING x_ref AND ref_sol
##############################
if (nametest.startswith("composite_wave")):
    solfoldName="reference_solutions/composite_wave"
elif (nametest.startswith("woodward_colella")):
    solfoldName="reference_solutions/woodward_colella_200000_theta1.3_char_LDCU_CFL0.25"
elif (nametest.startswith("RP2")):
    solfoldName="reference_solutions/RP2"
elif (nametest.startswith("shock_turbulence_interaction_shu_osher")):
    solfoldName="reference_solutions/shock_turbulence_interaction_shu_osher_200000_theta1.3_char_LDCU_CFL0.25"
elif (nametest.startswith("sod0.3")):
    solfoldName="reference_solutions/sod0.3_2000"
elif (nametest.startswith("sod0.6")):
    solfoldName="reference_solutions/sod0.6_2000"
elif (nametest.startswith("sod0.9")):
    solfoldName="reference_solutions/sod0.9_2000"
elif (nametest.startswith("sod1.0")):
    solfoldName="reference_solutions/sod1.0_2000"
elif (nametest.startswith("sod2d0.3")):
    solfoldName="reference_solutions/sod2d0.3_2000"
elif (nametest.startswith("sod2d0.6")):
    solfoldName="reference_solutions/sod2d0.6_2000"
elif (nametest.startswith("sod2d0.9")):
    solfoldName="reference_solutions/sod2d0.9_2000"
elif (nametest.startswith("sod2d1.0")):
    solfoldName="reference_solutions/sod2d1.0_2000"
else:
    solfoldName="reference_solutions/"+nametest

    
if (1==1): #IN CASE YOU DO NOT WANT THE REFERENCE
    if os.path.isdir(solfoldName):  #CONDITION: Is it a folder? If yes go on
        namefile=solfoldName+'/'+ref_file
        filename = solfoldName+'/'+"MESH_Coordinate"+coordinate+".dat"
        if os.path.isfile(filename):

            delimiter_in = ' '
            headerlines_in = 1
            mydata_mesh = np.loadtxt(filename, delimiter=delimiter_in, skiprows=headerlines_in)
            z = mydata_mesh
            z_ini = min(z)
            z_end = max(z)
            Nz    = len(z)

            # Importing the SOLUTION DATA
            if os.path.isfile(namefile):

                mydata_solution = np.loadtxt(namefile, skiprows=headerlines_in)

                ro   = mydata_solution[:, 0]
                rou    = mydata_solution[:, 1]
                rov    = mydata_solution[:, 2]
                Energy    = mydata_solution[:, 3]
                phi  = mydata_solution[:, 4]

                RO     = ro[:Nz]
                ROU      = rou[:Nz]
                ROV      = rov[:Nz]
                ENERGY      = Energy[:Nz]
                PHI    = phi[:Nz]


            x_ref=z
            ro_ref=RO
            q_ref=ROU
            E_ref=ENERGY
            v_ref=ROU/RO
            p_ref=(gmm-1)*(ENERGY-0.5*ROU**2/RO*eps**2)

            ax_ro.plot(x_ref,ro_ref,linestyle="-",linewidth=lw,label=solution_type,color="k")
            ax_v.plot( x_ref, v_ref,linestyle="-",linewidth=lw,label=solution_type,color="k")
            ax_p.plot( x_ref, p_ref,linestyle="-",linewidth=lw,label=solution_type,color="k")
            # ax_q.plot( x_ref, q_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            # ax_E.plot( x_ref, E_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            index_plot=index_plot+1


name_base_folder="IMEX_HYPERBOLIC_TRICK_RHO_P_STAR_STAGE_DEPENDENT"


for indt, test in enumerate(teststocompare): #Loop on the schemes
    epsilon=test[0]
    space_reconstruction = test[1]
    time_scheme          = test[2]
    K_coefficient        = test[3]
    CFL                  = test[4]

    foldName=name_base_folder+"/"+nametest+"/eps"+str(epsilon)+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/K"+str(K_coefficient)+"/CFL"+str(CFL)
    print(foldName)
    if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
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
                u    = mydata_solution[:, 1]
                v    = mydata_solution[:, 2]
                p    = mydata_solution[:, 3]
                phi  = mydata_solution[:, 4]

                RO     = ro[:Nz]
                U      = u[:Nz]
                V      = v[:Nz]
                P      = p[:Nz]
                PHI    = phi[:Nz]


                ax_ro.plot(z,RO,marker=".",linestyle="-",linewidth=lw, markersize=ms,label="AP sheme",color="red",alpha=0.7,zorder=3)
                ax_v.plot( z, U,marker=".",linestyle="-",linewidth=lw, markersize=ms,label="AP sheme",color="red",alpha=0.7,zorder=3)
                ax_p.plot( z, P,marker=".",linestyle="-",linewidth=lw, markersize=ms,label="AP sheme",color="red",alpha=0.7,zorder=3)
                ax_ro.set_xlim(xspan)
                ax_v.set_xlim(xspan)
                ax_p.set_xlim(xspan)
                index_plot=index_plot+1
            




# ax_ro.grid()
# ax_v.grid()
# ax_p.grid()
# ax_q.grid()
# ax_E.grid()
# pl.grid()

# Capturing handles and labels
handles, labels = pl.gca().get_legend_handles_labels()
# Define the desired order
desired_order=np.arange(index_plot)
# print(np.flip(desired_order))
desired_order=np.flip(desired_order)

# Reorder handles and labels
handles[:] = [handles[i] for i in desired_order[:]]
labels[:] = [labels[i] for i in desired_order[:]]



# Create the legend with reordered handles and labels
pl.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)

pl.savefig(nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_V_CFL"+str(CFL)+".pdf", format="pdf", bbox_inches="tight")
# pl.savefig(nametest+"_CFL"+str(CFL)+"_compare_char_plot.pdf", format="pdf", bbox_inches="tight")
# pl.show()

fig = pl.gcf()
fig.canvas.draw()  # Necessary for accurate bbox

def save_axis(ax, filename):
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, bbox_inches=extent)


ax_ro.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
ax_v.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
ax_p.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)

quit()

# save_axis(ax_ro, nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_V_CFL"+str(CFL)+"_rho.pdf")
# save_axis(ax_v, nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_V_CFL"+str(CFL)+"_u.pdf")
# save_axis(ax_p, nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_V_CFL"+str(CFL)+"_p.pdf")


save_axis(ax_ro, nametest+"_V_CFL"+str(CFL)+"_rho.pdf")
save_axis(ax_v, nametest+"_V_CFL"+str(CFL)+"_u.pdf")
save_axis(ax_p, nametest+"_V_CFL"+str(CFL)+"_p.pdf")