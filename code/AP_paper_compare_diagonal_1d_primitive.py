teststocompare=[] 

# eps=1.0
# nametest="sod2d1.0"
# solution_file="WC_SOLUTION_OCT_0000000.2500000.dat"
# ref_file="SOLUTION_OCT_0000000.2500000.dat"
# coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
# time_scheme=-82
# CFL_test=0.475 #OK
# solution_type="reference solution"

# eps=1.0
# nametest="sod2d1.0_with_timestep_relaxation"
# solution_file="WC_SOLUTION_OCT_0000000.2500000.dat"
# ref_file="SOLUTION_OCT_0000000.2500000.dat"
# coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
# time_scheme=-52
# CFL_test=0.475 #OK
# solution_type="reference solution"


eps=0.3 #0.3, 0.6, 0.9
nametest="sod2d"+str(eps)
solution_file="WC_SOLUTION_OCT_0000000.0500000.dat"
ref_file="SOLUTION_OCT_0000000.0500000.dat"
coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
time_scheme=-82
CFL_test=0.1 #OK
solution_type="reference solution"


# eps=0.3 #0.3, 0.6, 0.9
# nametest="sod2d"+str(eps)+"_with_timestep_relaxation"
# solution_file="WC_SOLUTION_OCT_0000000.0500000.dat"
# ref_file="SOLUTION_OCT_0000000.0500000.dat"
# coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
# time_scheme=-82
# CFL_test=0.475 #OK
# solution_type="reference solution"


gmm=1.4


which_plot="diagonal" #"diagonal", "semi_diagonal"

teststocompare=[] 
#epsilon, space_reconstruction, time_scheme, K_coefficient, CFL
teststocompare.append([eps, 27, time_scheme, 0.0, CFL_test]) 
name_time_scheme={-52:"IMEX_DeC2_Prim_Crank_Nicolson_Cons",-82:"IMEX_DeC2_Prim_Explicit_Cons"}


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

####################################################
label_reconstructed_variable = {0:"cons.",    1:"char."}
label_riemann_solver         = {-2:"Central", -1:"LxF", 0:"Rus", 1:"Ex.RS", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:r"FORCE-$\alpha$"}
label_space_reconstruction   =   { 1: "FV-o1"          ,
                                   10:"FV-o1-PWL in PC",
                                   2: "FV-o2-MUSCL"    ,
                                   20:"FV-o2-MUSCL"    ,
                                   21:"FV-o2-k-MUSCL"  ,
                                   22:"FV-o2-CO-MUSCL" ,
                                   23:"FV-o2-VL-MUSCL" ,
                                   24:"FV-o2-M-MUSCL"  ,
                                   25:"FV-o2-VA-MUSCL" ,
                                  -1: "AF-o1"          ,
                                  -10:"AF-o1-PWL in PC",
                                  -2: "AF-o2-MUSCL"    ,
                                  -20:"AF-o2-MUSCL"    ,
                                  -21:"AF-o2-k-MUSCL"  ,
                                  -22:"AF-o2-CO-MUSCL" ,
                                  -23:"AF-o2-VL-MUSCL" ,
                                  -24:"AF-o2-M-MUSCL"  ,
                                  -25:"AF-o2-VA-MUSCL" ,
                                  -26:"AF-o2-SBM",
                                  -27:"FV AF" }
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

nVar=4


params = {'mathtext.default': 'regular' } #parameters plot


import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl

#test="../AnotherFolder";
test="./"; #Actual folder, not used but in case you want to run it from a specific folder you can use it to implement it

fig, axs = pl.subplots(1,3, figsize=(13, 4)) #Array of subplots
# fig.suptitle("")
firsttime=False

col=0
#Plotting
##############################
# READ AND PLOT EXACT SOLUTION
##############################
# FILLING x_ref AND ref_sol
##############################
if (nametest.startswith("composite_wave")):
    foldName="reference_solutions/composite_wave"
elif (nametest.startswith("woodward_colella")):
    foldName="reference_solutions/woodward_colella_200000_theta1.3_char_LDCU_CFL0.25"
elif (nametest.startswith("RP2")):
    foldName="reference_solutions/RP2"
elif (nametest.startswith("shock_turbulence_interaction_shu_osher")):
    foldName="reference_solutions/shock_turbulence_interaction_shu_osher_200000_theta1.3_char_LDCU_CFL0.25"
elif (nametest.startswith("sod0.3")):
    foldName="reference_solutions/sod0.3_2000"
elif (nametest.startswith("sod0.6")):
    foldName="reference_solutions/sod0.6_2000"
elif (nametest.startswith("sod0.9")):
    foldName="reference_solutions/sod0.9_2000"
elif (nametest.startswith("sod1.0")):
    foldName="reference_solutions/sod1.0_2000"
elif (nametest.startswith("sod2d0.3")):
    foldName="reference_solutions/sod2d0.3_2000"
elif (nametest.startswith("sod2d0.6")):
    foldName="reference_solutions/sod2d0.6_2000"
elif (nametest.startswith("sod2d0.9")):
    foldName="reference_solutions/sod2d0.9_2000"
elif (nametest.startswith("sod2d1.0")):
    foldName="reference_solutions/sod2d1.0_2000"
else:
    foldName="reference_solutions/"+nametest


delimiter_in = ' '
headerlines_in = 1


index_plot=0

fs=10
lwex=2.0

axs[0].tick_params(axis='both', which='major', labelsize=15)
axs[1].tick_params(axis='both', which='major', labelsize=15)
axs[2].tick_params(axis='both', which='major', labelsize=15)

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
    filename=ref_file
    if os.path.isfile(foldName+"/"+ filename):

        mydata_solution = np.loadtxt(foldName+"/"+ filename, skiprows=headerlines_in)

        ro     = mydata_solution[:, 0]
        rou    = mydata_solution[:, 1]
        rov    = mydata_solution[:, 2]
        energy = mydata_solution[:, 3]
        phi    = mydata_solution[:, 4]

        RO     = np.reshape(ro,  (X.shape[0], Y.shape[1]))
        ROU    = np.reshape(rou,   (X.shape[0], Y.shape[1]))
        ROV    = np.reshape(rov,   (X.shape[0], Y.shape[1]))
        ENERGY = np.reshape(energy,   (X.shape[0], Y.shape[1]))
        PHI    = np.reshape(phi, (X.shape[0], Y.shape[1]))

        U=ROU/RO
        V=ROV/RO

        P=(gmm-1.0)*(ENERGY-0.5*RO*(U**2+V**2)*eps**2)

        RO=np.diag(RO)
        U=np.diag(U)
        V=np.diag(V)
        P=np.diag(P)
        modV=np.sqrt(U**2+V**2)
        z=np.diag(X)


        modROV=np.sqrt(ROU**2+ROV**2)

        if which_plot=="semi_diagonal":
            RO   = RO[z > 0]
            U    = U[z > 0]
            V    = V[z > 0]
            P    = P[z > 0]
            modV = modV[z > 0]
            z    = z[z>0]

        #passing to the diagonal
        z    = np.sqrt(2)*z

        lw=0.5
        ms=1

        axs[0].plot(z,RO,label=solution_type,linewidth=lwex,color="k")
        axs[0].set_title(r"$\rho$",fontsize=fs+5)
        # axs[0].set_xlabel("r:=$\sqrt{2}x$=$\sqrt{2}y$")
        # axs[0].grid()
        # axs[0].set_ylabel("y")
        # axs[0].legend("H")
        axs[0].set_xticks([])  # Remove all x-axis tick labels

        #u
        axs[1].plot(z,modV,label=solution_type,linewidth=lwex,color="k")
        axs[1].set_title(r"$\sqrt{u^2+v^2}$",fontsize=fs+5)
        # axs[1].set_xlabel("r:=$\sqrt{2}x$=$\sqrt{2}y$")
        # axs[1].grid()
        axs[1].set_xticks([])  # Remove all x-axis tick labels


        #q
        axs[2].plot(z,P,label=solution_type,linewidth=lwex,color="k")
        axs[2].set_title("p",fontsize=fs+5)
        # axs[2].set_xlabel("r:=$\sqrt{2}x$=$\sqrt{2}y$")
        # axs[2].grid()
        axs[2].set_xticks([])  # Remove all x-axis tick labels

        index_plot=index_plot+1
# pl.show()
# quit()

# name_base_folder="IMEX_DEFINITIVE"
# name_base_folder="IMEX_NO_HYPERBOLIC_TRICK"
name_base_folder="IMEX_HYPERBOLIC_TRICK_RHO_P_STAR_STAGE_DEPENDENT"

for indt, test in enumerate(teststocompare): #Loop on the schemes
    epsilon=test[0]
    space_reconstruction = test[1]
    time_scheme          = test[2]
    K_coefficient        = test[3]
    CFL                  = test[4]

    foldName=name_base_folder+"/"+nametest+"/eps"+str(epsilon)+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/K"+str(K_coefficient)+"/CFL"+str(CFL)
    print(foldName)


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
            u    = mydata_solution[:, 1]
            v    = mydata_solution[:, 2]
            p = mydata_solution[:, 3]
            phi  = mydata_solution[:, 4]

            RO  = np.reshape(ro , (X.shape[0], Y.shape[1]))
            U   = np.reshape(u  , (X.shape[0], Y.shape[1]))
            V   = np.reshape(v  , (X.shape[0], Y.shape[1]))
            P   = np.reshape(p  , (X.shape[0], Y.shape[1]))
            PHI = np.reshape(phi, (X.shape[0], Y.shape[1]))

            ROU=RO*U
            ROV=RO*V

            ENERGY=P/(gmm-1)+0.5*RO*eps**2*(U**2+V**2)

            RO=np.diag(RO)
            U=np.diag(U)
            V=np.diag(V)
            ROU=np.diag(ROU)
            ROV=np.diag(ROV)
            ENERGY=np.diag(ENERGY)
            P=np.diag(P)
            modV=np.sqrt(U**2+V**2)
            z=np.diag(X)
            modROV=np.sqrt(ROU**2+ROV**2)


            if which_plot=="semi_diagonal":
                RO   = RO[z > 0]
                U    = U[z > 0]
                V    = V[z > 0]
                P    = P[z > 0]
                modV = modV[z > 0]
                z    = z[z>0]


            #passing to the diagonal
            z    = np.sqrt(2)*z

            lw=1.
            ms=2

            axs[0].plot(z,RO,marker=".",linestyle="-",linewidth=lw, markersize=ms,label="AP sheme",color="red",alpha=0.7,zorder=3)
            axs[0].set_title(r"$\rho$",fontsize=fs)
            # axs[0].set_xlabel("r:=$\sqrt{2}x$=$\sqrt{2}y$")
            axs[0].set_xlim([np.min(z),np.max(z)])
            # axs[0].grid()
            # axs[0].set_ylabel("y")
            # axs[0].legend("H")

            #u
            axs[1].plot(z,modV,marker=".",linestyle="-",linewidth=lw, markersize=ms,label="AP sheme",color="red",alpha=0.7,zorder=3)
            axs[1].set_title(r"$\sqrt{u^2+v^2}$",fontsize=fs)
            # axs[1].set_xlabel("r:=$\sqrt{2}x$=$\sqrt{2}y$")
            axs[1].set_ylim([-0.05,np.max(modV)*1.3])
            axs[1].set_xlim([np.min(z),np.max(z)])
            # axs[1].grid()


            #q
            axs[2].plot(z,P,marker=".",linestyle="-",linewidth=lw, markersize=ms,label="AP sheme",color="red",alpha=0.7,zorder=3)
            axs[2].set_title("p",fontsize=fs)
            # axs[2].set_xlabel("r:=$\sqrt{2}x$=$\sqrt{2}y$")
            axs[2].set_xlim([np.min(z),np.max(z)])
            axs[2].set_ylim([0.0,np.max(P)*1.2])
            # axs[2].grid()

            index_plot=index_plot+1
 
fig.tight_layout()
# for row in axs: 
#     row.grid(which='major', color='#666666', linestyle='-',alpha=0.2)
#     row.minorticks_on()  
#     row.grid(which='minor', color='#666666', linestyle='-',alpha=0.2) 

# Move the legend outside the plot

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
pl.legend(handles, labels,fontsize=fs-1,loc='upper right')#loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs+2)


params = {'mathtext.default': 'regular' }   
pl.savefig(nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_CFL"+str(CFL_test)+"_V_compare_diag.pdf", format="pdf", bbox_inches="tight")
# # pl.savefig("plotting_compare.png",dpi=600)
# pl.show()

#SAVE FIGURES INDIVIDUALLY
from matplotlib.transforms import Bbox

fig = pl.gcf()
fig.canvas.draw()  # Necessary for accurate bbox

def save_axis(ax, filename):
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, bbox_inches=extent)


axs[0].legend(handles, labels,fontsize=fs-1,loc='upper right')#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
axs[1].legend(handles, labels,fontsize=fs-1,loc='upper right')#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
axs[2].legend(handles, labels,fontsize=fs-1,loc='upper right')#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)

quit()
save_axis(axs[0], nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_CFL"+str(CFL_test)+"_V_compare_diag_rho.pdf")
save_axis(axs[1], nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_CFL"+str(CFL_test)+"_V_compare_diag_mod_v.pdf")
save_axis(axs[2], nametest+"_timescheme_"+name_time_scheme[time_scheme]+"_CFL"+str(CFL_test)+"_V_compare_diag_p.pdf")
