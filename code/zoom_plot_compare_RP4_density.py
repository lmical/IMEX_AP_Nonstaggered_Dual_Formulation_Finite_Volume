#####################
#CFL 0.95
#100 elements
#####################
#NB: Results with reconstruction of conserved variables are too oscillatory to deserve to be zoomed
#Using the same zoom would be too much.
#####################
nametest="RP4"
file="SOLUTION_OCT_0000000.0350000.dat"
CFL_test=0.95 #0.5 #0.95 #OK
solution_type="exact solution"



gmm=1.4


teststocompare=[] 
#Reconstructed variable, riemann_solver, speed_estimate, order, CFL, RelaxedCFL, NRelaxedTimeSteps
teststocompare.append([-25,12, 1, CFL_test]) 
teststocompare.append([-24,12, 1, CFL_test]) 
teststocompare.append([-23,12, 1, CFL_test]) 
teststocompare.append([-22,12, 1, CFL_test]) 
teststocompare.append([-21,12, 1, CFL_test]) 
teststocompare.append([-20,12, 1, CFL_test]) 
teststocompare.append([-1,  12, 1, CFL_test]) 
teststocompare.append([25,  12, 1, CFL_test]) 
teststocompare.append([24,  12, 1, CFL_test]) 
teststocompare.append([23,  12, 1, CFL_test]) 
teststocompare.append([22,  12, 1, CFL_test]) 
teststocompare.append([21,  12, 1, CFL_test]) 
teststocompare.append([20,  12, 1, CFL_test]) 
teststocompare.append([1,   12, 1, CFL_test]) 



import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


linestyles_reconstruction = {0:"--",1:"-"}
markers_riemann_solvers   = {-2:"p",-1:"v",0:"*",1:"s",2:"+",3:"o",4:"^",5:"x",6:"d",7:"1"}
fs=12
lw=1.0
lwex=2.0
ms=2.0

####################################################
name_reconstructed_variable = {0:"cons",    1:"char"}
name_riemann_solver         = {0:"rusanov", 1:"exact"}
name_speed_estimate         = {0:"standard",1:"riemann"}
order_of_plotting           = {3:7,5:6,7:5,9:4,11:3,13:2}
####################################################
nVar=3

####################################################
label_reconstructed_variable = {0:"cons.",    1:"char."}
label_riemann_solver         = {0:"Rusanov", 1:"exact RS"}
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
                                  -25:"AF-o2-VA-MUSCL" }
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
                                  -1: "o",
                                  -10:"o",
                                  -2: "o",
                                  -20:"o",
                                  -21:"o",
                                  -22:"o",
                                  -23:"o",
                                  -24:"o",
                                  -25:"o"}


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


params = {'mathtext.default': 'regular' } #parameters plot




fig, ax = pl.subplots(1)
pl.xlabel("x",fontsize=fs)
pl.title("density",fontsize=fs)
pl.grid()
# Place the legend next to the main plot

# Create an inset axes object with specific dimensions and position it outside the main plot area
axins1 = inset_axes(ax, width="80%", height="80%", loc='upper right', bbox_to_anchor=(-0.4, 0.48, 0.4, 0.4), bbox_transform=fig.transFigure)
# Set the limits for the zoomed-in region
x1, x2, y1, y2 = 0.4, 0.73, 12.5, 16  # Specify the limits of the zoomed-in area
axins1.set_xlim(x1, x2)
axins1.set_ylim(y1, y2)
axins1.set_title("A",fontsize=fs)
# Add a rectangle to highlight the zoomed area in the main plot
ax.indicate_inset_zoom(axins1, edgecolor="black")
mark_inset(ax, axins1, loc1=2, loc2=4, fc="none", ec="0.5")

# Create an inset axes object with specific dimensions and position it outside the main plot area
axins2 = inset_axes(ax, width="80%", height="80%", loc='upper right', bbox_to_anchor=(-0.4, 0.05, 0.4, 0.4), bbox_transform=fig.transFigure)
# Set the limits for the zoomed-in region
x1, x2, y1, y2 = 0.405, 0.45, 3.5, 12.  # Specify the limits of the zoomed-in area
axins2.set_xlim(x1, x2)
axins2.set_ylim(y1, y2)
axins2.set_title("B",fontsize=fs)
# Add a rectangle to highlight the zoomed area in the main plot
ax.indicate_inset_zoom(axins2, edgecolor="black")
mark_inset(ax, axins2, loc1=2, loc2=4, fc="none", ec="0.5")

# Create an inset axes object with specific dimensions and position it outside the main plot area
axins3 = inset_axes(ax, width="80%", height="80%", loc='upper right', bbox_to_anchor=(0.95, 0.05, 0.4, 0.4), bbox_transform=fig.transFigure)
# Set the limits for the zoomed-in region
x1, x2, y1, y2 = 0.82, 0.86, 5.5, 12.0  # Specify the limits of the zoomed-in area
axins3.set_xlim(x1, x2)
axins3.set_ylim(y1, y2)
axins3.set_title("C",fontsize=fs)
# Add a rectangle to highlight the zoomed area in the main plot
ax.indicate_inset_zoom(axins3, edgecolor="black")
mark_inset(ax, axins3, loc1=2, loc2=4, fc="none", ec="0.5")

# Create an inset axes object with specific dimensions and position it outside the main plot area
axins4 = inset_axes(ax, width="80%", height="80%", loc='upper right', bbox_to_anchor=(0.95, 0.48, 0.4, 0.4), bbox_transform=fig.transFigure)
# Set the limits for the zoomed-in region
x1, x2, y1, y2 = 0.69, 0.84, 28.0, 34.5  # Specify the limits of the zoomed-in area
axins4.set_xlim(x1, x2)
axins4.set_ylim(y1, y2)
axins4.set_title("D",fontsize=fs)
# Add a rectangle to highlight the zoomed area in the main plot
ax.indicate_inset_zoom(axins4, edgecolor="black")
mark_inset(ax, axins4, loc1=2, loc2=4, fc="none", ec="0.5")



firsttime=True
index_plot=0

##############################
# READ AND PLOT EXACT SOLUTION
##############################
# FILLING x_ref AND ref_sol
##############################
solfoldName="reference_solutions/"+nametest
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

        ax.plot(x_ref,ro_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
        # ax_v.plot( x_ref, v_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
        # ax_p.plot( x_ref, p_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
        # ax_q.plot( x_ref, q_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
        # ax_E.plot( x_ref, E_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")

        # Plot the same data on the inset axes
        axins1.plot(x_ref,ro_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
        axins2.plot(x_ref,ro_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
        axins3.plot(x_ref,ro_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
        axins4.plot(x_ref,ro_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")


        index_plot=index_plot+1


name_base_folder="advanced_tests"
coordinate="X"

for indt, test in enumerate(teststocompare): #Loop on the schemes
    space_reconstruction=test[0]
    time_scheme         =test[1]
    post_processing     =test[2]
    CFL                 =test[3]

    foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/PP"+str(post_processing)+"/CFL"+str(CFL)

    if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
        namefile=foldName+'/'+file
        if os.path.isfile(namefile): #If the file exists, read it
            print("You are in the folder "+foldName+" dealing with "+file)
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


                    ax.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction]+" "+label_time_scheme[time_scheme],color=colors[space_reconstruction],alpha=0.7,zorder=3)
                    index_plot=index_plot+1

                    # Plot the same data on the inset axes
                    axins1.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction]+" "+label_time_scheme[time_scheme],color=colors[space_reconstruction],alpha=0.7,zorder=3)
                    axins2.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction]+" "+label_time_scheme[time_scheme],color=colors[space_reconstruction],alpha=0.7,zorder=3)
                    axins3.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction]+" "+label_time_scheme[time_scheme],color=colors[space_reconstruction],alpha=0.7,zorder=3)
                    axins4.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction]+" "+label_time_scheme[time_scheme],color=colors[space_reconstruction],alpha=0.7,zorder=3)


axins1.grid()
axins2.grid()
axins3.grid()
axins4.grid()

# ax_v.grid()
# ax_p.grid()
# ax_q.grid()
# ax_E.grid()


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
pl.legend(handles, labels,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)

pl.savefig(nametest+"_density_zoom_compare_plot_CFL"+str(CFL)+".pdf", format="pdf", bbox_inches="tight")
# pl.show()
