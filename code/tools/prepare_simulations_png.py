import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

variable_name = "W_X"

parser = OptionParser()
parser.add_option("-d", "--dir", action="store", type="string", dest="d", default="tsunami")
parser.add_option("-q", "--quality", action="store", type="int", dest="quality", default=200, 
                  help="Set the quality of the output images (DPI). Default is 200.")
options, args = parser.parse_args()

folder = options.d
output_format = "png"  # Save PNG images
quality = options.quality

aspect_ratio = "quad"

if aspect_ratio == "quad":
    aspect_ratio_dimensions = [4, 4]
elif aspect_ratio == "rect_X":
    aspect_ratio_dimensions = [14, 4]
else: 
    print("Aspect ratio undefined")
    sys.exit(1)

gmm1 = 1.249
gmm2 = 1.4
pinf1 = 0.0
pinf2 = 0.0

print(folder)
if not os.path.isdir(folder):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

def plotting(folder_name, levels, output_folder, quality):
    folder_name = folder_name + "/"
    
    delimiter_in = ' '
    headerlines_in = 1

    # Importing the MESH DATA
    import_mesh = [file for file in os.listdir(folder_name) if file.startswith(variable_name + '_MESH_')]
    import_mesh.sort()
    nfiles_mesh = len(import_mesh)
    for i in range(nfiles_mesh):
        filename = import_mesh[i]
        mydata_mesh = np.loadtxt(folder_name + filename, delimiter=delimiter_in, skiprows=headerlines_in)
        if i == 0:
            x = mydata_mesh
            x_ini = min(x)
            x_end = max(x)
        elif i == 1:
            y = mydata_mesh
            y_ini = min(y)
            y_end = max(y)

    X, Y = np.meshgrid(x, y)

    # Importing the SOLUTION DATA
    import_solution = [file for file in os.listdir(folder_name) if file.startswith(variable_name + '_SOLUTION_')]
    import_solution.sort()
    nfiles_solution = len(import_solution)

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    fig, axv = plt.subplots(figsize=(8, 6))

    for i in range(nfiles_solution):
        filename = import_solution[i]
        mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)
        ro = mydata_solution[:, 0]
        rou = mydata_solution[:, 1]
        rov = mydata_solution[:, 2]
        energy = mydata_solution[:, 3]
        levelset = mydata_solution[:, 4]
        phi = mydata_solution[:, 5]

        RO = np.reshape(ro, (X.shape[0], Y.shape[1]))
        ROU = np.reshape(rou, (X.shape[0], Y.shape[1]))
        ROV = np.reshape(rov, (X.shape[0], Y.shape[1]))
        ENERGY = np.reshape(energy, (X.shape[0], Y.shape[1]))
        LEVELSET = np.reshape(levelset, (X.shape[0], Y.shape[1]))
        PHI = np.reshape(phi, (X.shape[0], Y.shape[1]))

        U = ROU / RO
        V = ROV / RO

        gmm = np.where(LEVELSET >= 0, gmm1, gmm2)
        pinf = np.where(LEVELSET >= 0, pinf1, pinf2)

        P = (gmm - 1.0) * (ENERGY - 0.5 * RO * (U**2 + V**2)) - gmm * pinf

        VORT = np.zeros_like(RO)
        NX = X.shape[0]
        NY = Y.shape[1]

        dx = X[1, 1] - X[1, 0]
        dy = Y[1, 1] - Y[0, 1]

        for indi in range(1, NX - 1):
            for indj in range(1, NY - 1):
                VORT[indi, indj] = np.sqrt(((RO[indi + 1, indj] - RO[indi - 1, indj]) / (2 * dx))**2 + ((RO[indi, indj + 1] - RO[indi, indj - 1]) / (2 * dy))**2)

        maxgrad = np.max(VORT[:,:])
        K = 80
        VORT = np.exp(-K * VORT / maxgrad)

        # Clear the axis and update the plot with new data
        axv.clear()
        cont_vort = axv.contourf(X, Y, VORT, cmap='Greys_r', levels=levels)

        axv.set_xlabel('X')
        axv.set_ylabel('Y')
        axv.set_aspect(aspect_ratio_dimensions[0] / aspect_ratio_dimensions[1])  # Set aspect ratio

        # Save the frame as a PNG file with adjustable quality
        output_filename = os.path.join(output_folder, f"frame_{i + 1}.png")  # Save as frame_1, frame_2, etc.
        plt.savefig(output_filename, dpi=quality, bbox_inches='tight')

if __name__ == '__main__':
    output_folder = "simulation_frames"
    plotting(folder, 30, output_folder, quality)
