"""
This script compare results (in VTK format) from SRH-2D and HEC-RAS 2D
"""

import numpy as np

from matplotlib import pyplot as plt

import pyHMT2D
from pyHMT2D.Misc.tools import setNumpyArrayValueToNaN

plt.rc('text', usetex=False)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def plot_1D_profile():
    vtk_handler = pyHMT2D.Misc.vtkHandler()

    # Define the filename of VTK file
    srh_filename = 'SRH-2D/backwater_curve_N_0024.vtk'
    hec_filename = 'HEC-RAS-2D/RAS2D_channel_0012.vtk'

    # Set the points between which the line is constructed.
    p1 = [259651, 4.52109e+06, 0]
    p2 = [269601, 4.52104e+06, 0]

    # Define the numer of interpolation points
    numPoints = 100

    # for SRH-2D
    reader = vtk_handler.readVTK_UnstructuredGrid(srh_filename)  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points, U, elev_srh_2d = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader,
                                                                          'Water_Elev_m')  # interpolate the data over the line
    # print(points)
    # print(U)

    max_x = np.max(points[:, 0])
    min_x = np.min(points[:, 0])

    U = setNumpyArrayValueToNaN(U, 0)  # Set the zero's to NaN's

    plt.plot((points[:, 0] - max_x) / 1000, U, label='SRH-2D')  # plot the data

    # plot the bottom
    plt.plot((points[:, 0] - max_x) / 1000, elev_srh_2d, color='k')

    # for HEC-RAS
    reader = vtk_handler.readVTK_UnstructuredGrid(hec_filename)  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points, U, elev_ras_2d = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader,
                                                                          'WSE')  # interpolate the data over the line
    # print(points)
    # print(U)

    max_x = np.max(points[:, 0])
    min_x = np.min(points[:, 0])

    U = setNumpyArrayValueToNaN(U, 0)  # Set the zero's to NaN's

    plt.scatter(((points[:, 0] - max_x) / 1000)[::5], U[::5], marker='o', facecolors='none', edgecolors='r', s=10,
                label='HEC-RAS 2D')

    # plot the bottom
    plt.plot((points[:, 0] - max_x) / 1000, elev_ras_2d, color='k')

    # 1D backwater calculation
    slope = 1e-5
    ManningN = 0.03
    startx = 269601
    startH = 2.0
    startZ = 0.0
    riverLength = 10000
    nGrid = 1001
    specificDischarge = 0.48  # = 300 cms/625 m (the inlet discharge in HEC-RAS and SRH-2D is 300 cms and the inlet
    # width is 625 m

    normalDepth, criticalDepth, x, waterDepth, WSE = pyHMT2D.Misc.ocf_1D_backwater_curve(slope,
                                                                                         ManningN, startx, startH,
                                                                                         startZ, riverLength, nGrid,
                                                                                         specificDischarge)

    negX = -(x - startx) / 1000
    print("negX, WSE = ", negX, WSE)

    plt.scatter(negX[::50], WSE[::50], marker='^', facecolor='none', edgecolor='blue', s=10, linewidths=1,
                label='1D Backwater')
    # plt.plot(negX,WSE)

    plt.xlabel('x (km)', fontsize=16)
    plt.ylabel('Elevation (m)', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=12)

    plt.title('Comparison between SRH-2D and HEC-RAS 2D')

    plt.xlim(-11, 1)
    plt.ylim(-0.3, 4);

    plt.legend(loc='upper right', fontsize=14, frameon=False)

    plt.show()

def calculate_2D_difference():
    """Calculate the difference in 2D result; the difference for each variable is saved to a vtk file

    Pay attention to where the flow variables are stored (nodal or cell center). It is probably more convenient to
    specify "XMDFC" in SRHHYDRO so all SRH-2D results are at cell center. But it seems SMS does not provide this
    option in its GUI interface. The SRHHYDRO file can be modified manually.

    As of now, all calculations will be at cell centers. If a specified variable in a VTk file is nodal, it will be
    converted to cell center first.

    Returns
    -------

    """

    vtk_handler = pyHMT2D.Misc.vtkHandler()

    vtkFileName1 = "SRH-2D/backwater_curve_C_0023.vtk"
    vtkFileName2 = "HEC-RAS-2D/RAS2D_channel_0012.vtk"

    #water depth
    varName1 = "Water_Depth_m"
    varName2 = "Water_Depth_m"
    varNameDiff = "water_depth_diff_cell"

    vtkFileNameDiff = "diff_depth_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff)

    #WSE
    varName1 = "Water_Elev_m"
    varName2 = "Water_Elev_m"
    varNameDiff = "wse_diff_cell"

    vtkFileNameDiff = "diff_wse_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff)

    #bed elevation
    varName1 = "Bed_Elev_m"
    varName2 = "Bed_Elev_m"
    varNameDiff = "bed_elevation_diff_cell"

    vtkFileNameDiff = "diff_bed_elevation_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff)

    #velocity
    varName1 = "Velocity_m_p_s"
    varName2 = "Velocity_m_p_s"
    varNameDiff = "velocity_diff_cell"

    vtkFileNameDiff = "diff_velocity_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff)


if __name__ == "__main__":

    calculate_2D_difference()

    #plot_1D_profile()

    print("All done!")