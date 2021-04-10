"""
This script compare results (in VTK format) from SRH-2D and HEC-RAS 2D
"""

import pyHMT2D

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

    vtkFileName1 = "SRH-2D/SRH2D_Muncie_C_0006.vtk"
    vtkFileName2 = "HEC-RAS/RAS2D_Muncie_Geom_0024.vtk"

    #water depth
    varName1 = "Water_Depth_ft"
    varName2 = "Water_Depth_ft"
    varNameDiff = "water_depth_diff_cell"

    vtkFileNameDiff = "diff_depth_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff)

    #WSE
    varName1 = "Water_Elev_ft"
    varName2 = "Water_Elev_ft"
    varNameDiff = "wse_diff_cell"

    vtkFileNameDiff = "diff_wse_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff)

    #bed elevation
    varName1 = "Bed_Elev_ft"
    varName2 = "Bed_Elev_ft"
    varNameDiff = "bed_elevation_diff_cell"

    vtkFileNameDiff = "diff_bed_elevation_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff, diffNodal=True)

    #velocity
    varName1 = "Velocity_ft_p_s"
    varName2 = "Velocity_ft_p_s"
    varNameDiff = "velocity_diff_cell"

    vtkFileNameDiff = "diff_velocity_cell.vtk"

    vtk_handler.vtk_diff_consistent(vtkFileName1, vtkFileName2, vtkFileNameDiff,
                                    varName1, varName2, varNameDiff)


if __name__ == "__main__":

    calculate_2D_difference()

    print("All done!")