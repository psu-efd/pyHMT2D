"""
A module with the function to call HEC-RAS. This function has to be in a separate Python file for the parallelization
to work.
"""

import multiprocessing
import os
import shutil

from pathlib import Path

import pyHMT2D

def run_one_HEC_RAS_case(case_ID, ManningN_MaterialID, ManningN, ManningN_MaterialName, bDeleteCaseDir=True):
    """
    Run a sinlge HEC-RAS case with the specified case_ID. The ManningN value is used to set the Manning's n for the channel.

    Parameters
    ----------
    case_ID : int
        ID of the case to run
    ManningN_MaterialID : int
        ID of the Manning's n zone. Note: In HEC-RAS, the ID is 1-based. We need -1 to make it 0-based in Python.
    ManningN : float
        Manning's n value for this case
    ManningN_MaterialName : str
        Name of the Manning's n zone
    bDeleteCaseDir : bool
        whether to delete case directory when simulation is done. Default is True. If it is set to False, the case directory
        will be kept. However, be cautious on this. If the number of simulations is large, all cases together will use
        very large space on the hard disk.

        In this example, the HEC-RAS result is saved into a VTK file (thus the case directory is not needed after
        simulation is done).

    Returns
    -------

    """

    processID = multiprocessing.current_process()
    print(processID, ": running HEC-RAS case: ", case_ID)

    #base HEC-RAS case directory
    base_case_dir = 'base_case'

    #name of directory for the new case, e.g., case_000001, case_000002, etc.
    case_name = 'case_'+str(case_ID).zfill(6)
    new_case_dir = 'cases/'+case_name

    #copy the base HEC-RAS case folder
    destination = shutil.copytree(base_case_dir, new_case_dir)

    #chdir to the new case directory
    os.chdir(destination)

    #create a HEC-RAS model instance
    #my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="5.0.7",faceless=False)
    my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="6.1.0", faceless=False)

    #initialize the HEC-RAS model
    my_hec_ras_model.init_model()

    #print("Hydraulic model name: ", my_hec_ras_model.getName())
    print("Hydraulic model version: ", my_hec_ras_model.getVersion())

    #open a HEC-RAS project (This is hard-code; needs to be changed for a specific case)
    my_hec_ras_model.open_project("Muncie2D.prj", "Terrain/TerrainMuncie_composite.tif")

    #get the HEC-RAS case data
    my_hec_ras_data = my_hec_ras_model.get_simulation_case()

    #modify the Manning's n. In this example, we only modify one zone: channel.
    ManningN_MaterialIDs = [ManningN_MaterialID - 1]    #-1 to make it 0-based for Python
    ManningNs = [ManningN]
    ManningN_MaterialNames = [ManningN_MaterialName]
    my_hec_ras_data.modify_ManningsN(ManningN_MaterialIDs, ManningNs, ManningN_MaterialNames)

    # update the time stamp of the Manning's n GeoTiff file (to force HEC-RAS to re-compute 2D flow area's
    # properties table. (No need? The above Manning's modification already updated the time stamp.)
    if os.path.dirname(my_hec_ras_data.hdf_filename) == '':
        fileBase = b''
    else:
        fileBase = str.encode(os.path.dirname(my_hec_ras_data.hdf_filename) + '/')

    full_landcover_filename = (fileBase + my_hec_ras_data.landcover_filename).decode("ASCII")

    Path(full_landcover_filename).touch()

    # save the current project before run it
    my_hec_ras_model.save_project()

    #run the HEC-RAS model's current project
    bRunSucessful = my_hec_ras_model.run_model()

    #close the HEC-RAS project
    my_hec_ras_model.close_project()

    #quit HEC-RAS
    my_hec_ras_model.exit_model()

    #convert HEC-RAS result to VTK (This is hard-coded; needs to be changes for a specific case)
    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("Muncie2D.p01.hdf",
                                                "Terrain/TerrainMuncie_composite.tif")

    vtkFileNameList = my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)

    #copy the vtk result file (only the last time step) to "cases" directory (one level above)
    shutil.copy(vtkFileNameList[-1], "../"+case_name+".vtk")

    # go back to the root
    os.chdir("../..")

    # delete the case folder
    if bDeleteCaseDir:
        shutil.rmtree(destination)

    #if successful, return case_ID; otherwise, return -case_ID
    if bRunSucessful:
        return  case_ID
    else:
        return -case_ID