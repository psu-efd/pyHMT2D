"""
A module with the function to call SRH-2D. This function has to be in a separate Python file for the parallelization
to work.
"""

import multiprocessing
import os
import shutil

from pathlib import Path

import pyHMT2D

def run_one_SRH_2D_case(case_ID, srhcontrol_file, ManningN_MaterialIDs, ManningNs, ManningN_MaterialNames, bDeleteCaseDir=True):
    """
    Run a single SRH-2D case with the specified case_ID. The ManningN value is used to set the Manning's n for the zones specified by ManningN_MaterialIDs.

    Parameters
    ----------
    case_ID : int
        ID of the case to run
    srhcontrol_file : str
        Name of the SRH-2D control file (either srhhydro or SIF file), e.g., "case_SIF.dat" or "case.srhhydro"
    ManningN_MaterialIDs : list of int
        IDs of the Manning's n zones. Note: In HEC-RAS, the ID is 1-based. We need -1 to make it 0-based in Python.
    ManningNs : list of float
        Manning's n values for the ManningN_MaterialIDs
    ManningN_MaterialNames : list of str
        Names of the Manning's n zones
    bDeleteCaseDir : bool
        whether to delete case directory when simulation is done. Default is True. If it is set to False, the case directory
        will be kept. However, be cautious on this. If the number of simulations is large, all cases together will use
        very large space on the hard disk.

    Returns
    -------

    """

    processID = multiprocessing.current_process()
    print(processID, ": running SRH-2D case: ", case_ID)

    #base SRH-2D case directory
    base_case_dir = 'base_case'

    #name of directory for the new case, e.g., case_000001, case_000002, etc.
    case_name = 'case_'+str(case_ID).zfill(6)
    new_case_dir = 'cases/'+case_name

    #copy the base SRH-2D case folder
    destination = shutil.copytree(base_case_dir, new_case_dir, dirs_exist_ok=True)

    #chdir to the new case directory
    os.chdir(destination)

    #create a SRH-2D model instance
    version = "3.6.5"
    srh_pre_path = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
    srh_path = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console.exe"
    extra_dll_path = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe"

    #create a SRH-2D model instance
    my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(version, srh_pre_path,
                       srh_path, extra_dll_path, faceless=False)

    #initialize the SRH-2D model
    my_srh_2d_model.init_model()

    print("Hydraulic model name: ", my_srh_2d_model.getName())
    print("Hydraulic model version: ", my_srh_2d_model.getVersion())

    #open a SRH-2D project
    my_srh_2d_model.open_project(srhcontrol_file)

    #get the SRH-2D data
    my_srh_2d_data = my_srh_2d_model.get_simulation_case()

    #modify the Manning's n if the list is not empty
    if len(ManningN_MaterialIDs) > 0:
        my_srh_2d_data.modify_ManningsNs(ManningN_MaterialIDs, ManningNs, ManningN_MaterialNames)

    #save the srhcontrol file after the updates of parameters
    my_srh_2d_data.srhhydro_obj.save_as()  #without any argument, the original filename will be used

    #run SRH-2D Pre to preprocess the case
    bRunSuccessful = my_srh_2d_model.run_pre_model()

    #run the SRH-2D model's current project if SRH-2D Preprocessing is successful
    if bRunSuccessful:
        bRunSuccessful = my_srh_2d_model.run_model()

    #close the SRH-2D project
    my_srh_2d_model.close_project()

    #quit SRH-2D
    my_srh_2d_model.exit_model()

    #do postprocessing only if bRunSuccessful is true
    if bRunSuccessful:
        #convert SRH-2D result to VTK (This is hard-coded; needs to be changed for a specific case)
        my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(srhcontrol_file)

        #wether the result is nodal or cell center
        bNodal = False

        my_srh_2d_data.readSRHXMDFFile(my_srh_2d_data.get_case_name() + "_XMDFC.h5", bNodal)

        #export the SRH-2D result to VTK: lastTimeStep=True means we only want to deal with the last time step.
        # See the code documentation of outputXMDFDataToVTK(...) for more options. It returns a list of vtk file names.
        vtkFileNameList = my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir='')

        #copy the vtk result file (only the last time step) to "cases" directory (one level above)
        shutil.copy(vtkFileNameList[-1], "../"+case_name+".vtk")

    # go back to the root
    os.chdir("../..")

    # delete the case directory
    if bDeleteCaseDir:
        shutil.rmtree(destination)

    #if successful, return case_ID; otherwise, return -case_ID
    if bRunSuccessful:
        return  case_ID
    else:
        return -case_ID
