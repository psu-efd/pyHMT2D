"""
Perform Monte Carlo simulation with the HEC_RAS_Model and HEC_RAS_Data classes

It demonstrates how to use pyHMT2D to control the run of HEC-RAS model, extract information, and parallelization.
"""

#if run in the cloud, need to specify the location of pyHMT2D. If pyHMT2D is installed
#with regular pip install, then the following is not necessary.
#import sys
#sys.path.append(r"C:\Users\Administrator\python_packages\pyHMT2D")

import pyHMT2D

import numpy as np
from scipy.stats import truncnorm
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

import vtk
from vtk.numpy_interface import dataset_adapter as dsa

import multiprocessing
import os
import json

from HEC_RAS_solver_module import run_one_HEC_RAS_case

plt.rc('text', usetex=False)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def Monte_Carlo_simulations_parallel(samples, ManningN_MaterialID, ManningN_MaterialName, nProcess, bDeleteCaseDir=True):
    """
    Run Monte Carlo simulations nSamples times. Run all cases in parallel using multiple CPU cores. Each core
    simulates one case at a time.

    In HEC-RAS case setup, we should use single core ("Unsteady Flow Analysis"->"Options"->
    "Computation Options and Tolerances ..."->"General" and "2D Flow Options", select "1 cores").
    It does not make sense to use all available cores because we use these cores to run multiple cases in parallel.

    In this demo example, only one Manning's n is varied for a zone in the domain.

    Parameters
    ----------
    samples : numpy array
        list of Manning's n values
    ManningN_MaterialID : int
        ID of the Manning's n zone. Note: in HEC-RAS, the ID is 1-based.
    ManningN_MaterialName : str
        Name of the Manning's n zone
    nProcess : int
        number of parallel processes to use
    bDeleteCaseDir : bool
        whether to delete case directory when simulation is done. Default is True. If it is set to False, the case directory
        will be kept. However, be cautious on this. If the number of simulations is large, all cases together will use
        very large space on the hard disk.

        In this example, the HEC-RAS result is saved into a VTK file (thus the case directory is not needed after
        simulation is done).

    Returns
    -------

    """

    print("Running Monte Carlo simulations in parallel with", nProcess, "processes ...")

    #number of samples
    nSamples = samples.shape[0]

    #create a pool of parallel workers
    pool = multiprocessing.Pool(processes=nProcess)

    #make the directory "cases". If it already exists, reminder user to remove it first and exit.
    if os.path.isdir("cases"):
        print("The directory cases already exists. Make sure to remove it before run this script.")
        quit()
    else:
        os.mkdir("cases")

    #create the list of case IDs to be assigned to the parallel workers
    caseIDs = [*range(0, nSamples, 1)]
    #caseIDs = [*range(0, 4, 1)]

    #prepare the arguments for calling the "run_one_HEC_RAS_case(...)" function
    args = [(caseID, ManningN_MaterialID, ManningN, ManningN_MaterialName, bDeleteCaseDir) for caseID, ManningN in zip(caseIDs, samples)]

    # do the work
    outputs = pool.starmap(run_one_HEC_RAS_case, args)

    # save the returned outputs as a json file: outputs contains the caseIDs (positive if the case run was successful,
    # negative otherwise)
    with open('parallel_run_RAS_2D_outputs.json', 'w') as f:
        json.dump(outputs, f)

    print("Parallel run of cases done. Outputs = ", outputs)

def Monte_Carlo_simulations_serial(samples, ManningN_MaterialID, ManningN_MaterialName, bDeleteCaseDir=True):
    """
    Run Monte Carlo simulations in serial

    Parameters
    ----------
    samples : numpy array
        list of Manning's n values
    ManningN_MaterialID : int
        ID of the Manning's n zone. Note: in HEC-RAS, the ID is 1-based.
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

    print("Running Monte Carlo simulations in serial ...")

    #number of samples
    nSamples = samples.shape[0]

    #make the directory "cases". If it already exists, reminder user to remove it first and exit.
    if os.path.isdir("cases"):
        print("The directory cases already exists. Make sure to remove it before run this script.")
        quit()
    else:
        os.mkdir("cases")

    #list to collect outputs
    outputs = []

    for i in range(2):  #nSamples
        print("Running Monte Carlo simulation # ", i, "out of", nSamples)

        ManningN = samples[i]

        output = run_one_HEC_RAS_case(i, ManningN_MaterialID, ManningN, ManningN_MaterialName, bDeleteCaseDir)

        outputs.append(output)

    # save the returned outputs as a json file: outputs contains the caseIDs (positive if the case run was successful,
    # negative otherwise)
    with open('serial_run_RAS_2D_outputs.json', 'w') as f:
        json.dump(outputs, f)

    print("Serial run of cases done. Outputs = ", outputs)

def postprocess_results(nSamples):
    """
    Postprocess the Monte Carlo simulation results, which are in the format of VTK files in "cases" directory.

    - Statistics at a given a point
    - Innundation map with certain confidence intervals

    Returns
    -------

    """

    #probe point location (hard-coded; need to change for a specific case)
    probe_x = 411100.0
    probe_y = 1.80345e+06

    #array to store the probed WSE at the probing point
    probed_WSE = np.zeros(nSamples)

    #array to store water depth field from all VTK files
    water_depth_all = None

    vtk_handler = pyHMT2D.Misc.vtkHandler()

    #loop over all VTK result files
    for i in range(nSamples):
        #VTK file name
        vtk_filename = "cases/case_"+str(i).zfill(6)+".vtk"

        print("Processing ", vtk_filename)

        #create VTK unstructured grid reader
        reader = vtk_handler.readVTK_UnstructuredGrid(vtk_filename)

        #create the probe point in vtkPoints format
        pointVTK = vtk.vtkPoints()

        #There is only one point to probe at
        pointVTK.SetNumberOfPoints(1)

        #set the coordinates of the probing point
        pointVTK.InsertPoint(0, probe_x, probe_y, 0.0)

        #probe WSE at the probing point
        points, WSE, elev = vtk_handler.probeUnstructuredGridVTKOnPoints(pointVTK, reader, 'Water_Elev_ft')

        print("WSE = ", WSE)

        probed_WSE[i] = WSE[0]

        #get the water depth at cell centers from the VTK
        water_depth = vtk_handler.get_uGRid_cell_field_with_name(reader, 'Water_Depth_ft')

        if water_depth_all is None:
            water_depth_all = np.zeros((water_depth.shape[0], nSamples))

        water_depth_all[:,i] = water_depth

    #calculate the WSE at some specific exceedance probabilities
    exceedance_probabilities = [99, 90, 50, 10, 1]
    indices = [min(int(pi *(nSamples+1)/100), nSamples-1) for pi in exceedance_probabilities]

    print("probed_WSE = ", probed_WSE)
    print("indices = ", indices)

    probe_WSE_sorted = np.sort(probed_WSE)[::-1]

    #print out the WSE values corresponding to the exceedance probabilities
    WSE_exceedance = []

    for i in indices:
        WSE_exceedance.append(probe_WSE_sorted[i])

    print("WSE_exceedance = ", WSE_exceedance)

    #plot WSE exceedance of probability at monitoring point
    plot_WSE_exceedance_probability(exceedance_probabilities, probe_WSE_sorted)

    #process the water depth at each cell center in the mesh
    water_depth_exceedance_of_probabilities = np.zeros((water_depth_all.shape[0], len(exceedance_probabilities)))

    #loop over all cells
    for iCell in range(water_depth_all.shape[0]):
        depths_at_cell = water_depth_all[iCell,:]

        dephts_at_cell_sorted = np.sort(depths_at_cell)[::-1]

        for i, indx in enumerate(indices):
            water_depth_exceedance_of_probabilities[iCell, i] = dephts_at_cell_sorted[indx]

    #save the cell center water depth exceedance of probability to VTK file
    # Read in base vtk (just use "cases/case_000000.vtk")
    fileName = "cases/case_000000.vtk"
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    mesh = reader.GetOutput()

    # Add data set and write VTK file
    meshNew = dsa.WrapDataObject(mesh)

    for i, exceedance_probability in enumerate(exceedance_probabilities):
        meshNew.CellData.append(water_depth_exceedance_of_probabilities[:,i], "h_"+str(exceedance_probability)+"_percent")

    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName("water_depth_exceedance_of_probabilities.vtk")
    writer.SetInputData(meshNew.VTKObject)
    writer.Write()

    print("Done.")

def plot_WSE_exceedance_probability(exceedance_probabilities, probe_WSE_sorted):
    """
    Plot the WSE exceedance of probability

    Parameters
    ----------
    exceedance_probabilities
    probe_WSE_sorted

    Returns
    -------

    """

    exceedance = np.arange(1.0, len(probe_WSE_sorted)+1) / len(probe_WSE_sorted)

    plt.plot(exceedance*100, probe_WSE_sorted)

    # set the limit for the x and y axes
    plt.ylim([0, 100])
    plt.ylim([min(probe_WSE_sorted), max(probe_WSE_sorted)])

    # set x and y axes label and font size
    plt.xlabel(r'Probability of Exceedence (\%)', fontsize=16)
    plt.ylabel(r'WSE (ft)', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=12)

    # set axis label format
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))

    # show title and set font size
    plt.title('WSE exceedance of probability at monitoring point', fontsize=16)

    #save the plot to a file
    plt.savefig("WSE_exceedance_of_probability.png", dpi=300, bbox_inches='tight')

    # show legend, set its location, font size, and turn off the frame
    # plt.legend(loc='lower left', fontsize=14, frameon=False)
    #plt.show()


if __name__ == "__main__":   

    #or load saved Manning's n values from file
    samples = np.loadtxt("sampledManningN_2022_03_10-10_24_34_PM.dat", delimiter=',')

    # number of samples
    nSamples = samples.shape[0]

    #ID and name of the Manning's n zone whose n value is changes in the Monte Carlo simulations.
    #The ID and name are hard-coded here. Need to modify for a specific case.
    ManningN_MaterialID = 2           #Note: In HEC-RAS, the ID is 1-based.
    ManningN_MaterialName = 'channel'

    #run Monte Carlo simulations in parallel: nProcess is the number of cores to use.
    nProcess = 4   #number of cores to use
    #Monte_Carlo_simulations_parallel(samples, ManningN_MaterialID, ManningN_MaterialName, nProcess=nProcess, bDeleteCaseDir=True)

    #or run Monte Carlo simulations in serial
    #Monte_Carlo_simulations_serial(samples, ManningN_MaterialID, ManningN_MaterialName, bDeleteCaseDir=True)

    #postprocessing the results
    postprocess_results(nSamples)

    print("All done!")

