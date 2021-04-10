import numpy as np
import sys
import csv

import vtk

import pyHMT2D

from pyHMT2D import gVerbose

from matplotlib import pyplot as plt

plt.rc('text', usetex=False)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def read_sampling_points_from_csv_file(fileName):
    """Read name and coordiantes (x,y) of sampling points from csv file

    The format of the csv file
    Name, x, y
    Point1, xx, xx
    Point2, xx, xx
    ...


    Returns
    -------

    """

    with open(fileName, 'r') as csvfile:
        # creating a csv reader object
        csvreader = csv.reader(csvfile)

        # extracting field names in first row
        fields = next(csvreader)

        # only 3 columns allowed
        if len(fields) != 3:
            raise Exception("The number of columns in the sampling point file needs to be 3. Exit ...")

        rows = []

        # loop and extracting each data row
        for row in csvreader:
            rows.append(row)

        # get total number of sampling points
        if gVerbose: print("Total number of sampling points: %d" % (csvreader.line_num - 1))

    # printing all rows
    if gVerbose: print('\nSampling points are:\n')
    for row in rows:
        if gVerbose: print(', '.join(col for col in row))

    pointNames = []

    # build the measurement points numpy array
    samplingPoints = np.zeros((len(rows), 2))

    for pointI in range(len(rows)):
        pointNames.append(rows[pointI][0])
        samplingPoints[pointI, 0] = float(rows[pointI][1])
        samplingPoints[pointI, 1] = float(rows[pointI][2])

    return pointNames, samplingPoints

def demo_hec_ras_model_data():

    #set HEC-RAS
    version = "5.0.7"

    #whether to run HEC-RAS faceless
    faceless = False

    #create a HEC-RAS model instance
    my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version, faceless)

    #initialize the HEC-RAS model
    my_hec_ras_model.init_model()

    print("Hydraulic model name: ", my_hec_ras_model.getName())
    print("Hydraulic model version: ", my_hec_ras_model.getVersion())

    #open a HEC-RAS project
    my_hec_ras_model.open_project("Muncie2D.prj", terrainFileName = "Terrain/TerrainMuncie_composite.tif")

    #run the HEC-RAS model's current project
    #my_hec_ras_model.run_model()

    #read the HEC-RAS simulation result

    #get the RAS_2D_Data object in the HEC-RAS model
    my_ras_2d_data = my_hec_ras_model.get_simulation_case()

    #read in the result
    my_ras_2d_data.load2DAreaSolutions()

    #save the HEC-RAS simulation result to VTK. It returns a list of VTK files
    vtkFileNameList = my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)

    #sample some results as manufactured solution

    #load the sampling points from file
    pointNames, samplingPoints_temp = read_sampling_points_from_csv_file("sample_points.dat")

    # Set the sampling points as vtkPoints
    sampling_points = vtk.vtkPoints()

    for pointI in range(samplingPoints_temp.shape[0]):
        sampling_points.InsertNextPoint(samplingPoints_temp[pointI,0],
                                        samplingPoints_temp[pointI,1],
                                        0.0)

    vtk_handler = pyHMT2D.Misc.vtkHandler()

    #we only use the last time step result vtk file
    vtkUnstructuredGridReader = vtk_handler.readVTK_UnstructuredGrid(vtkFileNameList[-1])

    # sample WSE on the sampling points
    points, sampled_WSE, bed_elev = vtk_handler.probeUnstructuredGridVTKOverLine(
        sampling_points, vtkUnstructuredGridReader,
        "Water_Elev_ft")

    # sample velocity on the sampling points
    points, sampled_U, bed_elev = vtk_handler.probeUnstructuredGridVTKOverLine(
        sampling_points, vtkUnstructuredGridReader,
        "Velocity_ft_p_s")

    #output the true solution and sampled WSE and velocity
    fid = open("sampled_wse.csv", "w")
    fid.write("Name,x,y,wse\n")
    for pointI in range(len(samplingPoints_temp)):
       fid.write("%s, %f, %f, %f\n" % (pointNames[pointI],samplingPoints_temp[pointI,0],
                                        samplingPoints_temp[pointI,1],sampled_WSE[pointI]))
    fid.close()

    fid = open("sampled_velocity.csv", "w")
    fid.write("Name,x,y,angle,velocity\n")
    for pointI in range(len(samplingPoints_temp)):
        # azimuth angle
        angle = np.degrees(np.arctan(sampled_U[pointI,0] / (sampled_U[pointI,1] + 1e-8)))
        # velocity magnitude
        mag = np.sqrt(np.power(sampled_U[pointI,0], 2) + np.power(sampled_U[pointI,1], 2))

        fid.write("%s, %f, %f, %f, %f\n" % (pointNames[pointI],samplingPoints_temp[pointI,0],\
                                           samplingPoints_temp[pointI,1],angle,mag))
    fid.close()


    #close the HEC-RAS project
    my_hec_ras_model.close_project()

    #quit HEC-RAS
    my_hec_ras_model.exit_model()


def demo_calibrator():
    """ Test the Calibrator class

    Returns
    -------

    """

    #my_calibrator = pyHMT2D.Calibration.Calibrator("calibration.json")
    my_calibrator = pyHMT2D.Calibration.Calibrator("calibration_brute_force.json")
    #my_calibrator = pyHMT2D.Calibration.Calibrator("calibration_Nelder-Mead.json")
    #my_calibrator = pyHMT2D.Calibration.Calibrator("calibration_L-BFGS-B.json")

    my_calibrator.calibrate()

    my_calibrator.close_and_cleanup()


if __name__ == "__main__":

    #run the HEC-RAS model and sample some data for calibration (manufactured solution)
    #demo_hec_ras_model_data()

    demo_calibrator()

    print("All done!")