import numpy as np
import csv
import vtk
from matplotlib import pyplot as plt

import pyHMT2D

from pyHMT2D import gVerbose

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

def demo_srh_2d_model_data():

    #set SRH-2D
    version = "3.3"
    srh_pre_path = r"C:\Program Files\SMS 13.1 64-bit\Python36\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
    srh_path = r"C:\Program Files\SMS 13.1 64-bit\Python36\Lib\site-packages\srh2d_exe\SRH-2D_V330_Console.exe"
    extra_dll_path = r"C:\Program Files\SMS 13.1 64-bit\Python36\Lib\site-packages\srh2d_exe"

    #create a SRH-2D model instance
    my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(version, srh_pre_path,
                       srh_path, extra_dll_path, faceless=False)

    #initialize the SRH-2D model
    my_srh_2d_model.init_model()

    print("Hydraulic model name: ", my_srh_2d_model.getName())
    print("Hydraulic model version: ", my_srh_2d_model.getVersion())

    #open a SRH-2D project
    my_srh_2d_model.open_project("Muncie.srhhydro")

    #run SRH-2D Pre to preprocess the case
    my_srh_2d_model.run_pre_model()

    #run the SRH-2D model's current project
    my_srh_2d_model.run_model()

    my_srh_2d_data = my_srh_2d_model.get_simulation_case()

    # read SRH-2D result in XMDF format (*.h5)
    # Whether the XMDF result is nodal or cell center. In SRH-2D's ".srhhydro" file,
    # the output option for "OutputFormat" can be manually changed before simulation.
    # Options are "XMDF" (results at at nodes), "XMDFC" (results are at cell centers), etc.
    # For example, "OutputFormat XMDFC EN". The following lines show that the SRH-2D simulation
    # was run with "XMDFC" as output format (see the "XMDFC" part of the result file name) and thus
    # we set "bNodal = False".
    bNodal = False

    my_srh_2d_data.readSRHXMDFFile(my_srh_2d_data.get_case_name() + "_XMDFC.h5", bNodal)

    # export the SRH-2D result to VTK: lastTimeStep=True means we only want to deal with the last time step.
    # See the code documentation of outputXMDFDataToVTK(...) for more options. It returns a list of vtk file names.
    vtkFileNameList = my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir='')

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

    #close the SRH-2D project
    my_srh_2d_model.close_project()

    #quit SRH-2D
    my_srh_2d_model.exit_model()


def demo_calibrator():
    """ Demo the Calibrator class

    Returns
    -------

    """

    my_calibrator = pyHMT2D.Calibration.Calibrator("calibration.json")

    my_calibrator.calibrate()


if __name__ == "__main__":

    #run the SRH-2D model and sample some data for calibration (manufactured solution)
    #demo_srh_2d_model_data()

    #run the calibration
    demo_calibrator()

    print("All done!")