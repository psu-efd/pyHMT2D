"""
Some tools
"""

import vtk
import sys

def read_sampling_points(num_sampling_files, \
                         sampling_file_names, \
                         output_vtp_file_name \
                         ):
    """Read sampling points from a list of sampling points files

    To be used to sample 2D hydraulic modeling results with VTK.

    Parameters
    ----------
    num_sampling_files: number of sampling points files
    sampling_file_names: list of sampling points file names
    output_vtp_file_name: output vtp file names

    Returns
    -------
    vtkPolyData
    """

    sampling_points = vtk.vtkPoints()

    sample_ID = vtk.vtkIntArray()
    sample_ID.SetName("sample_ID")

    # read in the coordinates of sample points from file (only x and y are relevant due to 2D)
    z = 0.0

    # load all sample points files and each file is assigned a sample_ID
    for i in range(0, num_sampling_files):
        with open(sampling_file_names[i]) as f:
            for line in f:
                x, y = line.split()
                # print(x,y,z)
                sampling_points.InsertNextPoint(float(x), float(y), z)
                sample_ID.InsertNextValue(i)

    # Create a polydata object
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(sampling_points)

    polydata.GetPointData().SetScalars(sample_ID)

    # print(polydata)
    # print(polydata.GetPointData())

    # write out the combined sample points into vtp format -> load into Paraview for inspection.
    unstr_writer = vtk.vtkXMLPolyDataWriter()  # this can save as vtu format
    unstr_writer.SetFileName(output_vtp_file_name)
    unstr_writer.SetInputData(polydata)
    unstr_writer.Write()

    # print("Done.")

    return polydata


def read_sampling_input_file(sampling_input_file_name):
    """Read sampling input file

    The format of the sampling input file is as follows:
    ___________________________
    #one comment line start with "#"
    number of sampling point file names
    #one comment line start with "#"
    filename_1.dat
    filename_2.dat
    ...
    filename_n.dat
    #one comment line start with "#"
    number of sample variables
    #one comment line start with "#"
    varname_1
    varname_2
    ...
    varname_m
    ___________________________

    An example is as follows:
    ___________________________

    #number of sampling point files
    4
    #list of sampling point file names
    sampling_section_1.dat
    sampling_section_2.dat
    sampling_section_3.dat
    sampling_long_profile_1.dat
    #number of sampling variables
    3
    #list of sampling variables (need to be consistent with 2D hydraulic modeling output names)
    WATER_DEPTH_FT
    VEL_X_FT_P_S
    VEL_Y_FT_P_S

    ____________________________

    Parameters
    ----------
    sampling_input_file_name: sampling input file name

    Returns
    -------
    num_sampling_files: number of sampling files
    sampling_file_names: names of the sampling files
    num_sampling_variables: number of sampling variables
    sampling_variable_names: names of the sampling variables

    """
    file_sampling_input = open(sampling_input_file_name, 'r')

    # readout one comment line
    line = file_sampling_input.readline()

    # read the number of sampling point files
    num_sampling_files = int(file_sampling_input.readline())
    # print("num_sampling_files = ", num_sampling_files)

    # readout one comment line
    line = file_sampling_input.readline()

    # read the list of sampling point file names
    sampling_file_names = []
    for i in range(0, num_sampling_files):
        sampling_file_names.append(file_sampling_input.readline().rstrip())
    # print("sampling_file_names = ", sampling_file_names)

    # readout one comment line
    line = file_sampling_input.readline()

    # read the number of sampling variables
    num_sampling_variables = int(file_sampling_input.readline())
    # print("num_sampling_variables = ", num_sampling_variables)

    # readout one comment line
    line = file_sampling_input.readline()

    # read the list of sampling variable names
    sampling_variable_names = []
    for i in range(0, num_sampling_variables):
        sampling_variable_names.append(file_sampling_input.readline().rstrip())
    # print("sampling_variable_names = ", sampling_variable_names)

    file_sampling_input.close()

    # print("Done.")

    return [num_sampling_files, sampling_file_names, num_sampling_variables, sampling_variable_names]

def sampling_on_vtu(vtuFileName, samplingPointsFileName, fieldName):
    """Sampling on vtu

    Parameters
    ----------
    vtuFileName: vtu file name
    samplingPointsFileName: sampling points file name
    fieldName: what field variable to sample

    Returns
    -------

    """

    return NotImplemented

    resultVtu = vtu(vtuFileName)

    sampling_points = vtk.vtkPoints()

    #read in the coordinates of sample points from file (only x and y are relevant due to 2D)
    z = 0.0

    with open(samplingPointsFileName) as f:
        for line in f:
            x, y = line.split()
            #print(x,y,z)
            sampling_points.InsertNextPoint(float(x),float(y),z)

    #print(sampling_points)

    probe = vtktools.VTU_Probe(resultVtu.ugrid, GetLocations(sampling_points))

    probedData = probe.GetField(fieldName)

    print(probedData)

    return probedData

    print("Done.")

def printProgressBar(i,total,postText):
    """A simple function to print a progress bar (no need to use other libraries)

    Note: any print inbetween the calling of this function will create a new progress bar. So
    it is better not print anything in between the calling.

    Reference:
    https://stackoverflow.com/questions/3002085/python-to-print-out-status-bar-and-percentage

    Parameters
    ----------
    i
    total
    postText

    Returns
    -------

    """
    n_bar =20 #size of progress bar
    j= (i+1)/total
    j=max(j,0.01001)  #at least print something at the beginning
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%  {postText}")
    sys.stdout.flush()
