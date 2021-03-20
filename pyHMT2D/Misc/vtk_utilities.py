"""
Some utility tools to handle vtk

References:
https://stackoverflow.com/questions/21630987/probing-sampling-interpolating-vtk-data-using-python-tvtk-or-mayavi
https://vtk.org/Wiki/Demystifying_the_vtkProbeFilter

https://stackoverflow.com/questions/21630987/probing-sampling-interpolating-vtk-data-using-python-tvtk-or-mayavi
"""

import vtk
from vtk.util import numpy_support as VN

import sys

import numpy as np

class vtkHandler:
    """

    """
    def __init__(self, name=''):
        self.name = name   # name of the handler

    def readVTK_UnstructuredGrid(self, vtkFilename):
        """ read the vtk file with an unstructured grid

        Attributes
        -------
        vtkFilename: the name of the VTK file

        Returns
        -------
        reader: a vtk'vtkUnstructuredGridReader with all vectors and scalars read in

        """


        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtkFilename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        return reader

    def writeVTK_UnstructuredGrid(self, uGrid, vtkFileName):
        """ write unstructured grid to file

        Parameters
        ----------
        uGrid: vtkUnstructuredGrid object
        vtkFileName: file name to write to

        Returns
        -------

        """

        unstr_writer = vtk.vtkUnstructuredGridWriter()  # this can save as vtk format
        unstr_writer.SetFileName(vtkFileName)
        unstr_writer.SetInputData(uGrid)
        unstr_writer.Write()

    def read_sampling_points(self, num_sampling_files, \
                         sampling_file_names, \
                         output_vtp_file_name \
                         ):
        """Read sampling points from a list of sampling points files

        To be used to sample 2D hydraulic modeling results with VTK. The points are also written out to a vtp file
        for inspection in Paraview.

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


    def read_sampling_input_file(self, sampling_input_file_name):
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

        # read out one comment line
        line = file_sampling_input.readline()

        # read the number of sampling point files
        num_sampling_files = int(file_sampling_input.readline())
        # print("num_sampling_files = ", num_sampling_files)

        # read out one comment line
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

    def sampling_on_vtu(self, vtuFileName, samplingPointsFileName, fieldName):
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

    def createVtkLine(self, p1, p2, numPoints):
        # Create the line along which you want to sample along
        line = vtk.vtkLineSource()
        line.SetResolution(numPoints)
        line.SetPoint1(p1)
        line.SetPoint2(p2)
        line.Update()
        return line

    def probeUnstructuredGridVTKOverLine(self, lineVTK, readerUnstructuredGridVTK, varName):
        """ Interpolate the data from the Unstructured Grid VTK onto a line (profile).

        The unstructured grid VTK is supposed to be a 2D surface in 3D space, such as the mesh used in 2D hydraulics
        models.

        To probe on them, the surface has to be flattened first.

        Parameters
        ----------
        lineVTK: {vtkLineSource} -- coordinates of points in the vtkLineSource; the points don't need to ordered,
        thus they can be just a bunch of points
        readerUnstructuredGridVTK: {vtkUnstructuredGridReader} -- Unstructured Grid VTK reader
        varName: name of the variable to be probed

        Returns
        -------
        points: numpy arrays [number of points, 3]; points on the profile
        probed result array:
        elev: elevation (z) of points in the profile

        """

        # Get data from the Unstructured Grid VTK reader
        data = readerUnstructuredGridVTK.GetOutput()
        bounds = data.GetBounds()

        print("Unstructured Grid VTK bounds = ", bounds)

        ### make a transform to set all Z values to zero ###
        flattener = vtk.vtkTransform()
        flattener.Scale(1.0, 1.0, 0.0)

        ### flatten the input in case it's not already flat ###
        i_flat = vtk.vtkTransformFilter()
        i_flat.SetInputConnection(lineVTK.GetOutputPort())
        i_flat.SetTransform(flattener)

        ### transfer z elevation values to the source's point scalar data ###
        s_elev = vtk.vtkElevationFilter()
        s_elev.SetInputData(data)
        s_elev.SetHighPoint(0, 0, bounds[5])
        s_elev.SetLowPoint(0, 0, bounds[4])
        s_elev.SetScalarRange(bounds[4], bounds[5])
        s_elev.Update()

        # print("s_elev = ", s_elev.GetUnstructuredGridOutput())

        ### flatten the source data; the Z elevations are already in the scalars data ###
        s_flat = vtk.vtkTransformFilter()
        s_flat.SetInputConnection(s_elev.GetOutputPort())
        s_flat.SetTransform(flattener)

        # vtkProbeFilter, the probe line is the input, and the underlying dataset is the source.
        probe = vtk.vtkProbeFilter()
        probe.SetInputConnection(i_flat.GetOutputPort())
        probe.SetSourceConnection(s_flat.GetOutputPort())
        probe.Update()

        # get the data from the VTK-object (probe) to an numpy array
        print("varName =", varName)

        varProbedValues = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray(varName))

        numPoints = probe.GetOutput().GetNumberOfPoints()  # get the number of points on the line

        # get the elevation from the VTK-object (probe) to an numpy array
        elev = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray("Elevation"))

        # intialise the points on the line
        x = np.zeros(numPoints)
        y = np.zeros(numPoints)
        z = np.zeros(numPoints)
        points = np.zeros((numPoints, 3))

        # get the coordinates of the points on the line
        for i in range(numPoints):
            x[i], y[i], z[i] = probe.GetOutput().GetPoint(i)
            points[i, 0] = x[i]
            points[i, 1] = y[i]
            points[i, 2] = z[i]
        return points, varProbedValues, elev

    def vtk_diff_consistent(self, vtkFileName1, vtkFileName2, vtkFileNameDiff,
                            varName1, varName2, varNameDiff):
        """ Calculate the difference of a variable in two VTK files with consistent Unstructured Grid.

        vrkNameDiff = vtkFileName1 - vtkFileName2

        Both vtkFileName1 and vtkFileName2 should be Unstructured Grid VTK files.

        Parameters
        ----------
        vtkFileName1: file name of first VTK
        vtkFileName2: file name of second VTK
        vtkFileNameDiff: file name of VTK to save the difference
        varName1: variable name in first VTK
        varName2: variable name in second VTK
        varNameDiff: variable name for the difference
        (not implemented) diffNodal: whether varNameDiff should be calculated at node or cell center.

        Returns
        -------

        """

        vtkFileReader1 = self.readVTK_UnstructuredGrid(vtkFileName1)
        vtkFileReader2 = self.readVTK_UnstructuredGrid(vtkFileName2)

        # Get data from the Unstructured Grid VTK reader
        data1 = vtkFileReader1.GetOutput()
        data2 = vtkFileReader2.GetOutput()

        #check consistency
        if (data1.GetNumberOfCells()!=data2.GetNumberOfCells()) or \
                (data1.GetNumberOfPoints()!=data2.GetNumberOfPoints()):
           print("Vtk files 1 and 2 are not consistent. They have different number of cells and/or points. Exiting ...")
           sys.exit()

        #check whether the request varName exists in the VTK's cell data. If not check its point data.
        #If varName exists in the point data, interpolate it to cell center. Otherwise, report error.
        pointDataNames1, cellDataNames1 = self.get_uGRid_all_field_names(data1)

        print("pointDataNames, cellDataNames in first VTK file =", pointDataNames1, cellDataNames1)

        if (not varName1 in cellDataNames1): #varName1 does not exist in VTK file1's cell data:
            if (varName1 in pointDataNames1): #varName1 exists in VTK file1's point data; interpolate that to cell center
                self.vtk_point_to_cell_interpolation(data1)
            else:
                print("varName1 does not exist in vtkFileName1's cell data or point data. Exiting...")
                sys.exit()

        pointDataNames2, cellDataNames2 = self.get_uGRid_all_field_names(data2)

        print("pointDataNames, cellDataNames in second VTK file =", pointDataNames2, cellDataNames2)

        if (not varName2 in cellDataNames2):  # varName2 does not exist in VTK file1's cell data:
            if (varName2 in pointDataNames2):  # varName1 exists in VTK file1's point data; interpolate that to cell center
                self.vtk_point_to_cell_interpolation(data2)
            else:
                print("varName2 does not exist in vtkFileName2's cell data or point data. Exiting...")
                sys.exit()


        varName1_data = VN.vtk_to_numpy(data1.GetCellData().GetArray(varName1))
        varName2_data = VN.vtk_to_numpy(data2.GetCellData().GetArray(varName2))

        varDiff_data = varName1_data - varName2_data

        varDiff_data_vtk = VN.numpy_to_vtk(varDiff_data)
        varDiff_data_vtk.SetName(varNameDiff)

        #print(varDiff_data)


        #create a copy of the unstructured grid vtk
        uGrid = vtk.vtkUnstructuredGrid()
        uGrid.DeepCopy(data1)

        #get the names of all point and cell data in uGrid
        pointDataNames_temp, cellDataNames_temp = self.get_uGRid_all_field_names(uGrid)

        #remove all uGrid's point and cell data
        vtkPointData = uGrid.GetPointData()
        for i in range(len(pointDataNames_temp)):
            vtkPointData.RemoveArray(pointDataNames_temp[i])

        vtkCellData = uGrid.GetCellData()
        for i in range(len(cellDataNames_temp)):
            vtkCellData.RemoveArray(cellDataNames_temp[i])

        #add the diff field to the cell data
        vtkCellData.AddArray(varDiff_data_vtk)

        # write out the diff result to vtk
        self.writeVTK_UnstructuredGrid(uGrid, vtkFileNameDiff)

    def get_uGRid_all_field_names(self, uGrid):
        """ Get the field names of point data and cell data

        Parameters
        ----------
        uGrid: vtkUnstructuredGrid object

        Returns
        -------
        pointDataNames: {list} -- list of point data names
        cellDataNames: {list} -- list of cell data names

        """
        pointDataNames = []
        cellDataNames = []

        pointdata = uGrid.GetPointData()

        for i in range(pointdata.GetNumberOfArrays()):
            pointDataNames.append(pointdata.GetArrayName(i))

        celldata = uGrid.GetCellData()

        for i in range(celldata.GetNumberOfArrays()):
            cellDataNames.append(celldata.GetArrayName(i))

        return pointDataNames, cellDataNames

    def vtk_cell_to_point_interpolation(self, uGrid, varName=''):
        """Interpolate all cell data to point in an unstructured grid

        In VTK, it is simply an avarege of all data of cells sharing a point

        Parameters
        ----------
        uGrid: vtkUnstructuredGrid object
        varName: {string} -- optional variable name; if specified, only this variable will be processed (not
        implemented yet)


        Returns
        -------

        """

        #first check whether there is already a point data named "varName"
        pointDataNames, cellDataNames = self.get_uGRid_all_field_names(uGrid)

        print("pointDataNames, cellDataNames =", pointDataNames, cellDataNames)

        if varName in pointDataNames:
            print("varName = ", varName, "is already in the point data. Nothing to be done.")
            return

        c2p_interpolator = vtk.vtkCellDataToPointData()
        c2p_interpolator.SetInputData(uGrid)
        c2p_interpolator.Update()

        pointdata = c2p_interpolator.GetOutput().GetPointData()

        #print("Interpolated point data = ", pointdata)

        uGridPointData = uGrid.GetPointData()

        for i in range(pointdata.GetNumberOfArrays()):
            uGridPointData.AddArray(pointdata.GetArray(i))



    def vtk_point_to_cell_interpolation(self, uGrid, varName=''):
        """Interpolate all point data to cell in an unstructured grid

        In VTK, it is a simple average of data for all points defining a cell.

        Parameters
        ----------
        uGrid: vtkUnstructuredGrid object
        varName: {string} -- optional variable name; if specified, only this variable will be processed (not
        implemented yet)


        Returns
        -------

        """

        #first check whether there is already a cell data named "varName"
        pointDataNames, cellDataNames = self.get_uGRid_all_field_names(uGrid)

        print("pointDataNames, cellDataNames =", pointDataNames, cellDataNames)

        if varName in cellDataNames:
            print("varName = ", varName, "is already in the point data. Nothing to be done.")
            return

        p2c_interpolator = vtk.vtkPointDataToCellData()
        p2c_interpolator.SetInputData(uGrid)
        p2c_interpolator.Update()

        celldata = p2c_interpolator.GetOutput().GetCellData()

        #print("Interpolated cell data = ", celldata)

        uGridCellData = uGrid.GetCellData()

        for i in range(celldata.GetNumberOfArrays()):
            uGridCellData.AddArray(celldata.GetArray(i))


