import numpy as np
import csv
import vtk

from vtk.util import numpy_support as VN

import pyHMT2D

from ..__common__ import pyHMT2D_SCALAR, pyHMT2D_VECTOR

class Measurement(object):
    """ Measurement data base class

    Attributes
    ----------
        name : str
            name of the measurement data
        type : str, optional
            type of the measurement, e.g., "point", "line"

    """

    def __init__(self, name, type=""):
        self.name = name
        self.type = type

    def getName(self):
        return self.name

    def getType(self):
        return self.type

    def __str__(self):
        return "Measurement data name: %s " % (self.name)

    def init_data(self):
        """ Initialize the data

        Returns
        -------

        """
        pass


class PointMeasurement(Measurement):
    """Point measurement data

    Point measurement contains measured data at a list of points. Each point measurement
    can only have one variable, such as velocity (angle + magnitude). For different
    variables, different point measurement should be created.

    The point measurement data should have the following csv format for velocity (vector)
    Name, x, y, angle, velocity
    point1, xxx, xxx, xxx, xxx
    point2, xxx, xxx, xxx, xxx
    ...

    or the following for scalar

    Name, x, y, wse
    point1, xxx, xxx, xxx
    point2, xxx, xxx, xxx
    ...


    Attributes
    ----------

    """

    def __init__(self, name, pointMeasurement_filename):
        """PointMeasurement class constructor

        Parameters
        ----------
        name  : str
            name of the point measurement, e.g., "stage", "velocity"
        pointMeasurement_filename : str
            name of the file that contains point measurement.

        """

        Measurement.__init__(self, name, "point")

        self.pointMeasurement_filename = pointMeasurement_filename

        #measurement data type (default scalar)
        self.data_type = pyHMT2D_SCALAR

        #header of the point measurement data
        self.field = []

        #rows for each point (name, x, y, data)
        self.rows = []

        #load the measurement data
        self.load_measurement_data()

        #debug
        self.outputToVTK()

    def load_measurement_data(self):
        """Load the measurement data in the specified file (in csv format)

        Returns
        -------

        """

        print("Load point measurement data file", self.pointMeasurement_filename)

        with open(self.pointMeasurement_filename, 'r') as csvfile:
            # creating a csv reader object
            csvreader = csv.reader(csvfile)

            # extracting field names in first row
            self.fields = next(csvreader)

            # measurement data type: scalar or vector
            # it depends on the number of columns (4-scalar, 5-vector)
            if len(self.fields) == 4:
                self.data_type = pyHMT2D_SCALAR
            elif len(self.fields) == 5:
                self.data_type = pyHMT2D_VECTOR
            else:
                raise Exception("The number of columns in the point measurement data file needs to be either 4 or 5. Exit ...")

            # loop and extracting each data row
            for row in csvreader:
                self.rows.append(row)

            # get total number of measurement points
            print("Total number of measurement points: %d" % (csvreader.line_num-1))

        # printing the field names in the file
        print('Header names are:' + ', '.join(field for field in self.fields))

        #  printing all rows
        print('\nPoint measurement data:\n')
        for row in self.rows:
            print(', '.join(col for col in row))

    def get_measurement_points(self):
        """Get the measurement points in the form of Numpy 2D array

        Returns
        -------

        """

        measurement_points = np.zeros((len(self.rows), 2))

        for pointI in range(len(self.rows)):
            measurement_points[pointI, 0] = self.rows[pointI][1]
            measurement_points[pointI, 1] = self.rows[pointI][2]

        return measurement_points

    def get_measurement_points_as_vtkPoints(self):
        """Get the measurement points in the form of vtkPoints

        Returns
        -------

        """
        points = vtk.vtkPoints()

        for pointI in range(len(self.rows)):
            points.InsertNextPoint(float(self.rows[pointI][1]),
                                   float(self.rows[pointI][2]),
                                   0.0)

        return points

    def get_measurement_data(self):
        """ Get the measurement data in numpy array format

        Returns
        -------

        """

        if self.data_type == pyHMT2D_SCALAR:
            measurement_data = np.zeros(len(self.rows))

            for pointI in range(len(self.rows)):
                measurement_data[pointI] = self.rows[pointI][3]

        elif self.data_type == pyHMT2D_VECTOR:
            measurement_data = np.zeros((len(self.rows), 2))

            for pointI in range(len(self.rows)):
                angle = self.rows[pointI][3]
                mag = self.rows[pointI][4]

                measurement_data[pointI, 0] = mag * np.sin(np.deg2rad(angle))
                measurement_data[pointI, 1] = mag * np.cos(np.deg2rad(angle))

        return  measurement_data

    def outputToVTK(self, dir=''):
        """Output point measurement data to vtkUnstructuredGrid


        Parameters
        ----------
            dir : str, optional
                directory to write to

        Returns
        -------
        vtkFileName : str
            name of the output vTK file

        """
        print("Output point measurement data", self.name, "to VTK ...")

        vtkFileName = ''

        if dir!='':
            vtkFileName = dir + '/' + 'PointMeasurement_' + self.name + '.vtk' #use the case name
        else:
            vtkFileName = 'PointMeasurement_' + self.name + '.vtk'
        print("vtkFileName = ", vtkFileName)

        #build result variable name
        resultVarName = self.fields[3] if self.data_type == pyHMT2D_SCALAR else self.fields[4]

        print("Measurement variable name: ", resultVarName)

        # numpy array for measurement data
        resultData_scalar = np.zeros(len(self.rows), dtype="float32")

        resultData_vector = np.zeros((len(self.rows),3), dtype="float32")

        if self.data_type == pyHMT2D_SCALAR:
            for pointI in range(len(self.rows)):
                resultData_scalar[pointI] = float(self.rows[pointI][3])

        elif self.data_type == pyHMT2D_VECTOR:
            for pointI in range(len(self.rows)):
                angle = float(self.rows[pointI][3])
                mag = float(self.rows[pointI][4])

                resultData_vector[pointI, 0] = mag * np.sin(np.deg2rad(angle))
                resultData_vector[pointI, 1] = mag * np.cos(np.deg2rad(angle))
                resultData_vector[pointI, 2] = 0.0

        # build VTK object:
        # points
        pointsVTK = self.get_measurement_points_as_vtkPoints()

        polydataVTK = vtk.vtkPolyData()
        polydataVTK.SetPoints(pointsVTK)

        point_data = polydataVTK.GetPointData()  # This holds point data

        if self.data_type == pyHMT2D_SCALAR:
            temp_point_data_array = VN.numpy_to_vtk(resultData_scalar)
            temp_point_data_array.SetName(resultVarName)
            point_data.AddArray(temp_point_data_array)
        elif self.data_type == pyHMT2D_VECTOR:
            temp_point_data_array = VN.numpy_to_vtk(resultData_vector)
            temp_point_data_array.SetName(resultVarName)
            point_data.AddArray(temp_point_data_array)


        # write to vtk file
        unstr_writer = vtk.vtkPolyDataWriter()
        unstr_writer.SetFileName(vtkFileName)
        unstr_writer.SetInputData(polydataVTK)
        unstr_writer.Write()

        return vtkFileName
