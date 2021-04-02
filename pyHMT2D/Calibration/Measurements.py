import numpy as np
import csv
import copy

import vtk

from vtk.util import numpy_support as VN


from ..__common__ import pyHMT2D_SCALAR, pyHMT2D_VECTOR, gVerbose

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

        #rows for measurement data at each point (name, x, y, data)
        self.rows = []

        # measurement points in numpy array
        self.measurementPoints = None

        # measurement data in numpy array
        self.measurementData = None

        #simulation data in numpy array
        #for comparision and output
        self.simulationData = None

        #load the measurement data
        self.load_measurement_data()

        #debug
        self.outputToVTK()

    def load_measurement_data(self):
        """Load the measurement data in the specified file (in csv format)

        Returns
        -------

        """

        if gVerbose: print("Load point measurement data file", self.pointMeasurement_filename)

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
            if gVerbose: print("Total number of measurement points: %d" % (csvreader.line_num-1))

        # printing the field names in the file
        if gVerbose: print('Header names are:' + ', '.join(field for field in self.fields))

        # printing all rows
        if gVerbose: print('\nPoint measurement data:\n')
        for row in self.rows:
            if gVerbose: print(', '.join(col for col in row))


        # build the measurement points numpy array
        self.measurementPoints = np.zeros((len(self.rows), 2))

        for pointI in range(len(self.rows)):
            self.measurementPoints[pointI, 0] = float(self.rows[pointI][1])
            self.measurementPoints[pointI, 1] = float(self.rows[pointI][2])

        # build the measurement data numpy array
        if self.data_type == pyHMT2D_SCALAR:
            self.measurementData = np.zeros(len(self.rows))

            for pointI in range(len(self.rows)):
                self.measurementData[pointI] = float(self.rows[pointI][3])

        elif self.data_type == pyHMT2D_VECTOR:
            self.measurementData = np.zeros((len(self.rows), 3))  # 3 here; be z-velocity =0 and does not matter.

            for pointI in range(len(self.rows)):
                angle = float(self.rows[pointI][3])
                mag = float(self.rows[pointI][4])

                self.measurementData[pointI, 0] = mag * np.sin(np.deg2rad(angle))
                self.measurementData[pointI, 1] = mag * np.cos(np.deg2rad(angle))

        # make simulationData a copy of measurementData (which is used to store
        # corresponding simulation data in calibration)
        self.simulationData = copy.deepcopy(self.measurementData)

    def get_measurement_points(self):
        """Get the measurement points in the form of Numpy 2D array

        Returns
        -------

        """

        return self.measurementPoints

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

        return  self.measurementData

    def set_simulation_results(self, simulationData):
        """ Set the simulation results at sampling points to rows_simulation

        Parameters
        ----------
        simulationData : numpy array
            numpy array for the simulation results at measurement points

        Returns
        -------

        """

        #check the dimension of simulationData and its consistency with measurement data
        if simulationData.shape != self.simulationData.shape:
            raise  Exception("In set_simulation_results, the simulation data have different shape. Exiting ...")

        self.simulationData = copy.deepcopy(simulationData)

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

        if gVerbose: print("Output point measurement data", self.name, "to VTK ...")

        vtkFileName = ''

        if dir!='':
            vtkFileName = dir + '/' + 'PointMeasurement_' + self.name + '.vtk' #use the case name
        else:
            vtkFileName = 'PointMeasurement_' + self.name + '.vtk'
        if gVerbose: print("vtkFileName = ", vtkFileName)

        #build result variable name
        resultVarName = self.fields[3] if self.data_type == pyHMT2D_SCALAR else self.fields[4]

        if gVerbose: print("Measurement variable name: ", resultVarName)

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


    def outputSimulationResultsToCSV(self, dir=''):
        """Output simulation results at measurement points to csv (for postprocessing)


        Parameters
        ----------
            dir : str, optional
                directory to write to

        Returns
        -------
        csvFileName : str
            name of the output csv file

        """

        if gVerbose: print("Output simulated results at measurement points", self.name, "to csv ...")

        csvFileName = ''

        if dir!='':
            csvFileName = dir + '/' + 'PointSimulation_' + self.name + '.csv' #use the case name
        else:
            csvFileName = 'PointSimulation_' + self.name + '.csv'
        if gVerbose: print("csvFileName = ", csvFileName)

        #build result variable name
        resultVarName = self.fields[3] if self.data_type == pyHMT2D_SCALAR else self.fields[4]

        if gVerbose: print("Simulation variable name: ", resultVarName)

        if self.data_type == pyHMT2D_VECTOR:
            resultData_vector = np.zeros((len(self.rows),2), dtype="float32")

            #azimuth angle
            resultData_vector[:,0] = np.degrees(np.arctan(self.simulationData[:,0]/(self.simulationData[:,1]+1e-8)))
            #velocity magnitude
            resultData_vector[:,1] = np.sqrt(np.power(self.simulationData[:,0],2) + np.power(self.simulationData[:,1],2))

        fid = open(csvFileName, "w")

        if self.data_type == pyHMT2D_SCALAR: #scalar
            fid.write("Name,x,y,%s\n"%resultVarName)

            #loop through each point
            for pointI in range(len(self.rows)):
                fid.write("%s, %f, %f, %f\n" % (self.rows[pointI][0],
                                                float(self.rows[pointI][1]),
                                                float(self.rows[pointI][2]),
                                                self.simulationData[pointI]
                                                ))
        elif self.data_type == pyHMT2D_VECTOR:
            fid.write("Name,x,y,%s,%s\n" % (self.fields[3], self.fields[4]))

            #loop through each point
            for pointI in range(len(self.rows)):
                fid.write("%s, %f, %f, %f, %f\n" % (self.rows[pointI][0],
                                                float(self.rows[pointI][1]),
                                                float(self.rows[pointI][2]),
                                                resultData_vector[pointI,0],
                                                resultData_vector[pointI,1]
                                                ))

        fid.close()

        return csvFileName
