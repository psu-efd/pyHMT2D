import sys
import numpy as np
import vtk
from vtk.util import numpy_support as VN
from scipy import interpolate
import os.path

import json

from pyHMT2D.Hydraulic_Models_Data import HydraulicData
from pyHMT2D.Misc.vtk_utilities import vtkHandler

from pyHMT2D.__common__ import *


class Backwater_1D_Data(HydraulicData):
    """A class for Backwater-1D data I/O, manipulation, and format conversion

    Attributes
    ----------


    Methods
    -------


    """

    def __init__(self, config_json_file_name):
        """Backwater_1D_Data class constructor

        Parameters
        ----------
        config_json_file_name : str
            name of the JSON configuration file
        """

        HydraulicData.__init__(self, "Backwater-1D")

        self.config_json_file_name = config_json_file_name

        #dictionary to hold information from the configuration JSON file
        self.configuration = {}

        #case set up parameters (read in from configuration JSON file)
        self.case_name = ""     #case name
        self.units = ""         #units: SI or EN
        self.startx = 0.0       #start point x
        self.startH = 0.0       #start water depth H
        self.startZ = 0.0       #start bed elevation Z
        self.riverLength = 0.0  #total river length
        self.slope = 0.0        #river slope
        self.ManningNZones = {} #dictionary to hold Manning's n zone information
        self.nGrid = 1          #number of grid
        self.specificDischarge = 0.0  #specific discharge

        #load configuration data
        self.load_configuration_data()

        #Manning's n interpolator
        self.ManningNFunc = None

        #build Manning's n interpolator
        self.build_ManningNFunc()

        #result block
        self.normalDepth = 0.0
        self.criticalDepth = 0.0
        self.gridx = None
        self.waterDepth = None
        self.WSE = None
        self.gridManningN = None

    def load_configuration_data(self):
        """Load only the configuration data (not result) from the JSON file

        Returns
        -------

        """

        print("Load Backwater-1D model configuration from file", self.config_json_file_name)

        with open(self.config_json_file_name) as f:
            self.configuration = json.load(f)

        #print(self.configuration)

        #print some information about the calibration configuration
        if gVerbose: print("Configuration for the Backwater-1D model:")
        if gVerbose: print(json.dumps(self.configuration["Backwater-1D"], indent=4, sort_keys=False))

        self.case_name = self.configuration["Backwater-1D"]["case_name"]
        self.units = self.configuration["Backwater-1D"]["units"]
        self.startx = self.configuration["Backwater-1D"]["startx"]
        self.startH = self.configuration["Backwater-1D"]["startH"]
        self.startZ = self.configuration["Backwater-1D"]["startZ"]
        self.riverLength = self.configuration["Backwater-1D"]["riverLength"]
        self.slope = self.configuration["Backwater-1D"]["slope"]
        self.ManningNZones = self.configuration["Backwater-1D"]["ManningNZones"]
        self.nGrid = self.configuration["Backwater-1D"]["nGrid"]
        self.specificDischarge = self.configuration["Backwater-1D"]["specificDischarge"]

        # build the Manning's n interpolation function
        self.build_ManningNFunc()


    def build_ManningNFunc(self):
        """ Define the Manning's n interpolation function based on the Manning's n zone information

        Returns
        -------

        """

        # create the lists of x and Manning's n value
        x = []
        ManningN = []

        # loop over Manning's N zones
        for ManningNZone in self.ManningNZones:
            #check if the x coordinate is already in the list. Make sure
            #the zones are defined with some small gap in between.
            if ManningNZone["startx"] in x or ManningNZone["endx"] in x:
                raise Exception("startx and endx of Manning's zones should not overlap. Leave a small gap"
                                "between different zones so the interpolation works properly. Exiting ...")

            x.append(ManningNZone["startx"])
            ManningN.append(ManningNZone["n"])

            x.append(ManningNZone["endx"])
            ManningN.append(ManningNZone["n"])

        self.ManningNFunc = interpolate.interp1d(x, ManningN, bounds_error=False, fill_value="extrapolate")


    def write_configuration_data_to_JSON(self, json_file_name):
        """Write the configuration data only (not result) to JSON file

        Parameters
        ----------
            json_file_name : str
                name of the JSON file to write to

        Returns
        -------

        """

        with open(json_file_name, "w") as json_file:
            json.dump(self.configuration, json_file, indent=4, sort_keys=False)

    def modify_ManningsN(self, materialID, newManningsNValue):
        """Modify materialID's Manning's n value to new value

        Parameters
        ----------
        materialID : int
            material ID
        newManningsNValue : float
            new Manning's n value

        Returns
        -------

        """

        if gVerbose: print("Modify Manning's n value ...")

        if not isinstance(materialID, int):
            print("Material ID has to be an integer. The type of materialID passed in is ", type(materialID),
                  ". Exit.\n")

        if not isinstance(newManningsNValue, float):
            print("Manning's n has to be a float. The type of newManningsNValue passed in is ", type(newManningsNValue),
                  ". Exit.\n")

        bFound = False

        #loop through all ManningNZones to find the materialID
        for zoneI in range(len(self.ManningNZones)):
            if materialID == self.ManningNZones[zoneI]["materialID"]:
                bFound = True

                if gVerbose: print("    Old Manning's n value =", self.ManningNZones[zoneI]["n"], "for material ID = ",
                                  materialID, "zone name = ", self.ManningNZones[zoneI]["name"])

                self.ManningNZones[zoneI]["n"] = newManningsNValue

                if gVerbose: print("    New Manning's n value =", self.ManningNZones[zoneI]["n"], "for material ID = ",
                                  materialID, "zone name = ", self.ManningNZones[zoneI]["name"])

        #if didn't find the specified materialID, something is wrong
        if not bFound:
            raise Exception("The specified materialID", materialID, "is not in the Manning's n list. Please check.")

        #update the Manning's n interpolator
        self.build_ManningNFunc()


    def outputResultToVTK(self, dir=''):
        """Output result data to vtkUnstructuredGrid

        This has to be called after the result has been obtained by running the model

        Parameters
        ----------
            dir : str, optional
                directory to write to

        Returns
        -------
        vtkFileName : str
            name of the output vTK file

        """

        if gVerbose: print("Output result data to VTK ...")

        if (len(self.waterDepth) == 0):
            print("Empty solution arrays. Run the model first. Exiting ...")
            sys.exit()

        vtkFileName = ''

        if dir!='':
            vtkFileName = dir + '/' + 'Backwater1D_' + self.case_name + '.vtk' #use the case name
        else:
            vtkFileName = 'Backwater1D_' + self.case_name + '.vtk'
        if gVerbose: print("vtkFileName = ", vtkFileName)

        #build result variable names
        resultVarNames = []
        if self.units == "SI":
            resultVarNames = ["Water_Elev_m","Velocity_m_p_s","ManningN","Bed_Elev_m"]
        elif self.units == "EN":
            resultVarNames = ["Water_Elev_ft", "Velocity_ft_p_s","ManningN","Bed_Elev_ft"]
        else:
            raise Exception("Units system is not set.")

        #print("All nodal solution variable names: ", resultVarNames)

        # numpy array for all solution variables (except velocity)
        resultData = np.zeros((self.nGrid, len(resultVarNames)-1), dtype="float32")

        # numpy array for velocity
        resultVelocity = np.zeros((self.nGrid, 3), dtype="float32")

        #add WSE
        resultData[:,0] = self.WSE

        #add Velocity
        resultVelocity[:,0] = self.specificDischarge/self.waterDepth
        resultVelocity[:,1] = 0.0
        resultVelocity[:,2] = 0.0

        #add ManningN
        resultData[:,1] = self.gridManningN

        #add bed elevation
        resultData[:,2] = self.WSE - self.waterDepth

        # build VTK object:
        # points
        pointsVTK = vtk.vtkPoints()

        gridCordinates = np.zeros((self.nGrid, 3))
        for gridI in range(self.nGrid):
            gridCordinates[gridI,0] = self.gridx[gridI]
            gridCordinates[gridI,1] = 0.0
            gridCordinates[gridI,2] = (self.WSE - self.waterDepth)[gridI]

        pointsVTK.SetData(VN.numpy_to_vtk(gridCordinates))

        # cell topology information list: [num. of nodes, node0, node1, .., num. of nodes, nodexxx]
        # the list start with the number of nodes for a cell and then the list of node indexes
        connectivity_list = []

        # type of cells (contains the number of points)
        cellFPCounts = np.zeros(self.nGrid-1, dtype=np.int64)

        # loop over all elements (lines)
        for k in range(self.nGrid-1):
            connectivity_list.append(2)  #VTK_LINE

            connectivity_list.append(k)
            connectivity_list.append(k+1)

            cellFPCounts[k] = 2

        connectivity = np.array(connectivity_list, dtype=np.int64)

        # convert cell's number of face points to VTK cell type
        vtkHandler_obj = vtkHandler()
        cell_types = vtkHandler_obj.number_of_nodes_to_vtk_celltypes(cellFPCounts)

        cellsVTK = vtk.vtkCellArray()
        cellsVTK.SetCells(self.nGrid-1, VN.numpy_to_vtkIdTypeArray(connectivity))

        uGrid = vtk.vtkUnstructuredGrid()
        uGrid.SetPoints(pointsVTK)
        uGrid.SetCells(cell_types, cellsVTK)

        point_data = uGrid.GetPointData()  # This holds point data

        # add WSE
        temp_point_data_array = VN.numpy_to_vtk(resultData[:,0])
        temp_point_data_array.SetName(resultVarNames[0])
        point_data.AddArray(temp_point_data_array)

        # add velocity
        temp_point_data_array = VN.numpy_to_vtk(resultVelocity)
        temp_point_data_array.SetName(resultVarNames[1])
        point_data.AddArray(temp_point_data_array)

        # add ManningN
        temp_point_data_array = VN.numpy_to_vtk(resultData[:,1])
        temp_point_data_array.SetName(resultVarNames[2])
        point_data.AddArray(temp_point_data_array)

        # add bed elevation
        temp_point_data_array = VN.numpy_to_vtk(resultData[:,2])
        temp_point_data_array.SetName(resultVarNames[3])
        point_data.AddArray(temp_point_data_array)


        # write to vtk file
        unstr_writer = vtk.vtkUnstructuredGridWriter()
        unstr_writer.SetFileName(vtkFileName)
        unstr_writer.SetInputData(uGrid)
        unstr_writer.Write()

        return vtkFileName

