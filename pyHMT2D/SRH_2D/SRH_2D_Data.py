# -*- coding: utf-8 -*-
""" A set of classes to handle SRH-2D data

SRH-2D data include information in SRHHydro, SRHGeom, SRHMat, and simulation results.

Author: Xiaofeng Liu, PhD, PE
Penn State University
"""

import sys
import os
import copy
import numpy as np
import h5py
from scipy import interpolate
from osgeo import gdal
from os import path
import shlex
import vtk
from vtk.util import numpy_support as VN

from .helpers import *
from ..__common__ import *
from ..Misc.vtk_utilities import vtkCellTypeMap
from ..Misc import printProgressBar, vtkHandler

class SRH_2D_SRHHydro:
    """A class to handle srhhydro file for SRH-2D

    Attributes:


    Methods:

    """

    def __init__(self, srhhydro_filename):
        """ Constructor for SRH_2D_SRHHydro

        The SRH_2D_SRHHydro file contains all information for a SRH-2D run.

        Parameters
        ----------
        srhhydro_filename: str
            The name of the SRHHydro file
        """

        self.srhhydro_filename = srhhydro_filename #srhhydro file name

        #dict to hold the content of the srhhydro file
        self.srhhydro_content = {}

        #parse the srhhydro file and build srhhydro_content
        self.parse_srhhydro_file()

    def parse_srhhydro_file(self):
        """ Parse the SRHHydro file

        It will read the SRHHydro file and build the dictionary self.srhhydro_content

        Returns
        -------

        """

        res_all = {}
        res_ManningsN = {}  # dict for ManningsN (there cuould be multiple entries)
        res_BC = {}  # dict for BC (there could be multiple entries)
        res_MONITORING = {} #dict for monitoring lines (it is under BC in SRHHYDRO)
        res_IQParams = {}  # dict for subcritical inlet discharge boundary condition
        res_ISupCrParams = {}  # dict for supercritical inlet, the same as IQParams, with the additon of WSE
        res_EWSParamsC = {}  # dict for stage exit boundary condition (constant)
        res_EWSParamsRC = {}  # dict for stage exit boundary condition (rating curve)
        res_EQParams = {}  # dict for exit discharge boundary condition
        res_NDParams = {}  # dict for normal depth outlet boundary condition

        #check whether the srhhydro file exists
        if not path.isfile(self.srhhydro_filename):
            print("The SRHHYDRO file", self.srhhydro_filename, "does not exists. Exitting ...")
            sys.exit()

        for line in open(self.srhhydro_filename):
            #parts = line.strip().split(' ')
            parts = shlex.split(line.strip())

            #print(parts)

            map(str.strip, parts)

            #print(parts)

            if len(parts) <= 1:  # if there is only one word, assume it is a comment; do nothing
                continue

            if parts[0] == 'ManningsN':
                res_ManningsN[int(parts[1])] = float(parts[2])
            elif parts[0] == 'BC': #contains both boundary conditions and monitoring lines ("MONITORING")
                if parts[2] == 'MONITORING': #this is a monitoring line, not a BC
                    res_MONITORING[int(parts[1])] = parts[2]
                else:                        #this is a real BC
                    res_BC[int(parts[1])] = parts[2]
            elif parts[0] == 'IQParams':
                res_IQParams[int(parts[1])] = [parts[2], parts[3], parts[4]]
            elif parts[0] == 'ISupCrParams': #need to check these
                res_ISupCrParams[int(parts[1])] = parts[2]
            elif parts[0] == 'EWSParamsC':
                res_EWSParamsC[int(parts[1])] = [parts[2], parts[3], parts[4]]
            elif parts[0] == 'EWSParamsRC':
                res_EWSParamsRC[int(parts[1])] = [parts[2], parts[3], parts[4]]
            elif parts[0] == 'EQParams': #need to check these
                res_EQParams[int(parts[1])] = parts[2]
            elif parts[0] == 'NDParams': #need to check these
                res_NDParams[int(parts[1])] = parts[2]
            elif parts[0] == 'OutputFormat':
                res_all[parts[0]] = [parts[1],parts[2]]
            elif parts[0] == 'SimTime':
                res_all[parts[0]] = [float(parts[1]),float(parts[2]),float(parts[3])]
            elif parts[0] == 'ParabolicTurbulence':
                res_all[parts[0]] = float(parts[1])
            elif parts[0] == 'OutputOption':
                res_all[parts[0]] = int(parts[1])
            elif parts[0] == 'OutputInterval':
                res_all[parts[0]] = float(parts[1])
            else:
                #res_all[parts[0].lstrip()] = parts[1].strip('"').split(' ')
                res_all[parts[0].lstrip()] = parts[1]

        # add ManningsN, BC and all other sub-dict to res_all
        if res_ManningsN:
            res_all['ManningsN'] = res_ManningsN

        if res_BC:
            res_all['BC'] = res_BC

        if res_MONITORING:
            res_all['MONITORING'] = res_MONITORING

        if res_IQParams:
            res_all['IQParams'] = res_IQParams

        if res_ISupCrParams:
            res_all['ISupCrParams'] = res_ISupCrParams

        if res_EWSParamsC:
            res_all['EWSParamsC'] = res_EWSParamsC

        if res_EWSParamsRC:
            res_all['EWSParamsRC'] = res_EWSParamsRC

        if res_EQParams:
            res_all['EQParams'] = res_EQParams

        if res_NDParams:
            res_all['NDParams'] = res_NDParams


        if False:
            print(res_all)

            if res_ManningsN: print(res_all['ManningsN'])
            if res_BC: print(res_all['BC'])
            if res_MONITORING: print(res_all['MONITORING'])
            if res_IQParams: print(res_all['IQParams'])
            if res_ISupCrParams: print(res_all['ISupCrParams'])
            if res_EWSParamsC: print(res_all['EWSParamsC'])
            if res_EWSParamsRC: print(res_all['EWSParamsRC'])
            if res_EQParams: print(res_all['EQParams'])
            if res_NDParams: print(res_all['NDParams'])

        self.srhhydro_content = res_all

    def modify_ManningsN(self, materialID, newManningsNValue):
        """Modify materialID's Manning's n value to newManningsNValue

        Parameters
        ----------
        materialID: {int} -- material ID
        newManningsNValue: {float} -- new Manning's n value

        Returns
        -------

        """

        print("Modify Manning's n value ...")

        if not isinstance(materialID, int):
            print("Material ID has to be an integer. The type of materialID passed in is ", type(materialID),
                  ". Exit.\n")

        if not isinstance(newManningsNValue, float):
            print("Manning's n has to be a float. The type of newManningsNValue passed in is ", type(newManningsNValue),
                  ". Exit.\n")


        nDict = self.srhhydro_content["ManningsN"]

        if materialID in nDict:
            print("    Old Manning's n value =", nDict[materialID], "for material ID = ", materialID)
            nDict[materialID] = newManningsNValue
            print("    New Manning's n value =", nDict[materialID], "for material ID = ", materialID)
        else:
            print("The specified materialID", materialID, "is not in the Manning's n list. Please check.")



    def write_to_file(self, new_srhhydro_file_name):
        """Write to a SRHHydro file (useful for modification of the SRHHydro file)

        Returns
        -------

        """

        print("Wring the SRHHYDRO file %s \n" % new_srhhydro_file_name)

        try:
            fid = open(new_srhhydro_file_name, 'w')
        except IOError:
            print('srhhydro file open error')
            sys.exit()

        for key, value in self.srhhydro_content.items():
            #print("key, value = ", key, value)
            if type(value) is dict: #if the current item is dictionary itself
                for subkey, subvalue in value.items():
                    if "ManningsN" in key:
                        fid.write("ManningsN " + str(subkey) + ' ' + str(subvalue) + '\n')
                    elif "MONITORING" in key:
                        fid.write("BC " + str(subkey) + ' ' + str(subvalue) + '\n')
                    elif "BC" in key:
                        fid.write("BC " + str(subkey) + ' ' + str(subvalue) + '\n')
                    elif "IQParams" in key:
                        fid.write("IQParams " + str(subkey) + ' ' + str(subvalue[0]) + ' ' + str(subvalue[1])+ ' ' +
                                  str(subvalue[2]) + '\n' )
                    elif "EWSParamsC" in key:
                        fid.write("EWSParamsC " + str(subkey) + ' ' + str(subvalue[0])+ ' ' + str(subvalue[1])+ ' ' +
                                  str(subvalue[2]) + '\n' )
                    elif "EWSParamsRC" in key:
                        fid.write("EWSParamsRC " + str(subkey) + ' \"' + str(subvalue[0])+ '\" ' + str(subvalue[1])+
                                  ' ' + str(subvalue[2]) + '\n' )
                    elif "ISupCrParams" in key:  #need to check these (no reference)
                        fid.write("ISupCrParams " + str(subkey) + ' ' + str(subvalue) + '\n' )
                    elif "EQParams" in key:      #need to check these (no reference)
                        fid.write("EQParams " + str(subkey) + ' ' + str(subvalue) + '\n' )
                    elif "NDParams" in key:      #need to check these (no reference)
                        fid.write("NDParams " + str(subkey) + ' ' + str(subvalue) + '\n' )
            else:
                if "Case" in key or "Description" in key or "Grid" in key \
                        or "HydroMat" in key or "MonitorPtFile" in key:
                    fid.write(str(key) + ' \"' + str(value) + '\"\n')
                elif "OutputFormat" in key:
                    fid.write(str(key) + ' ' + str(value[0]) + ' ' + str(value[1]) + '\n')
                elif "SimTime" in key:
                    fid.write(str(key) + ' ' + str(value[0]) + ' ' + str(value[1]) + ' ' + str(value[2]) + '\n')
                else:
                    #print("last, key, value", key, value, str(key), str(value))
                    fid.write(str(key) + ' ' + str(value) + '\n')

        fid.close()

    def get_simulation_start_end_time(self):
        """Get the start and end time of simulation (in hours, from the "SimTime" entry)

        Returns
        -------
        startTime (hr)
        endTime (hr)

        """

        value = self.srhhydro_content["SimTime"]

        return value[0],value[2]

    def get_simulation_time_step_size(self):
        """Get the time step size of simulation (in seconds, from the "SimTime" entry)

        Returns
        -------
        deltaT (s)

        """

        value = self.srhhydro_content["SimTime"]

        return value[1]

    def get_grid_file_name(self):
        """Get grid file name (shoudl be like "xxx.srhgeom")

        Returns
        -------
        str
            Name of the grid file

        """

        value = self.srhhydro_content["Grid"]

        return value

    def get_mat_file_name(self):
        """Get material file name (shoudl be like "xxx.srhmat")

        Returns
        -------
        str
            Name of the material file

        """

        value = self.srhhydro_content["HydroMat"]

        return value

class SRH_2D_SRHGeom:
    """A class to handle srhgeom file for SRH-2D

    Notes:
        1. In SMS (SRH-2D), monitoring lines (ML) and monitoring points (MP) are treated seperately and differently.
           ML is created in the same way as boundary lines. Nodes closest to the line are recorded as a NodeString
           in SRHGEOM file. MP is created as a feature point and its coordiantes are recorded in the "srhmpoint" file.
           Interpolation is done to probe the results at MPs.

           This potentially have some drawbacks:

           a. The MLs and MPs have to be defined before the simulation is run, not afterwards. One can use plot over a
              line or point, but it is not convenient and it is limited. For example, to calculate flow over a line,
              one has to cut a line and then integrate.
           b. The way MLs are defined if slightly limited by the mesh (node locations). To be exactly at a ML, one has
              to interpolate the simulated solution to the points on MLs.

        2. Becasue BC and ML are mixed together, we need to separate them. However, this can only be done with the BC
           information in the srhhydro file.

    Attributes
    ----------

    Methods
    ----------

    """

    def __init__(self, srhgeom_filename, bcDict):
        """Constructor for SRH_2D_SRHGeom

        Parameters
        ----------
        srhgeom_filename : str
            Name of the SRHGeom file
        bcDict : dict
            Boundary condtion dictionary
        """

        self.srhgeom_filename = srhgeom_filename

        #boundary condition dictionary. In srhgeom file, not all nodeStrings are boundary conditions. Some of them
        #could be monitoring lines. This bcDict is passed in to distinguish the two (srhgeom does not have this
        #information; only srhhydro has). An SRH_2D_SRHHydro object has to be previously created, and then
        #bcDict <- srhhydro_obj.srhhydro_content["BC"].
        self.bcDict = bcDict

        #mesh information:
        #   number of elements
        self.numOfElements = -1

        #   number of nodes
        self.numOfNodes = -1

        # number of NodeString
        self.numOfNodeStrings = -1

        # First read of the SRHGEOM file to get number of elements and nodes
        self.getNumOfElementsNodes()

        # list of nodes for all elements
        self.elementNodesList = np.zeros([self.numOfElements, gMax_Nodes_per_Element], dtype=int)

        # number of nodes for each element (3,4,...,gMax_Nodes_per_Element)
        self.elementNodesCount = np.zeros(self.numOfElements, dtype=int)

        # each element's vtk cell type
        self.vtkCellTypeCode = np.zeros(self.numOfElements, dtype=int)

        # each node's 3D coordinates
        self.nodeCoordinates = np.zeros([self.numOfNodes, 3], dtype=np.float64)

        # twoD mesh's boundingBox: [xmin, ymin, zmin, xmax, ymax, zmax]
        self.twoDMeshBoundingbox = []

        # each element (cell)'s bed elevation (z)
        self.elementBedElevation = np.zeros(self.numOfElements, dtype=np.float64)

        # each NodeString's list of nodes (stored in a dictionary)
        self.nodeStringsDict = {}

        # get the mesh information (elementNodesList, elementNodesCount, vtkCellTypeCode, and nodeCoordinates)
        # from reading the SRHGEOM file again
        self.readSRHGEOMFile()

        # list of elements for all nodes
        self.nodeElementsList = np.zeros([self.numOfNodes, gMax_Elements_per_Node], dtype=int)

        # number of elements for each node (1,2,...,gMax_Elements_per_Node)
        self.nodeElementsCount = np.zeros(self.numOfNodes, dtype=int)

        # build the element list for each node
        self.buildNodeElements()

        #edges (connecting two nodes). Not using "faces" because they are lines. When 2D mesh is extruded to 3D,
        #these lines will become faces.
        self.edges = {} #dictionary: {[node1, node2]: edgeID}
        self.edges_r = {} #dictionary: {edgeID: [node1, node2]}  #revise of self.edges

        self.edgeElements = {} #dictionary: {edgeID: [element list]}. Here since one edge can be shared by two
                               #elements at the most, the element list's length is either 1 or 2

        #edges of each boundary. Only contains real boundaries, not monitoring lines.
        #boundaryID = nodeString ID read in from srhgeom file
        self.boundaryEdges = {} #dictionary: {boundaryID: [list of edge IDs]}

        #list of all boundary edge IDs (all lumped to one list)
        self.allBoundaryEdgeIDs = []

        #build "lines" and "boundaryLines" dictionaries. SRH-2D has no such information in SRHGEOM. We need to build it.
        self.buildEdgesAndBoundaryEdges()


    def getNumOfElementsNodes(self):
        """ Get the number of elements and nodes in srhgeom mesh file

        Returns
        -------

        """
        print("Getting numbers of elements and nodes from the SRHGEOM file ...")

        # read the "srhgeom" mesh file
        try:
            srhgeomfile = open(self.srhgeom_filename, 'r')
        except:
            print('Failed openning the SRHGEOM file', self.srhgeom_filename)
            sys.exit()

        count = 0
        elemCount = 0
        nodeCount = 0
        nodeStringCount = 0

        while True:
            count += 1

            # Get next line from file
            line = srhgeomfile.readline()

            # if line is empty
            # end of file is reached
            if not line:
                break

            # print("Line{}: {}".format(count, line.strip()))

            search = line.split()
            # print(search)

            if len(search) != 0:
                if search[0] == "Elem":
                    elemCount += 1
                    # print("Elem # %d: %s" % (elemCount, line))
                elif search[0] == "Node":
                    nodeCount += 1
                    # print("Node # %d: %s" % (nodeCount, line))
                elif search[0] == "NodeString":
                    nodeStringCount += 1

        srhgeomfile.close()

        self.numOfElements = elemCount
        self.numOfNodes = nodeCount
        self.numOfNodeStrings = nodeStringCount

        print("There are %d elements, %d nodes, and %d node strings in the mesh." % (self.numOfElements,
                                                                                    self.numOfNodes, self.numOfNodeStrings))


    def readSRHGEOMFile(self):
        """ Get mesh information by reading the srhgeom file

        Parameters
        ----------

        Returns
        -------

        """

        print("Reading the SRHGEOM file ...")

        # read the "srhgeom" mesh file
        try:
            srhgeomfile = open(self.srhgeom_filename, 'r')
        except:
            print('Failed openning srhgeom file', self.srhgeom_filename)
            sys.exit()

        count = 0
        elemCount = 0
        nodeCount = 0
        nodeStringCound = 0

        # nodeString could be long enough to have more than one line. Use this to record the current nodeString during reading.
        currentNodeStringID = -1

        while True:
            count += 1

            # Get next line from file
            line = srhgeomfile.readline()

            # if line is empty
            # end of file is reached
            if not line:
                break

            # print("Line{}: {}".format(count, line.strip()))

            search = line.split()
            # print(search)

            if len(search) != 0:
                if search[0] == "Elem":
                    elemCount += 1
                    # print("Elem # %d: %s" % (elemCount, line))
                    self.elementNodesList[elemCount - 1][0:len(search[2:])] = [int(i) for i in search[2:]]
                    self.elementNodesCount[elemCount - 1] = len(search[2:])
                    if len(search[2:]) < 1 or len(search[2:]) > gMax_Nodes_per_Element:
                        sys.exit("Number of nodes for element %d is less than 1 or larger than the max of %d." % (elemCount,gMax_Nodes_per_Element))
                    self.vtkCellTypeCode[elemCount - 1] = vtkCellTypeMap[len(search[2:])]
                elif search[0] == "Node":
                    nodeCount += 1
                    # print("Node # %d: %s" % (nodeCount, line))
                    self.nodeCoordinates[nodeCount - 1] = [float(i) for i in search[2:]]
                elif search[0] == "NodeString":
                    nodeStringCound += 1
                    self.nodeStringsDict[int(search[1])] = [int(i) for i in search[2:]]
                    currentNodeStringID = int(search[1])
                elif ("SRHGEOM".lower() not in search[0].lower()) and \
                     ("Name".lower() not in search[0].lower()) and \
                     ("Grid".lower() not in search[0].lower()): #assume this is still nodeString
                    node_list = self.nodeStringsDict[currentNodeStringID]
                    for i in search:
                        node_list.append(int(i))
                    self.nodeStringsDict[currentNodeStringID] = node_list

        srhgeomfile.close()

        #calculate the bed elevation at cell center
        for cellI in range(self.numOfElements):

            elev_temp = 0.0

            #loop over all nodes of current element
            for nodeI in range(self.elementNodesCount[cellI]):
                elev_temp += self.nodeCoordinates[self.elementNodesList[cellI][nodeI]-1,2]  #z elevation; nodeI-1 because node number is 1-based for SRH-2D

            self.elementBedElevation[cellI] = elev_temp / self.elementNodesCount[cellI]

        #calculate the 2D mesh's bounding box
        xmin = np.min(self.nodeCoordinates[:,0])
        ymin = np.min(self.nodeCoordinates[:,1])
        zmin = np.min(self.nodeCoordinates[:,2])
        xmax = np.max(self.nodeCoordinates[:,0])
        ymax = np.max(self.nodeCoordinates[:,1])
        zmax = np.max(self.nodeCoordinates[:,2])

        self.twoDMeshBoundingbox = [xmin, ymin, zmin, xmax, ymax, zmax]

        print("2D mesh's bounding box = ", self.twoDMeshBoundingbox)

        if False:
            print("elementNodesList = ", self.elementNodesList)
            print("elementNodesCount = ", self.elementNodesCount)
            print("vtkCellTypeCode = ", self.vtkCellTypeCode)
            print("nodeCoordinates = ", self.nodeCoordinates)
            print("elementBedElevation = ", self.elementBedElevation)
            print("nodeStrings = ", self.nodeStringsDict)

    def buildNodeElements(self):
        """ Build node's element list for all nodes

        Returns
        -------

        """

        print("Building mesh's node, elements, and topology ...")

        #create an emptp list of lists for all nodes
        self.nodeElementsList = [[] for _ in range(self.numOfNodes)]

        #loop over all cells
        for cellI in range(self.numOfElements):
            #loop over all nodes of current cell
            for i in range(self.elementNodesCount[cellI]):
                #get the node ID of current node
                nodeID = self.elementNodesList[cellI][i]

                self.nodeElementsList[nodeID-1].append(cellI)  #nodeID-1 because SRH-2D is 1-based

        #count beans in each node's element list
        for nodeI in range(self.numOfNodes):
            self.nodeElementsCount[nodeI] = len(self.nodeElementsList[nodeI])

        #print("nodeElementsCount = ", self.nodeElementsCount)
        #print("nodeElementsList = ", self.nodeElementsList)

    def buildEdgesAndBoundaryEdges(self):
        """ Build edges and boundaryEdges dictionaries

        Returns
        -------

        """

        current_edgeID = 1  #note: SRH-2D is 1-based.

        #loop over all elements
        for cellI in range(self.numOfElements):
            # loop over all edges of current element
            for i in range(self.elementNodesCount[cellI]):
                # get the node ID of current and next nodes
                # connecting current edge

                nodeID_1 = 0
                nodeID_2 = 0

                if i != (self.elementNodesCount[cellI]-1):  #if not the last edge
                    nodeID_1 = self.elementNodesList[cellI][i]
                    nodeID_2 = self.elementNodesList[cellI][i+1]
                else:                                      #if the last edge
                    nodeID_1 = self.elementNodesList[cellI][i]
                    nodeID_2 = self.elementNodesList[cellI][0]

                curr_edge_node_IDs = tuple(sorted([nodeID_1, nodeID_2]))

                #check whether the current edge is in the edges dictionary or not
                if curr_edge_node_IDs not in self.edges:  #if no, added a new edge
                    self.edges[curr_edge_node_IDs] = current_edgeID
                    self.edges_r[current_edgeID] = curr_edge_node_IDs

                    #add the new edge to edgeElements dictionary too
                    self.edgeElements[current_edgeID] = [cellI+1]  #+1 because SRH-2D is 1-based

                    current_edgeID += 1
                else: # if yes, neighbor element already has that edge and no need to do anything, but need to add
                      # the current element to edgeElements list

                    #get the existing edge's ID
                    existingEdgeID = self.edges[curr_edge_node_IDs]

                    #in this case, there should be one and only one element in the list
                    assert len(self.edgeElements[existingEdgeID]) == 1

                    self.edgeElements[existingEdgeID].append(cellI+1) #+1 because SRH-2D is 1-based


        #build the allBoundaryEdgeIDs list: boundary edges are those with only one element
        #loop over all edges in the mesh
        for edge in self.edgeElements:
            if len(self.edgeElements[edge]) == 1:
                self.allBoundaryEdgeIDs.append(edge)

        #list to flag whether a boundary edge in self.allBoundaryEdgeIDs is used or not
        #if not, it is a default boundary (wall)
        allBoundaryEdgeUsageFlag = {}
        for boundaryEdgeID in self.allBoundaryEdgeIDs:
            allBoundaryEdgeUsageFlag[boundaryEdgeID] = False

        #now all edges are built, we check the boundary edges
        #loop through all boundaries.
        for nodeString in self.nodeStringsDict:

            #not all NodeStrings are boundaries. Could be monitoring lines. Need to exclude them.
            if nodeString not in self.bcDict.keys():
                continue

            #list of nodes in current nodeString
            nodeString_nodeList = self.nodeStringsDict[nodeString]

            current_boundary_edge_list = []

            #loop through all edges in current boundary
            for i in range(len(nodeString_nodeList)-1):
                nodeID_1 = nodeString_nodeList[i]
                nodeID_2 = nodeString_nodeList[i + 1]

                curr_edge_node_IDs = tuple(sorted([nodeID_1, nodeID_2]))

                #check whether the edge is in the edge list. If not, we got a problem.
                if curr_edge_node_IDs not in self.edges:
                    print("Boundary edge ", curr_edge_node_IDs, "in NodeString", nodeString, "can not be found in edge list. Mesh is wrong. Exiting...")
                    sys.exit()

                current_boundary_edge_list.append(self.edges[curr_edge_node_IDs])

                allBoundaryEdgeUsageFlag[self.edges[curr_edge_node_IDs]] = True

            self.boundaryEdges[nodeString] = current_boundary_edge_list

        #build default boundary as wall (in SMS, the default boundary is not exported in SRHGEOM)
        #boundary edges are those who is used only by one element.
        unusedBoundaryEdgeList = []
        for boundaryEdgeID in self.allBoundaryEdgeIDs:
            if not allBoundaryEdgeUsageFlag[boundaryEdgeID]: #if not used
                unusedBoundaryEdgeList.append(boundaryEdgeID)

        defaultWallNodeStringID = self.numOfNodeStrings+1   #nodeString ID for default boundary
        self.boundaryEdges[defaultWallNodeStringID] = unusedBoundaryEdgeList

        #debug
        if False:
            print("edges", self.edges)
            print("edgeElements",self.edgeElements)
            print("allBoundaryEdgeIDs",self.allBoundaryEdgeIDs)
            print("boundaryEdges", self.boundaryEdges)

    def extrude_to_3D(self, layerHeights, mshFileName, bTerrainFollowing=False, dir=''):
        """ Extrude the 2D mesh to 3D with layers

        When extruded: node -> line, edge -> face, 2D element-> volume, boundary line -> boundary face



        Also writes to a GMSH MSH file. Tried MSH version 4. But it seems too complicated (entities etc.). They are not
        necessary for our purpose. Therefore, we will stick with version 2 and OpenFOAM takes this version well.
        https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029

        Attributes
        -------
        layerHeights: a list with strictly increasing heights for the layers. If bTerrainFollowing is True, only
                      the last value is used as the top elevation.
        mshFileName: name of the GMSH MSH file to be written
        bTerrainFollowing: whether to follow the terrain. If False, the bottom will be at z = 0, and the layer heights
                           are specified in layerHeights. If True, the bottom will key the terrain elevation, the top
                           will be at layerHeights[-1] and interpolation inbetween.

        Returns
        -------

        """

        #GMSH element type code
        MSHTRI = 2  # 3-node triangle
        MSHQUAD = 3  # 4-node quadrilateral

        MSHHEX = 5  # 8-node hexahedron
        MSHPRISM = 6  # 6-node prism

        allNodes = {}   #dict: {nodeID: [x,y,z]}
        allCells = {}   #dict: {cellID: [cell_type, [list of nodes]]}
        number_of_prism = 0 #number of prism in the extruded 3D mesh
        number_of_hex = 0 #number of hex in the extruded 3D mesh

        #some sanity check nlayers and layerHeights
        nlayers = len(layerHeights)

        if nlayers == 0:
            print("layerHeight is empty. Check. Exiting ...")
            sys.exit()

        #check the values in the list layerHeights are strictly increasing
        bHeighIncreasing = all(i < j for i, j in zip(layerHeights, layerHeights[1:]))

        if not bHeighIncreasing:
            print("Values in layerHeight are not strictly increasing. Check. Exiting ...")
            sys.exit()

        #a list for all side boundaries (not including top and bottom)
        allSideBoundaryFaces = {} #dict: {boundaryID: [list of faces]}. Here face is a list of nodes.

        #create all nodes
        #1. add the original nodes in 2D mesh
        for nodeI in range(self.nodeCoordinates.shape[0]):
            if bTerrainFollowing:
                allNodes[nodeI] = self.nodeCoordinates[nodeI,:]
            else:
                allNodes[nodeI] = self.nodeCoordinates[nodeI, :]
                allNodes[nodeI][2] = 0.0  # make the 2D mesh flat

        #2. add the extruded nodes
        for layerI in range (1, nlayers+1):
            for nodeI in range(self.nodeCoordinates.shape[0]):
                if bTerrainFollowing: #interpolation between bottom elevation and the top elevation in layerHeights[-1]
                    elev = allNodes[nodeI][2] + (layerHeights[-1] - allNodes[nodeI][2])/nlayers*layerI
                    allNodes[nodeI + layerI * self.numOfNodes] = np.array([self.nodeCoordinates[nodeI, 0],
                                                                           self.nodeCoordinates[nodeI, 1],
                                                                           elev])
                else:
                    # make the original 2D mesh flat; only add the layer height
                    allNodes[nodeI + layerI*self.numOfNodes] = np.array([self.nodeCoordinates[nodeI,0],
                                                                         self.nodeCoordinates[nodeI,1],
                                                                         layerHeights[layerI-1]])

        #create all cells

        #bottom and top boundary faces list
        bottomBoundaryFaces = [None] * self.numOfElements
        topBoundaryFaces = [None] * self.numOfElements

        cellID = 0

        #layer by layer
        for layerI in range(1, nlayers+1):
            #loop through each element in 2D mesh
            for elementI in range(self.numOfElements):

                if self.elementNodesCount[elementI] == 3: #triangle -> prism
                    cell_type = MSHPRISM
                    number_of_prism += 1
                elif self.elementNodesCount[elementI] == 4: #quad -> hex
                    cell_type = MSHHEX
                    number_of_hex += 1

                cellID += 1

                cell_list = []

                top_face_node_list = []
                bottom_face_node_list = []

                #loop over all nodes of the current elment in 2D mesh (at current layer's bottom)
                for nodeI in range(self.elementNodesCount[elementI]):
                    cell_list.append(self.elementNodesList[elementI, nodeI] + (layerI-1)*self.numOfNodes)

                    #if this is the first layer, record the bottom face to bottomBoundaryFaces list
                    if layerI == 1:
                        bottom_face_node_list.append(self.elementNodesList[elementI, nodeI])

                #now add the nodes at current layer's top.
                for nodeI in range(self.elementNodesCount[elementI]):
                    cell_list.append(self.elementNodesList[elementI, nodeI] + layerI*self.numOfNodes)

                    #if this is the last layer, record the top face to topBoundaryFaces list
                    if layerI == nlayers:
                        top_face_node_list.append(self.elementNodesList[elementI, nodeI] + layerI*self.numOfNodes)

                allCells[cellID] = [cell_type, cell_list]

                bottomBoundaryFaces[elementI] = bottom_face_node_list
                topBoundaryFaces[elementI] = top_face_node_list


        #create all boundaries

        #allSideBoundaryFaces = {} #dict: {boundaryID: [list of faces]}. Here face is a list of nodes.

        #loop through all boundary lines (nodeStrings) in 2D mesh. self.boundaryEdges only
        #contains real boundaries, not monitoring lines.
        for boundaryID in self.boundaryEdges:

            #list of all faces for current boundary
            #each face is a list four nodes (quad because we extrude a line in vertical direction)
            faceList = []

            #loop through all edges of current boundary
            for edgeID in self.boundaryEdges[boundaryID]:
                nodeID_1 = self.edges_r[edgeID][0]
                nodeID_2 = self.edges_r[edgeID][1]

                #loop through all layers
                for layerI in range(1, nlayers + 1):
                    #add the four nodes to the current face (ordered)
                    curFace = [nodeID_1+(layerI-1)*self.numOfNodes, nodeID_1+layerI*self.numOfNodes,
                               nodeID_2+ layerI*self.numOfNodes, nodeID_2+(layerI-1)*self.numOfNodes]

                    #add the current face to the faceList for current boundary
                    faceList.append(curFace)

            #now we have all faces in the current boundary
            allSideBoundaryFaces[boundaryID] = faceList


        #write to GMSH MSH file

        if dir!='':
            mshFileName_base = dir + '/' + mshFileName
        else:
            mshFileName_base = mshFileName
        print("Write to GMESH MSH file with name = ", mshFileName_base)


        fid = open(mshFileName_base, 'w')

        #write MeshFormat
        fid.write("$MeshFormat\n")
        fid.write("2.1 0 8\n")
        fid.write("$EndMeshFormat\n")

        #write PhysicalNames:
        fid.write("$PhysicalNames\n")
        fid.write("%d\n" % (len(allSideBoundaryFaces)+2+1))   #need to output all side boundaries,
                                                              #top + bottom (2), and the volume (1)

        # maximum side boundary ID. This is the starting point for IDs of top, bottom, and volume
        maxSideBoundaryID = -999

        #loop through all side boundaries
        for boundaryID in allSideBoundaryFaces:
            if boundaryID in self.bcDict.keys(): #if the current BC is defined in srhhydro
                boundary_name = "boundary_"+str(boundaryID)+"_"+self.bcDict[boundaryID]  #name = boundary_ID_BCType
            else:                                #else, this should be default wall boundary
                boundary_name = "boundary_"+str(boundaryID)+"_"+"wall"                   #name = boundary_ID_wall

            fid.write("2 %d \"%s\"\n" % (boundaryID, boundary_name))  #dimension, physicalTag, name

            maxSideBoundaryID = max(maxSideBoundaryID, boundaryID)

        #output the top and bottom boundaries
        boundary_name = "top"
        fid.write("2 %d \"%s\"\n" % (maxSideBoundaryID+1, boundary_name))  # dimension, physicalTag, name

        boundary_name = "bottom"
        fid.write("2 %d \"%s\"\n" % (maxSideBoundaryID+2, boundary_name))  # dimension, physicalTag, name


        #output the volume
        volume_name = "channel"  #just a generic name for the volume
        fid.write("3 %d \"%s\"\n" % (maxSideBoundaryID+2+1, volume_name))

        fid.write("$EndPhysicalNames\n")

        #write Nodes
        fid.write("$Nodes\n")
        #total_num_nodes = (nlayers+1)*self.numOfNodes
        total_num_nodes = len(allNodes)

        fid.write("%d\n" % total_num_nodes)

        #output node coordinates
        for nodeI in allNodes:
            #GMSH is 1-based
            fid.write("%d %f %f %f\n" %(nodeI+1, allNodes[nodeI][0], allNodes[nodeI][1], allNodes[nodeI][2]))

        fid.write("$EndNodes\n")

        #write Elements (boundary conditions and volume)
        #  numElements = number of all faces in all boundaries + number of cells in all volumes

        elementTag_counter = 0 #count the number of elements

        fid.write("$Elements\n")

        if (number_of_hex == 0) and (number_of_prism == 0):
            print("There is on hex or prism, the only supported cell types, in the extruded 3D mesh. Check the mesh. Exiting ...")
            sys.exit()

        #output total number of elements
        #                    num. of side boundary faces + num. of faces on top and bottom + num. of hex + num. of prism
        fid.write("%d \n" % (len(self.allBoundaryEdgeIDs)*nlayers + self.numOfElements*2 + number_of_hex + number_of_prism))

        #loop through all side boundaries
        for boundaryID in allSideBoundaryFaces:
            #output each quad face's
            for faceID in range(len(allSideBoundaryFaces[boundaryID])):
                elementTag_counter += 1
                fid.write("%d 3 1 %d %d %d %d %d\n" % (elementTag_counter, boundaryID, allSideBoundaryFaces[boundaryID][faceID][0],
                                    allSideBoundaryFaces[boundaryID][faceID][1], allSideBoundaryFaces[boundaryID][faceID][2],
                                    allSideBoundaryFaces[boundaryID][faceID][3]))

        #output top and bottom boundaries
        for faceI in topBoundaryFaces:
            elementTag_counter += 1

            if len(faceI) == 3: #triangle
                fid.write("%d 2 1 %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID+1,
                                                     faceI[0],faceI[1],faceI[2]))
            elif len(faceI) == 4: #quad
                fid.write("%d 3 1 %d %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID+1,
                                                     faceI[0], faceI[1], faceI[2], faceI[3]))

        for faceI in bottomBoundaryFaces:
            elementTag_counter += 1

            if len(faceI) == 3:  # triangle
                fid.write("%d 2 1 %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID + 2,
                                                     faceI[0], faceI[1], faceI[2]))
            elif len(faceI) == 4:  # quad
                fid.write("%d 3 1 %d %d %d %d %d \n" % (elementTag_counter, maxSideBoundaryID + 2,
                                                        faceI[0], faceI[1], faceI[2], faceI[3]))

        #output the volume's tag and node list. Hex and prism need to be seperate
        #1.Hex
        if number_of_hex > 0:
            #loop through all cells in the volume
            for cellID in allCells:
                if allCells[cellID][0] == MSHHEX:
                    elementTag_counter += 1
                    fid.write("%d 5 1 %d %d %d %d %d %d %d %d %d\n" % (elementTag_counter, maxSideBoundaryID+2+1,
                                                                       allCells[cellID][1][0], allCells[cellID][1][1],
                                            allCells[cellID][1][2], allCells[cellID][1][3], allCells[cellID][1][4],
                                            allCells[cellID][1][5], allCells[cellID][1][6], allCells[cellID][1][7]))

        #2. prism
        if number_of_prism > 0:
            #loop through all cells in the volume
            for cellID in allCells:
                if allCells[cellID][0] == MSHPRISM:
                    elementTag_counter += 1
                    fid.write("%d 6 1 %d %d %d %d %d %d %d\n" % (elementTag_counter, maxSideBoundaryID+2+1,
                                                                 allCells[cellID][1][0], allCells[cellID][1][1],
                                                                 allCells[cellID][1][2], allCells[cellID][1][3],
                                                                 allCells[cellID][1][4], allCells[cellID][1][5]))

        fid.write("$EndElements")

        fid.close()


    def output_2d_mesh_to_vtk(self,flatMeshVTKFileName, bFlat=False, dir=''):
        """output the 2D mesh to vtk


        Attributes
        -------
        bFlat: whether to make a flat mesh (z=0)

        Returns
        -------

        """

        vtkFileName = ''
        if len(dir) == 0:
            vtkFileName = flatMeshVTKFileName
        else:
            vtkFileName = dir + "/" + flatMeshVTKFileName

        # build VTK object:
        # points
        pointsVTK = vtk.vtkPoints()
        flatCoordinates = np.copy(self.nodeCoordinates)

        if bFlat:
            flatCoordinates[:,2] = 0.0

        pointsVTK.SetData(VN.numpy_to_vtk(flatCoordinates))

        # cell topology information list: [num. of nodes, node0, node1, .., num. of nodes, nodexxx]
        # the list start with the number of nodes for a cell and then the list of node indexes
        connectivity_list = []
        # type of cells (contains the number of face points
        cellFPCounts = np.zeros(self.elementNodesList.shape[0], dtype=np.int64)

        # loop over all elements
        for k in range(self.elementNodesList.shape[0]):
            connectivity_list.append(self.elementNodesCount[k])

            for nodeI in range(self.elementNodesCount[k]):
                connectivity_list.append(
                    self.elementNodesList[k][nodeI] - 1)  # -1 becasue SRH-2D is 1-based.

            cellFPCounts[k] = self.elementNodesCount[k]

        connectivity = np.array(connectivity_list, dtype=np.int64)

        # convert cell's number of face points to VTK cell type
        vtkHandler_obj = vtkHandler()
        cell_types = vtkHandler_obj.number_of_nodes_to_vtk_celltypes(cellFPCounts)

        cellsVTK = vtk.vtkCellArray()
        cellsVTK.SetCells(self.elementNodesList.shape[0], VN.numpy_to_vtkIdTypeArray(connectivity))

        uGrid = vtk.vtkUnstructuredGrid()
        uGrid.SetPoints(pointsVTK)
        uGrid.SetCells(cell_types, cellsVTK)

        # write to vtk file
        unstr_writer = vtk.vtkUnstructuredGridWriter()
        unstr_writer.SetFileName(vtkFileName)
        unstr_writer.SetInputData(uGrid)
        unstr_writer.Write()


    def output_nodeString_line_coordinates(self, nodeStringID, nodeStringFileName, dir=''):
        """ Output the nodeString line coordinates into a text file


        Parameters
        ----------
        nodeStringID: {integer} -- the ID of the nodeString to be written out
        nodeStringFileName: {string} -- the name of the file to be written to
        dir: {string} -- optional directory name

        Returns
        -------

        """

        nodeStringFileName_final = ''
        if len(dir) == 0:
            nodeStringFileName_final = nodeStringFileName
        else:
            nodeStringFileName_final = dir + "/" + nodeStringFileName

        #check whether the given nodeStringID is in the nodeStringDict
        if nodeStringID not in self.nodeStringsDict.keys():
            print("The given nodeStringID", nodeStringID, "is not valid. Valid nodeString IDs are: ",
                  self.nodeStringsDict.keys())

        try:
            fid = open(nodeStringFileName_final, 'w')
        except IOError:
            print('NodeString file open error')
            sys.exit()

        #write the header
        fid.write("station, x, y, z\n")

        # list of nodes in current nodeString
        nodeString_nodeList = self.nodeStringsDict[nodeStringID]

        station = [0] * len(nodeString_nodeList)
        x = [0] * len(nodeString_nodeList)
        y = [0] * len(nodeString_nodeList)
        z = [0] * len(nodeString_nodeList)

        #loop over all nodes in the list
        for nodeI in range(len(nodeString_nodeList)):
            x[nodeI] = self.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 0]
            y[nodeI] = self.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 1]
            z[nodeI] = self.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 2]

            if nodeI!=0: #if not the first node in the list, calculate the station
                station[nodeI] = station[nodeI-1] + np.sqrt( (x[nodeI]-x[nodeI-1])**2 + (y[nodeI]-y[nodeI-1])**2 )

        #write the data
        for nodeI in range(len(nodeString_nodeList)):
            fid.write("%f, %f, %f, %f\n" % (station[nodeI], x[nodeI], y[nodeI], z[nodeI]))

        fid.close()

    def output_as_gmsh(self, gmshFileName, dir=''):
        """ Output SRH-2D mesh in GMESH format

        Parameters
        ----------
        dir : str, optional
            directory name

        Returns
        ---------

        """

        return NotImplemented

        if dir!='':
            gmshFileName_base = dir + '/' + gmshFileName
        else:
            gmshFileName_base = gmshFileName
        print("GMESH file name = ", gmshFileName_base)




class SRH_2D_SRHMat:
    """A class to handle srhmat file for SRH-2D

    Attributes:


    Methods:

    """
    def __init__(self, srhmat_filename):
        self.srhmat_filename = srhmat_filename

        #number of materials (Manning's n zones)
        self.numOfMaterials = -1

        #list of Material names (in a dictionary: ID and name)
        self.matNameList = {}

        #each material zone's cell list (in a dictionary: material zone ID and cell list)
        self.matZoneCells = {}

        #build material zones data
        self.buildMaterialZonesData()

    def buildMaterialZonesData(self):
        """Build the data for material zones

        Returns
        -------

        """
        print("Reading the SRHMAT file ...")

        # read the "srhmat" material file
        try:
            srhmatfile = open(self.srhmat_filename, 'r')
        except:
            print('Failed openning SRHMAT file', self.srhmat_filename)
            sys.exit()

        res_MatNames = {} #dict to hold "MatName" entries for material names list: ID and name
        res_Materials = {} #dict to hold "Material" information: ID and list of cells

        current_MaterialID = -1
        current_MaterialCellList = []

        while True:
            # Get a line from file
            line = srhmatfile.readline()

            # if line is empty
            # end of file is reached
            if not line:
                break

            # print("Line{}: {}".format(count, line.strip()))

            search = line.split()
            # print(search)

#            print("search[0] == Material", (search[0] == "Material"))

            if len(search) != 0: #if it is not just an empty line
                if search[0] == "SRHMAT":
                    continue
                elif search[0] == "NMaterials":
                    self.numOfMaterials = int(search[1])
                elif search[0] == "MatName":
                    res_MatNames[search[1]] = search[2]
                elif search[0] == "Material": # a new material zone starts
                    #if this is not the first Material zone; save the previous Material zone
                    if ((current_MaterialID != -1) and (len(current_MaterialCellList)!=0)):
                        res_Materials[current_MaterialID] = copy.deepcopy(current_MaterialCellList)

                        #clear up current_MaterialID and current_MaterialCellList
                        current_MaterialID = -1
                        current_MaterialCellList.clear()

                    current_MaterialID = int(search[1])
                    current_MaterialCellList.extend(search[2:])
                else: #still in a material zone (assume there is no other things other than those listed above in the
                      #srhmat file)
                    current_MaterialCellList.extend(search)

            #add the last material zone
            res_Materials[current_MaterialID] = copy.deepcopy(current_MaterialCellList)

        srhmatfile.close()

        #convert material ID from str to int
        for zoneID, cellList in res_Materials.items():
            #print("cellList ", cellList)
            temp_list = [int(i) for i in cellList]
            res_Materials[zoneID] = temp_list

        self.matNameList = res_MatNames
        self.matZoneCells = res_Materials

        #print("matNameList = ", self.matNameList)
        #print("matZoneCells = ", self.matZoneCells)


    def find_cell_material_ID(self, cellID):
        """ Given a cell ID, find its material (Manning's n) ID

        Parameters
        ----------
        cellID:

        Returns
        -------

        """

        if not self.matZoneCells:
            self.buildMaterialZonesData()

        #whether the cell is found in the list
        bFound = False

        for matID, cellList in self.matZoneCells.items():
            if (cellID+1) in cellList:  #cellID+1 because SRH-2D is 1-based
                bFound = True
                return matID

        #if cellID is not found, report problem.
        #This could be due to several reasons. Most likely, the material coverage in SMS does not
        #cover certain area of the mesh. In this case, SMS will not report the material for the missed
        #cells in srhmat (no warning given either). Here, if we can't find the cell in material list,
        #simply use the default Manning's n ID.
        if not bFound:
            print("In find_cell_material_ID(cellID), cellID =", cellID, "is not found. Check mesh and material coverage. Default is used.")
            return 0 #return the default (?)
            #sys.exit()


class SRH_2D_Data:
    """
    A class for SRH-2D data I/O, manipulation, and format conversion
    
    This class is designed to read SRH-2D results in format other than VTK. It can 
    save SRH-2D results into VTK format for visualization in Paraview, parse
    SRH-2D mesh information and convert/save to other formats, query SRH-2D
    results (e.g., for calibration), etc.
    
    Attributes
    ----------
    hdf_filename : str
        The name for the HDF result file generated by HEC-RAS
    srhhydro_obj : SRH_2D_SRHHydro
        An object to hold information in the SRHHYDRO file
    srhgeom_obj : SRH_2D_SRHGeom
        An object to hold information in the SRHGEOM file
    srhmat_obj : SRH_2D_SRHMat
        An object to hold information in the SRHMAT file

    Methods
    -------

    
    """
    
    def __init__(self, srhhydro_filename):
        """Constructor for SRH_2D_Data

        srhgeom_filename and srhmat_filename are contained in the SRHHydro file.

        Parameters
        ----------
        srhhydro_filename : str
            Name of the SRHHydro file

        """

        self.srhhydro_filename = srhhydro_filename

        #extract path in srhhydro_filename. We assume the SRHGeom and SRHMat files
        #are located in the same directory as the SRHHydro file. In the SRHHydro file,
        #the file names for SRHGeom and SRHMat do not contain the directory.
        file_path, _, = os.path.split(srhhydro_filename)

        #read SRH_2D_SRHHydro file and build SRH_2D_SRHHydro object
        self.srhhydro_obj = SRH_2D_SRHHydro(self.srhhydro_filename)

        #get the srhgeom_filename and srhmat_filename from srhhydro_obj
        self.srhgeom_filename =      self.srhhydro_obj.get_grid_file_name() if len(file_path) == 0 \
                                else file_path+"/"+self.srhhydro_obj.get_grid_file_name()

        self.srhmat_filename =       self.srhhydro_obj.get_mat_file_name() if len(file_path) == 0 \
                                else file_path+"/"+self.srhhydro_obj.get_mat_file_name()

        #read and build SRH_2D_SRHHydro, SRH_2D_Geom, and SRH_2D_Mat objects
        self.srhhydro_obj = SRH_2D_SRHHydro(self.srhhydro_filename)
        self.srhgeom_obj = SRH_2D_SRHGeom(self.srhgeom_filename,self.srhhydro_obj.srhhydro_content["BC"])
        self.srhmat_obj = SRH_2D_SRHMat(self.srhmat_filename)

        #Manning's n value at cell centers and nodes
        self.ManningN_cell = np.zeros(self.srhgeom_obj.numOfElements, dtype=float)
        self.ManningN_node = np.zeros(self.srhgeom_obj.numOfNodes, dtype=float)

        #build Manning's n value at cell centers and nodes
        self.buildManningN_cells_nodes()

        #XMDF data: if desired, the result data in a XMDF file
        #can be read in and stored in the following arrays.
        #Depending on whether the XMDF is nodal or at cell centers, use different arrays.
        self.xmdfTimeArray_Nodal = None   #numpy array to store all time values
        self.xmdfAllData_Nodal = {}     #dict to store {varName: varValues (numpy array)}

        self.xmdfTimeArray_Cell = None   #numpy array to store all time values
        self.xmdfAllData_Cell = {}     #dict to store {varName: varValues (numpy array)}

    def get_case_name(self):
        """Get the "Case" name for this current case ("Case" in srhhyro file)

        Returns
        -------
        case_name : str
            Case name

        """
        return self.srhhydro_obj.srhhydro_content["Case"]

    def buildManningN_cells_nodes(self):
        """ Build Manning's n values at cell centers and nodes

        This calculation is based on the srhhydro and srhgeom files, not interpolation from a GeoTiff file. This is the
        Manning's n value used by SRH-2D.

        Returns
        -------

        """

        print("Building Manning's n values for cells and nodes in SRH-2D mesh ...")

        #Manning's n dictionary in srhhydro
        nDict = self.srhhydro_obj.srhhydro_content["ManningsN"]
        print("nDict = ", nDict)

        #loop over all cells in the mesh
        for cellI in range(self.srhgeom_obj.numOfElements):
            #get the material ID of current cell
            matID = self.srhmat_obj.find_cell_material_ID(cellI)

            #print("matID = ", matID)
            #print("nDict[matID] = ", nDict[matID])

            #get the Manning's n value for current cell
            self.ManningN_cell[cellI] = nDict[matID]

        #interpolate Manning's n from cell centers to nodes
        #loop over all nodes in the mesh
        for nodeI in range(self.srhgeom_obj.numOfNodes):
            n_temp = 0.0

            #loop over all cells that share the current node
            for i in range(self.srhgeom_obj.nodeElementsCount[nodeI]):
                cellI = self.srhgeom_obj.nodeElementsList[nodeI][i]
                n_temp += self.ManningN_cell[cellI]

            #take the average
            self.ManningN_node[nodeI] = n_temp/self.srhgeom_obj.nodeElementsCount[nodeI]

    def output_boundary_manning_n_profile(self, nodeStringID, nodeStringFileName, dir=''):
        """ Output Manning's n profile for a specified boundary

        Parameters
        ----------
        nodeStringID : int
            nodeString ID of specified boundary
        nodeStringFileName : str
            name of output file
        dir : str, optional
            directory for the output file

        Returns
        -------

        """

        nodeStringFileName_final = ''
        if len(dir) == 0:
            nodeStringFileName_final = nodeStringFileName
        else:
            nodeStringFileName_final = dir + "/" + nodeStringFileName

        #check whether the given nodeStringID is in the nodeStringDict
        if nodeStringID not in self.srhgeom_obj.nodeStringsDict.keys():
            print("The given nodeStringID", nodeStringID, "is not valid. Valid nodeString IDs are: ",
                  self.srhgeom_obj.nodeStringsDict.keys())

        try:
            fid = open(nodeStringFileName_final, 'w')
        except IOError:
            print('Boundary\'s Manning n file open error')
            sys.exit()

        #write the header
        fid.write("station, n\n")

        # list of node ID (1-based) in current nodeString
        nodeString_nodeList = self.srhgeom_obj.nodeStringsDict[nodeStringID]

        station = [0] * len(nodeString_nodeList)
        x = [0] * len(nodeString_nodeList)
        y = [0] * len(nodeString_nodeList)
        z = [0] * len(nodeString_nodeList)

        #loop over all nodes in the list
        for nodeI in range(len(nodeString_nodeList)):
            x[nodeI] = self.srhgeom_obj.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 0]
            y[nodeI] = self.srhgeom_obj.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 1]
            z[nodeI] = self.srhgeom_obj.nodeCoordinates[nodeString_nodeList[nodeI] - 1, 2]

            if nodeI!=0: #if not the first node in the list, calculate the station
                station[nodeI] = station[nodeI-1] + np.sqrt( (x[nodeI]-x[nodeI-1])**2 + (y[nodeI]-y[nodeI-1])**2 )

        #write the data
        for nodeI in range(len(nodeString_nodeList)):
            fid.write("%f, %f\n" % (station[nodeI], self.ManningN_node[nodeString_nodeList[nodeI] - 1])) #-1 because SRH-2D is 1-based

        fid.close()


    def readSRHXMDFFile(self, xmdfFileName, bNodal):
        """ Read SRH-2D result file in XMDF format (current version 13.1.6 of SMS only support data at node).

        Parameters
        ----------
        xmdfFileName : str
            file name for the XMDF data
        bNodal : bool
            whether it is nodal (True) or at cell center (False)

        Returns
        -------


        """

        #for debug
        #dumpXMDFFileItems(xmdfFileName)

        if bNodal:
            if "XMDFC" in xmdfFileName:
                print("Warning: bNodal indices the result is for nodal values. However, the XMDF file name contains "
                      "\"XMDFC\", which indicates it is for cell center values. Please double check.")
        else:
            if "XMDFC" not in xmdfFileName:
                print("Warning: bNodal indices the result is for cell centers. However, the XMDF file name "
                      "does not contain \"XMDFC\", which indicates it is for nodal values. Please double check.")

        print("Reading the XMDF file ...\n")

        xmdfFile = h5py.File(xmdfFileName, "r")

        #build the list of solution variables
        varNameList = []

        #copy of the original velocity vector variable name
        varNameVelocity = ''

        for ds in xmdfFile.keys():
            #print(ds)

            if ds != "File Type" and ds != "File Version":
                if "Velocity" in ds:
                    varNameVelocity = '%s' % ds
                    vel_x = ds.replace("Velocity","Vel_X",1)
                    vel_y = ds.replace("Velocity","Vel_Y",1)
                    varNameList.append(vel_x)
                    varNameList.append(vel_y)
                else:
                    varNameList.append(ds)

        if not varNameList: #empty list
            print("There is no solution varialbes in the XMDF file. Exiting ...")
            sys.exit()

        print("Variables in the XMDF file: ", varNameList)

        #build the 1d array of time for the solution
        #It seems all solution variables share the same time (although
        #each has their own recrod of "Time"). So only get one is sufficient.
        if bNodal:
            self.xmdfTimeArray_Nodal = np.array(xmdfFile[varNameList[0]]['Times'])
            print("Time values for the solutions: ", self.xmdfTimeArray_Nodal)
        else:
            self.xmdfTimeArray_Cell = np.array(xmdfFile[varNameList[0]]['Times'])
            print("Time values for the solutions: ", self.xmdfTimeArray_Cell)

        #result for all varialbes and at all times in a dictionary (varName : numpy array)
        for varName in varNameList:
            print("varName =", varName)

            varName_temp = varName

            if ("Vel_X" in varName) or ("Vel_Y" in varName):  # if it is velocity components
                varName = varNameVelocity                     # use the original velocity variable name

            if bNodal: #for nodal data
                if "Velocity" in varName:
                    if "Vel_X" in varName_temp: #vel-x
                        #print("varName_temp = ", varName_temp)
                        self.xmdfAllData_Nodal[varName_temp] = np.array(xmdfFile[varName]['Values'][:,:,0])
                    else:
                        #print("varName_temp = ", varName_temp)
                        self.xmdfAllData_Nodal[varName_temp] = np.array(xmdfFile[varName]['Values'][:,:,1])
                else:
                    self.xmdfAllData_Nodal[varName] = np.array(xmdfFile[varName]['Values'])

                    #fixe the water elevation = -999 in SRH-2D
                    if varName == "Water_Elev_ft" or varName == "Water_Elev_m":
                        for nodeI in range(len(self.xmdfAllData_Nodal[varName])):
                            if self.xmdfAllData_Nodal[varName][nodeI] == -999:
                                self.xmdfAllData_Nodal[varName][nodeI] = self.srhgeom_obj.nodeCoordinates[nodeI, 2] #use node elevation


                if np.array(xmdfFile[varName]['Values']).shape[1] != self.srhgeom_obj.numOfNodes:
                    print("The number of nodes in the XMDF file (%d) is different from that in the mesh (%d). "
                          "Abort readSRHXMDFFile(...) function."
                          % (np.array(xmdfFile[varName]['Values']).shape[1], self.srhgeom_obj.numOfNodes))
                    return

                #print(self.xmdfAllData_Nodal)
            else: #for cell center data
                if "Velocity" in varName:
                    if "Vel_X" in varName_temp: #vel-x
                        #print("varName_temp, varName = ", varName_temp, varName)
                        self.xmdfAllData_Cell[varName_temp] = np.array(xmdfFile[varName]['Values'])[:,:,0]
                    else: #vel-y
                        #print("varName_temp = ", varName_temp)
                        self.xmdfAllData_Cell[varName_temp] = np.array(xmdfFile[varName]['Values'])[:,:,1]
                else:
                    self.xmdfAllData_Cell[varName] = np.array(xmdfFile[varName]['Values'])

                    #fixe the water elevation = -999 in SRH-2D
                    if varName == "Water_Elev_ft" or varName == "Water_Elev_m":
                        #loop through time
                        for timeI in range(self.xmdfAllData_Cell[varName].shape[0]):
                            for cellI in range(self.xmdfAllData_Cell[varName].shape[1]):
                                if self.xmdfAllData_Cell[varName][timeI, cellI] == -999:
                                    self.xmdfAllData_Cell[varName][timeI, cellI] = self.srhgeom_obj.elementBedElevation[cellI]

                if np.array(xmdfFile[varName]['Values']).shape[1] != self.srhgeom_obj.numOfElements:
                    print("The number of elements in the XMDF file (%d) is different from that in the mesh (%d). "
                          "Abort readSRHXMDFFile(...) function."
                          % (np.array(xmdfFile[varName]['Values']).shape[1], self.srhgeom_obj.numOfElements))

                #print(self.xmdfAllData_Cell)

        xmdfFile.close()

    def outputXMDFDataToVTK(self, bNodal, timeStep=-1,lastTimeStep=False, dir=''):
        """Output XMDF result data to VTK

        This has to be called after the XMDF data have been loaded by calling readSRHXMDFFile(...).

        It calls outputVTK(...).

        Parameters
        ----------
            bNodal : bool
                whether export nodal data or cell center data. Currently, it can't output both.
            timeStep : int
                optionly only the specified time step will be saved
            lastTimeStep : bool
                optionally specify only the last time step
            dir : str, optional
                directory to write to

        Returns
        -------

        """
        print("Output all data in the XMDF file to VTK ...")

        if (bNodal and (not self.xmdfAllData_Nodal)) or ((not bNodal) and (not self.xmdfAllData_Cell)):
            print("Empty XMDF data arrays. Call readSRHXMDFFile() function first. Exiting ...")
            sys.exit()

        #check the sanity of timeStep
        if (timeStep != -1):
            if bNodal:
                if (not timeStep in range(len(self.xmdfTimeArray_Nodal))):
                    message = "Specified timeStep = %d not in range (0 to %d)." % (timeStep, len(self.xmdfTimeArray_Nodal))
                    sys.exit(message)
            else:
                if (not timeStep in range(len(self.xmdfTimeArray_Cell))):
                    message = "Specified timeStep = %d not in range (0 to %d)." % (timeStep, len(self.xmdfTimeArray_Cell))
                    sys.exit(message)

        vtkFileName_base = ''
        #print("self.srhhydro_obj.srhhydro_content[\"Case\"] = ", self.srhhydro_obj.srhhydro_content["Case"])

        if dir!='':
            vtkFileName_base = dir + '/' + 'SRH2D_' + self.srhhydro_obj.srhhydro_content["Case"] #use the case name of the
            # SRH-2D case
        else:
            vtkFileName_base = 'SRH2D_' + self.srhhydro_obj.srhhydro_content["Case"]  # use the case name of the SRH-2D case
        print("vtkFileName_base = ", vtkFileName_base)

        #build result variable names
        resultVarNames = list(self.xmdfAllData_Nodal.keys()) if bNodal else list(self.xmdfAllData_Cell.keys())

        #add Manning's n (at cell center regardless bNodal)
        ManningNVarName = "ManningN"

        if bNodal:
            print("All nodal solution variable names: ", resultVarNames[:-1], "and cell center ManningN")
        else:
            print("All cell center solution variable names: ", resultVarNames)

        #add bed elevation at nodes regardless whether bNodal is True or False. Nodal elevation is more accurate representation
        #of the terrain because cell center elevation is averaged from nodal elevations.
        #get the units of the results
        units = self.srhhydro_obj.srhhydro_content['OutputFormat'][1]
        print("SRH-2D result units:", units)
        bedElevationVarName = ''
        velocityVarName = ''
        if units == "SI":
            bedElevationVarName = "Bed_Elev_m"
            velocityVarName = "Velocity_m_p_s"
        else:
            bedElevationVarName = "Bed_Elev_ft"
            velocityVarName = "Velocity_ft_p_s"

        print("Nodal bed elevation name: ", bedElevationVarName)

        #loop through each time step
        timeArray = self.xmdfTimeArray_Nodal if bNodal else self.xmdfTimeArray_Cell
        for timeI in range(timeArray.shape[0]):

            if lastTimeStep:
                if timeI < (timeArray.shape[0] - 1):
                    continue

            if (timeStep != -1) and (timeI != timeStep):
                continue

            print("timeI and time = ", timeI, timeArray[timeI])

            #numpy array for all solution variables at one time
            resultData =      np.zeros((self.srhgeom_obj.numOfNodes, len(resultVarNames)), dtype="float32") if bNodal \
                         else np.zeros((self.srhgeom_obj.numOfElements, len(resultVarNames)), dtype="float32")

            str_extra = "_N_" if bNodal else "_C_"

            vtkFileName = vtkFileName_base + str_extra + str(timeI+1).zfill(4) + ".vtk"
            print("vtkFileName = ", vtkFileName)

            #loop through each solution variable (except Bed_Elev and ManningN, which will be added seperately)
            for varName, varI in zip(resultVarNames, range(len(resultVarNames))):
                print("varName = ", varName)
                #get the values of current solution varialbe at current time
                resultData[:,varI] =      self.xmdfAllData_Nodal[varName][timeI,:] if bNodal \
                                     else self.xmdfAllData_Cell[varName][timeI,:]

            #add Mannng's n to cell center
            #self.ManningN_cell

            #add nodal bed elevation to VTK's point data
            nodalElevation = self.srhgeom_obj.nodeCoordinates[:,2]

            #build VTK object:
            # points
            pointsVTK = vtk.vtkPoints()
            pointsVTK.SetData(VN.numpy_to_vtk(self.srhgeom_obj.nodeCoordinates))

            # cell topology information list: [num. of nodes, node0, node1, .., num. of nodes, nodexxx]
            # the list start with the number of nodes for a cell and then the list of node indexes
            connectivity_list = []

            # type of cells (contains the number of face points
            cellFPCounts = np.zeros(self.srhgeom_obj.elementNodesList.shape[0], dtype=np.int64)

            #loop over all elements
            for k in range(self.srhgeom_obj.elementNodesList.shape[0]):
                connectivity_list.append(self.srhgeom_obj.elementNodesCount[k])

                for nodeI in range(self.srhgeom_obj.elementNodesCount[k]):
                    connectivity_list.append(self.srhgeom_obj.elementNodesList[k][nodeI]-1) #-1 becasue SRH-2D is 1-based.

                cellFPCounts[k] = self.srhgeom_obj.elementNodesCount[k]

            connectivity = np.array(connectivity_list, dtype=np.int64)

            # convert cell's number of face points to VTK cell type
            vtkHandler_obj = vtkHandler()
            cell_types = vtkHandler_obj.number_of_nodes_to_vtk_celltypes(cellFPCounts)

            cellsVTK = vtk.vtkCellArray()
            cellsVTK.SetCells(self.srhgeom_obj.elementNodesList.shape[0], VN.numpy_to_vtkIdTypeArray(connectivity))

            uGrid = vtk.vtkUnstructuredGrid()
            uGrid.SetPoints(pointsVTK)
            uGrid.SetCells(cell_types, cellsVTK)

            cell_data = uGrid.GetCellData()  # This holds cell data
            point_data = uGrid.GetPointData()  # This holds point data

            # add solutions

            # column numbers for Vel_X and Vel_Y for vector assemble
            nColVel_X = -1
            nColVel_Y = -1

            # First output all solution variables as scalars
            print('The following solution variables are processed and saved to VTK file: \n')
            for k in range(len(resultVarNames)):
                print('     %s\n' % resultVarNames[k])

                #if it is a veloctiy component, only record its location
                #not output to VTK by itself, will be assembed to vector.
                if resultVarNames[k].find('Vel_X') != -1:
                    nColVel_X = k
                    continue
                elif resultVarNames[k].find('Vel_Y') != -1:
                    nColVel_Y = k
                    continue

                if bNodal:
                    temp_point_data_array = VN.numpy_to_vtk(resultData[:,k])
                    temp_point_data_array.SetName(resultVarNames[k])
                    point_data.AddArray(temp_point_data_array)
                else:
                    temp_cell_data_array = VN.numpy_to_vtk(resultData[:,k])
                    temp_cell_data_array.SetName(resultVarNames[k])
                    cell_data.AddArray(temp_cell_data_array)

            #add veloctiy by combining components
            currentTimePointV = np.zeros((self.srhgeom_obj.numOfNodes, 3)) if bNodal else \
                                np.zeros((self.srhgeom_obj.numOfElements,3))
            currentTimePointV[:, 0] = resultData[:,nColVel_X]
            currentTimePointV[:, 1] = resultData[:,nColVel_Y]
            currentTimePointV[:, 2] = 0.0

            if bNodal:
                temp_point_data_array = VN.numpy_to_vtk(currentTimePointV)
                temp_point_data_array.SetName(velocityVarName)
                point_data.AddArray(temp_point_data_array)
            else:
                temp_cell_data_array = VN.numpy_to_vtk(currentTimePointV)
                temp_cell_data_array.SetName(velocityVarName)
                cell_data.AddArray(temp_cell_data_array)

            #add nodal bed elevation
            currentBedElev = np.zeros(self.srhgeom_obj.numOfNodes) if bNodal else np.zeros(self.srhgeom_obj.numOfElements)

            currentBedElev = self.srhgeom_obj.nodeCoordinates[:,2]

            temp_point_data_array = VN.numpy_to_vtk(currentBedElev)
            temp_point_data_array.SetName(bedElevationVarName)
            point_data.AddArray(temp_point_data_array)

            # write to vtk file
            unstr_writer = vtk.vtkUnstructuredGridWriter()
            unstr_writer.SetFileName(vtkFileName)
            unstr_writer.SetInputData(uGrid)
            unstr_writer.Write()


    def outputVTK(self, vtkFileName, resultVarNames, resultData, bNodal):
        """ Output result to VTK file

        The supplied resultVarNames and resultData should be compatible with the mesh. If resultVarNames is empty,
        it only outputs the mesh with no data.

        Parameters
        ----------
        vtkFileName : str
            name for the output vtk file
        resultVarNames : str
            result variable names
        resultData : `numpy.ndarray`
            2D array containing result data
        bNodal : bool
            whether the data is nodal (True) or at cell center (False)


        Returns
        -------

        """

        print("Output to VTK ...")

        try:
            fid = open(vtkFileName, 'w')
        except IOError:
            print('vtk file open error')
            sys.exit()

        #get the units of the results
        units = self.srhhydro_obj.srhhydro_content['OutputFormat'][1]
        print("SRH-2D result units:", units)

        fid.write('# vtk DataFile Version 3.0\n')
        fid.write('Results from SRH-2D Modeling Run\n')
        fid.write('ASCII\n')
        fid.write('DATASET UNSTRUCTURED_GRID\n')
        fid.write('\n')

        # output points
        fid.write('POINTS %d double\n' % self.srhgeom_obj.nodeCoordinates.shape[0])

        point_id = 0  # point ID counter
        for k in range(self.srhgeom_obj.nodeCoordinates.shape[0]):
            point_id += 1
            fid.write(" ".join(map(str, self.srhgeom_obj.nodeCoordinates[k])))
            fid.write("\n")

        # output elements
        fid.write('CELLS %d %d \n' % (self.srhgeom_obj.elementNodesList.shape[0],
                                      self.srhgeom_obj.elementNodesList.shape[0] + np.sum(self.srhgeom_obj.elementNodesCount)))

        cell_id = 0  # cell ID counter
        for k in range(self.srhgeom_obj.elementNodesList.shape[0]):
            cell_id += 1
            fid.write('%d ' % self.srhgeom_obj.elementNodesCount[k])
            fid.write(" ".join(map(str, self.srhgeom_obj.elementNodesList[k][:self.srhgeom_obj.elementNodesCount[k]] - 1)))
            fid.write("\n")

        # output element types
        fid.write('CELL_TYPES %d \n' % self.srhgeom_obj.elementNodesList.shape[0])

        for k in range(self.srhgeom_obj.elementNodesList.shape[0]):
            fid.write('%d ' % self.srhgeom_obj.vtkCellTypeCode[k])
            if (((k + 1) % 20) == 0):
                fid.write('\n')

        fid.write('\n')

        # How many solution data entries to output depending on whehter it is nodal or not
        entryCount = -1
        if bNodal:
            entryCount = self.srhgeom_obj.nodeCoordinates.shape[0]
        else:
            entryCount = self.srhgeom_obj.elementNodesList.shape[0]

        # output solution variables: only if there is solution variable
        if len(resultVarNames) != 0:
            if not bNodal:
                print('Solution variables are at cell centers. \n')
                fid.write('CELL_DATA %d\n' % self.srhgeom_obj.elementNodesList.shape[0])
            else:
                print('Solution variables are at vertices. \n')
                fid.write('POINT_DATA %d\n' % self.srhgeom_obj.nodeCoordinates.shape[0])

            # column numbers for Vel_X and Vel_Y for vector assemble
            nColVel_X = -1
            nColVel_Y = -1

            # First output all solution variables as scalars
            #print('The following solution variables are processed: \n')
            for k in range(len(resultVarNames)):
                #print('     %s\n' % resultVarNames[k])

                if resultVarNames[k].find('Vel_X') != -1:
                    nColVel_X = k
                elif resultVarNames[k].find('Vel_Y') != -1:
                    nColVel_Y = k

                fid.write('SCALARS %s double 1 \n' % resultVarNames[k])
                fid.write('LOOKUP_TABLE default\n')

                for entryI in range(entryCount):
                    fid.write('%f ' % resultData[entryI][k])

                    if (((entryI + 1) % 20) == 0):
                        fid.write('\n')

                fid.write('\n \n')

            # Then output Vel_X and Vel_Y as velocity vector (Vel_Z = 0.0)
            # print('nColVel_X, nColVel_Y = %d, %d' % (nColVel_X, nColVel_Y))
            if (nColVel_X != -1) and (nColVel_Y != -1):
                if units == "SI":
                    fid.write('VECTORS Velocity_m_p_s double \n')
                else:
                    fid.write('VECTORS Velocity_ft_p_s double \n')

                for entryI in range(entryCount):
                    fid.write('%f %f 0.0   ' % (resultData[entryI][nColVel_X],
                                                resultData[entryI][nColVel_Y]))

                    if (((entryI + 1) % 20) == 0):
                        fid.write('\n')

                fid.write('\n \n')

        fid.close()

    def output_2d_mesh_to_vtk(self, flatMeshVTKFileName, bFlat=False, dir=''):
        """ output the flat mesh to vtk

        Parameters
        ----------
        flatMeshVTKFileName : str
            file name for the flat mesh
        bFlat : bool
            whether to make the 2d mesh flat (node's z coordinate -> 0)
        dir : str, optional
            directory to write to

        Returns
        -------

        """
        #just call srhgeom_obj's function
        self.srhgeom_obj.output_2d_mesh_to_vtk(flatMeshVTKFileName, bFlat, dir)

    def readSRHFile(self, srhFileName):
        """ Read SRH-2D result file in SRHC (cell center) or SRH (point) format.

        Note: SRH-2D outputs an extra "," to each line. As a result, Numpy's
        genfromtext(...) function adds a column of "nan" to the end. Need to take care of this.

        Parameters
        ----------
        srhFileName : str
            file name for the SRH result file

        Returns
        -------


        """

        print("Reading the SRH/SRHC result file ...")

        data = np.genfromtxt(srhFileName, delimiter=',', names=True)

        return data.dtype.names[:-1], data

    def readTECFile(self, tecFileName):
        """Read SRH-2D results in Tecplot format

        Parameters
        ----------
        tecFileName : str
            Tecplot file name

        Returns
        -------

        """

        readerTEC = vtk.vtkTecplotReader()
        readerTEC.SetFileName(tecFileName)
        # 'update' the reader i.e. read the Tecplot data file
        readerTEC.Update()

        polydata = readerTEC.GetOutput()

        # print(polydata)  #this is a multibock dataset
        # print(polydata.GetBlock(0))  #we only need the first block

        # If there are no points in 'vtkPolyData' something went wrong
        if polydata.GetBlock(0).GetNumberOfPoints() == 0:
            raise ValueError(
                "No point data could be loaded from '" + tecFileName)
            return None

        #return of the first block (there should be only one block)
        return polydata.GetBlock(0)


    def cell_center_to_vertex(self, cell_data, vertex_data):
        """Interpolate result from cell center to vertex

        Not implemented yet.

        Returns
        -------

        """
        pass

    def vertex_to_cell_center(self, vertex_data, cell_data):
        """Interpolate result from vertex to cell center

        Not implemented yet.

        Returns
        -------

        """
        pass
