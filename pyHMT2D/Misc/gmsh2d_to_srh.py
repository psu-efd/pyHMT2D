"""
A tool to convert Gmsh's 2D MSH file to srhgeom format (to be used by SRH-2D). Optionally, it can interpolate terrain
to the 2D mesh. Gmsh's 2D mesh usually do not have elevation information.
"""

import numpy as np
import sys
import meshio
import math


def gmsh2d_to_srh(gmsh2d_fileName, srh_caseName, shift_x=0.0, shift_y=0.0, units="Meters",
                  bAddMonitoringLines=False, monitoringLines=[], monitoringlineTol=0.001
                  ):
    """Convert Gmsh 2D mesh into SRH-2D format with the option to add monitoring lines.

    It generates two files: srhgeom for mesh and srhmat for Manning's n

    The srhhydro file has to be generated separately.

    For monitoring lines: SRH-2D treats the monitoring lines similarly as boundaries. Both are nodeStrings in the srhgeom
    file. However, Gmsh2d does not directly support internal boundary such as the monitoring lines. Thus, for any
    monitoring line, we need to add it to the nodeStrings list here.

    Parameters
    ----------
    gmsh2d_fileName : str
        file name of the Gmsh MSH file
    srh_caseName : str
        case name for SRH-2D.
    shift_x: float
        shift of coordinates in x direction
    shift_y: float
        shift of coordinates in y direction
    units : str, default Meters
        length units of Gmsh file
    bAddMonitoringLines: whether add monitoring lines
    monitoringLines: list of monitoring lines. Currently only straightlines are supported.
    monitoringlineTol: a float tolerance to check whether a node is close to the monitoring line

    Returns
    -------

    """

    print("Converting Gmsh's MSH mesh to SRH-2D format ...")

    #read in the Gmsh MSH with meshio
    mesh = meshio.read(gmsh2d_fileName)

    #output the mesh as vtk for checking
    mesh.write("check_mesh.vtk", file_format="vtk42")

    #build the nodeStrings (boundaries)
    nodeStrings = build_nodeStrings(mesh,bAddMonitoringLines,monitoringLines,monitoringlineTol)

    #write srhgeom file
    write_srhgeom(srh_caseName, mesh, nodeStrings, shift_x, shift_y, units)

    #build Manning's n zone information
    nManningNZones, ManningNZoneNames, cellsInManningZones = build_ManningNZones(mesh)

    #write srhmat file
    write_srhmat(srh_caseName, nManningNZones, ManningNZoneNames, cellsInManningZones)

    print("Finished converting Gmsh's MSH mesh to SRH-2D format.")


def build_nodeStrings(mesh,bAddMonitoringLines,monitoringLines,monitoringlineTol):
    """Build the nodeStrings (boundary lines)

    Parameters
    ----------
    mesh : meshio's mesh
        meshio's mesh object
    bAddMonitoringLines: whether there is any monitoring line
    monitoringLines: a list for the definitions of monitoring lines (straightlines only for now)
    monitoringlineTol: a float tolerance to check whether a node is close to monitoring line.

    Returns
    -------

    """

    #get all the lines's physical group number
    lines_physical_group_IDs = mesh.cell_data_dict['gmsh:physical']['line']

    #how many different line boundaries
    unique_line_physical_group_IDs = np.unique(lines_physical_group_IDs)

    #get all lines' nodes (assume there is only one line list. Sure?)
    lines_nodes = mesh.cells[0].data

    #build node strings for each line boundary
    nodeStrings = {}

    #loop through each line group
    for nodeStringID in unique_line_physical_group_IDs:
        #print("nodeStringID =", nodeStringID)

        current_node_list = []

        #loop through each line in lines_physical_group_IDs:
        for line_group_ID, line_ID in zip(lines_physical_group_IDs, range(len(lines_physical_group_IDs))):
            #print("line_group_ID, line_ID =", line_group_ID, line_ID)

            if nodeStringID == line_group_ID:
                # add the nodes of line to current_node_list (if they are not in already)
                if (lines_nodes[line_ID][0]+1) not in current_node_list:
                    current_node_list.append(lines_nodes[line_ID][0]+1)  #+1 because meshio is 0-based; gmsh is 1-based

                #For a closed loop, one node shows up in the list twice; thus it can repeat. Need to consider this scenario, i.e., 
                #if the node is already in the list, but not the last one, then we need to add the last one.
                if (lines_nodes[line_ID][0]+1) in current_node_list and current_node_list[-1] != (lines_nodes[line_ID][0]+1):
                    current_node_list.append(lines_nodes[line_ID][0]+1)  #+1 because meshio is 0-based; gmsh is 1-based

                if (lines_nodes[line_ID][1]+1) not in current_node_list:
                    current_node_list.append(lines_nodes[line_ID][1]+1)  #+1 because meshio is 0-based; gmsh is 1-based

                if (lines_nodes[line_ID][1]+1) in current_node_list and current_node_list[-1] != (lines_nodes[line_ID][1]+1):
                    current_node_list.append(lines_nodes[line_ID][1]+1)  #+1 because meshio is 0-based; gmsh is 1-based

        #add the nodeString to the nodeStrings dictionary
        nodeStrings[nodeStringID] = current_node_list

    #deal with monitoring lines
    if bAddMonitoringLines:

        #start ID for the monitoring lines
        nodeStringID = max(unique_line_physical_group_IDs) + 1

        #get the edge lists. For now, monitoring lines are only defined with internal edges.
        boundary_edges, internal_edges = get_msh_edges(mesh)

        for monitoringLine in monitoringLines:

            current_node_list = []

            xML_start = monitoringLine[0, 0]
            yML_start = monitoringLine[0, 1]
            xML_end = monitoringLine[1, 0]
            yML_end = monitoringLine[1, 1]

            #loop over all internal edges and check whether they are close to the monitoring lines
            for edge in internal_edges:
                nodeIDStart = edge[0]
                nodeIDEnd = edge[1]

                coordNodeStart = mesh.points[nodeIDStart]
                coordNodeEnd = mesh.points[nodeIDEnd]

                distanceNodeStart = point_to_segment_distance(coordNodeStart[0], coordNodeStart[1],
                                                              xML_start, yML_start,
                                                              xML_end, yML_end)

                distanceNodeEnd = point_to_segment_distance(coordNodeEnd[0], coordNodeEnd[1],
                                                              xML_start, yML_start,
                                                              xML_end, yML_end)

                #if (distanceNodeStart+distanceNodeEnd)/2 < monitoringlineTol:
                if distanceNodeStart < monitoringlineTol and distanceNodeEnd < monitoringlineTol:
                    if (nodeIDStart + 1) not in current_node_list:
                        current_node_list.append(nodeIDStart + 1)  # +1 because meshio is 0-based; gmsh is 1-based

                    if (nodeIDEnd + 1) not in current_node_list:
                        current_node_list.append(nodeIDEnd + 1)  # +1 because meshio is 0-based; gmsh is 1-based

            # add the nodeString to the nodeStrings dictionary
            nodeStrings[nodeStringID] = current_node_list

            nodeStringID += 1

    #print(type(nodeStrings))

    return nodeStrings

def point_to_segment_distance(x0, y0, x1, y1, x2, y2):
    # Calculate the squared length of the segment
    segment_len_squared = (x2 - x1) ** 2 + (y2 - y1) ** 2

    # Handle the case where the segment length is zero (start and end points are the same)
    if segment_len_squared == 0:
        return math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)

    # Calculate the projection of (x0, y0) onto the line (x1, y1) to (x2, y2)
    t = ((x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1)) / segment_len_squared

    # Clamp t to the range [0, 1] to find the closest point on the segment
    t = max(0, min(1, t))

    # Find the closest point on the segment to (x0, y0)
    closest_x = x1 + t * (x2 - x1)
    closest_y = y1 + t * (y2 - y1)

    # Calculate the distance from (x0, y0) to the closest point on the segment
    distance = math.sqrt((x0 - closest_x) ** 2 + (y0 - closest_y) ** 2)

    return distance

def get_msh_edges(mesh):
    """
    Get the boundary and internal lists of edges.

    :param mesh:
    :return:
    """

    # Initialize edge counter dictionary
    edge_count = {}

    # Process 2D elements (triangles and quadrilaterals)
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            for triangle in cell_block.data:
                # Extract triangle edges
                edges = [
                    tuple(sorted((triangle[0], triangle[1]))),
                    tuple(sorted((triangle[1], triangle[2]))),
                    tuple(sorted((triangle[2], triangle[0])))
                ]
                # Count each edge
                for edge in edges:
                    edge_count[edge] = edge_count.get(edge, 0) + 1

        elif cell_block.type == "quad":
            for quad in cell_block.data:
                # Extract quadrilateral edges
                edges = [
                    tuple(sorted((quad[0], quad[1]))),
                    tuple(sorted((quad[1], quad[2]))),
                    tuple(sorted((quad[2], quad[3]))),
                    tuple(sorted((quad[3], quad[0])))
                ]
                # Count each edge
                for edge in edges:
                    edge_count[edge] = edge_count.get(edge, 0) + 1

    # Separate boundary and internal edges based on their counts
    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
    internal_edges = [edge for edge, count in edge_count.items() if count == 2]

    return boundary_edges, internal_edges

def build_ManningNZones(mesh):
    """Build Manning's n zones

    Currently, only one zone is supported. All elements will be in a single zone. In future, it may be
    possible to use Gmsh with multiple physical zones.

    Parameters
    ----------
    mesh : meshio's mesh
        meshio's mesh object

    Returns
    -------

    """


    #Only one Manning's n zone is supported currently.
    nManningNZones = 1
    ManningNZoneNames = ['channel']

    cellsInManningZones = {}

    # all cells
    cellCounter = 0

    for cells_blockI in range(len(mesh.cells)):
        # don't need the lines
        if mesh.cells[cells_blockI].type == 'line':
            continue

        # triangles
        elif mesh.cells[cells_blockI].type == 'triangle':
            # loop over all triangles
            for triI in range(mesh.cells[cells_blockI].data.shape[0]):
                cellCounter += 1

        # quads
        elif mesh.cells[cells_blockI].type == 'quad':
            # loop over all triangles
            for quadI in range(mesh.cells[cells_blockI].data.shape[0]):
                cellCounter += 1
        else:
            raise Exception("Gmesh element type is not supported.")

    cellID_list = []

    for cellI in range(cellCounter):
        cellID_list.append(cellI+1)  #+1 because SRH-2D is 1-based

    cellsInManningZones[0] = cellID_list

    return nManningNZones, ManningNZoneNames, cellsInManningZones

def orientation_2D(xA, yA, xB, yB, xC, yC):
    """
    Calculate the orientation of the order of three points A, B, and C in 2D.

    Ref: https://en.wikipedia.org/wiki/Curve_orientation

    :param xA:
    :param yA:
    :param xB:
    :param yB:
    :param xC:
    :param yC:
    :return: < 0: clockwise, >0: counter clockwise
    """

    #construct the square matrix
    matrix = np.array([ [1, xA, yA],[1, xB, yB],[1, xC, yC] ])

    return np.linalg.det(matrix)

def write_srhgeom(srhmatFileName, mesh, nodeStrings, shift_x=0.0, shift_y=0.0, units="Meters"):
    """


    Parameters
    ----------
    srhmatFileName
    mesh
    nodeStrings
    units

    Returns
    -------

    """

    fname = srhmatFileName + ".srhgeom"

    try:
        fid = open(fname, 'w')
    except IOError:
        print('.srhgeom error')
        sys.exit()

    fid.write('SRHGEOM 30\n')
    fid.write('Name \"Converted from Gmsh 2D Mesh \"\n')

    fid.write('\n')

    fid.write('GridUnit \"%s\" \n' % units)

    # all cells

    cellI = 0

    for cells_blockI in range(len(mesh.cells)):
        #don't need the lines
        if mesh.cells[cells_blockI].type == 'line':
            continue

        #triangles
        elif mesh.cells[cells_blockI].type == 'triangle':
            #loop over all triangles
            for triI in range(mesh.cells[cells_blockI].data.shape[0]):
                cellI += 1

                #determine the orientation of the triangle
                xA = mesh.points[mesh.cells[cells_blockI].data[triI][0], 0]
                yA = mesh.points[mesh.cells[cells_blockI].data[triI][0], 1]

                xB = mesh.points[mesh.cells[cells_blockI].data[triI][1], 0]
                yB = mesh.points[mesh.cells[cells_blockI].data[triI][1], 1]

                xC = mesh.points[mesh.cells[cells_blockI].data[triI][2], 0]
                yC = mesh.points[mesh.cells[cells_blockI].data[triI][2], 1]

                sign = orientation_2D(xA, yA, xB, yB, xC, yC)

                fid.write("Elem ")

                if sign > 0: #nodes in counterclockwise direction
                    fid.write("%d %d %d %d \n" % (cellI,
                                              mesh.cells[cells_blockI].data[triI][0]+1,
                                              mesh.cells[cells_blockI].data[triI][1]+1,
                                              mesh.cells[cells_blockI].data[triI][2]+1
                                              ))
                else:        #nodes in clockwise direction; need to revise the order
                    fid.write("%d %d %d %d \n" % (cellI,
                                                  mesh.cells[cells_blockI].data[triI][0] + 1,
                                                  mesh.cells[cells_blockI].data[triI][2] + 1,
                                                  mesh.cells[cells_blockI].data[triI][1] + 1
                                                  ))


        #quads
        elif mesh.cells[cells_blockI].type == 'quad':
            # loop over all triangles
            for quadI in range(mesh.cells[cells_blockI].data.shape[0]):
                cellI += 1

                # determine the orientation of the quad (only need to use the first three nodes)
                xA = mesh.points[mesh.cells[cells_blockI].data[quadI][0], 0]
                yA = mesh.points[mesh.cells[cells_blockI].data[quadI][0], 1]

                xB = mesh.points[mesh.cells[cells_blockI].data[quadI][1], 0]
                yB = mesh.points[mesh.cells[cells_blockI].data[quadI][1], 1]

                xC = mesh.points[mesh.cells[cells_blockI].data[quadI][2], 0]
                yC = mesh.points[mesh.cells[cells_blockI].data[quadI][2], 1]

                sign = orientation_2D(xA, yA, xB, yB, xC, yC)

                fid.write("Elem ")

                if sign > 0:  # nodes in counterclockwise direction
                    fid.write("%d %d %d %d %d\n" % (cellI,
                                              mesh.cells[cells_blockI].data[quadI][0]+1,
                                              mesh.cells[cells_blockI].data[quadI][1]+1,
                                              mesh.cells[cells_blockI].data[quadI][2]+1,
                                              mesh.cells[cells_blockI].data[quadI][3]+1
                                              ))
                else:  # nodes in clockwise direction; need to revise the order
                    fid.write("%d %d %d %d %d\n" % (cellI,
                                                    mesh.cells[cells_blockI].data[quadI][0] + 1,
                                                    mesh.cells[cells_blockI].data[quadI][3] + 1,
                                                    mesh.cells[cells_blockI].data[quadI][2] + 1,
                                                    mesh.cells[cells_blockI].data[quadI][1] + 1
                                                    ))
        else:
            raise Exception("Gmesh element type is not supported.")


    # all points
    for pointI in range(mesh.points.shape[0]):
        fid.write("Node %d " % (pointI + 1))  # pointI+1 because SRH-2D is 1-based
        curr_point_coordinates = [mesh.points[pointI, 0] + shift_x,
                                  mesh.points[pointI, 1] + shift_y,
                                  mesh.points[pointI, 2]]

        fid.write(" ".join(map(str, curr_point_coordinates)))
        fid.write("\n")

    # NodeString
    boundary_id = 0  # boundary ID counter (not the ID used in Gmsh)



    # loop through all boundaries
    for nodeStringID, node_list in nodeStrings.items():
        boundary_id += 1

        fid.write("NodeString %d " % boundary_id)

        # line break counter (start a new line every 10 nodes)
        line_break_counter = 0

        # loop over each node ID in the current NodeString
        for nodeID in node_list:

            line_break_counter += 1

            fid.write(" %d" % (nodeID))

            # 10 numbers per line
            if ((line_break_counter % 10) == 0):

                fid.write("\n")

                line_break_counter = 0

        fid.write("\n")

    fid.close()


def write_srhmat(srhmatFileName, nManningNZones, ManningNZoneNames, cellsInManningZones):
    """Export the SRHMAT file

            Parameters
            ----------
            srhmatFileName : str
                name of the srhmat file to write to
            mesh:
                meshio object
            nManningNZones: int
                number of Manning's n zones
            ManningNZoneNames: list
                list of Manning's n zone names
            cellsInManningZones: dict
                cells in each Manning's n zones

            Returns
            -------

            """

    fname = srhmatFileName + '.srhmat'

    try:
        fid = open(fname, 'w')
    except IOError:
        print('.srhmat error')
        sys.exit()

    fid.write('SRHMAT 30\n')
    fid.write('NMaterials %d\n' % (
            nManningNZones + 1))  # +1 is because SRH-2D also counts the default Manning's n in srhhydro file.

    # output MatName
    for matID in range(nManningNZones):
        fid.write('MatName %d \"%s\" \n' % (
            matID + 1, ManningNZoneNames[matID]))  # +1 because SRH-2D is 1-based

    # output cells in different material categories
    for matID in range(nManningNZones):
        if not cellsInManningZones[matID]:  # this Manning's n zone has no cells
            continue

        fid.write('Material %d ' % (matID + 1))

        # loop over all cells in current Manning's n zone
        for cellI in range(len(cellsInManningZones[matID])):
            fid.write(" %d" % (cellsInManningZones[matID][cellI]))

            # 10 numbers per line or this is the last cell, start a new line
            if (((cellI + 1) % 10) == 0) or (cellI == (len(cellsInManningZones[matID]) - 1)):
                fid.write("\n")

    fid.close()
