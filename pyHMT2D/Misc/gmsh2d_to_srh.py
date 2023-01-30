"""
A tool to convert Gmsh's 2D MSH file to srhgeom format (to be used by SRH-2D). Optionally, it can interpolate terrain
to the 2D mesh. Gmsh's 2D mesh usually do not have elevation information.
"""

import numpy as np
import sys
import meshio


def gmsh2d_to_srh(gmsh2d_fileName, srh_caseName, shift_x=0.0, shift_y=0.0, units="Meters"):
    """Convert Gmsh 2D mesh into SRH-2D format

    It generates two files: srhgeom for mesh and srhmat for Manning's n

    The srhhydro file has to be generated separately.

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

    Returns
    -------

    """

    print("Converting Gmsh's MSH mesh to SRH-2D format ...")

    #read in the Gmsh MSH with meshio
    mesh = meshio.read(gmsh2d_fileName)

    #output the mesh as vtk for checking
    mesh.write("check_mesh.vtk")

    #build the nodeStrings (boundaries)
    nodeStrings = build_nodeStrings(mesh)

    #write srhgeom file
    write_srhgeom(srh_caseName, mesh, nodeStrings, shift_x, shift_y, units)

    #build Manning's n zone information
    nManningNZones, ManningNZoneNames, cellsInManningZones = build_ManningNZones(mesh)

    #write srhmat file
    write_srhmat(srh_caseName, nManningNZones, ManningNZoneNames, cellsInManningZones)

    print("Finished converting Gmsh's MSH mesh to SRH-2D format.")


def build_nodeStrings(mesh):
    """Build the nodeStrings (boundary lines)

    Parameters
    ----------
    mesh : meshio's mesh
        meshio's mesh object

    Returns
    -------

    """

    #get all the lines's physcial group number
    lines_physical_group_IDs = mesh.cell_data_dict['gmsh:physical']['line']

    #how many different line boundaries
    unique_line_physical_group_IDs = np.unique(lines_physical_group_IDs)

    #get all lines' nodes (assume there is only one line list. Sure?)
    lines_nodes = mesh.cells[0].data

    #build node strings for each line boundary
    nodeStrings = {}

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

                if (lines_nodes[line_ID][1]+1) not in current_node_list:
                    current_node_list.append(lines_nodes[line_ID][1]+1)  #+1 because meshio is 0-based; gmsh is 1-based

        #add the nodeString to the nodeStrings dictionary
        nodeStrings[nodeStringID] = current_node_list

    #print(type(nodeStrings))

    return nodeStrings

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
