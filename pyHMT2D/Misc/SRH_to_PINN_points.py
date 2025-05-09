"""
This set of tools accomplish this task: convert a SRH-2D case's 2D mesh to points needed for using NN to solve PDEs.
    - "points" include both equation points and boundary points
    - Save format: JSON
"""

import numpy as np
import sys
import meshio
import math
import json
from .tools import generate_random01_exclude_boundaries_with_center, point_on_triangle, point_on_line
import pyHMT2D

def srh_to_pinn_points(srhcontrol_file, refinement=1):
    """Convert SRH-2D case's 2D mesh to points for PINN training. The saved "mesh_points.json" file can be used for PINN training. The json file contains "equation_points" and "boundary_points".

    The saved points only have spatial information. The time information is not included.

    Parameters
    ----------
    srhcontrol_file : str
        Name of the SRH-2D case's control file, e.g. "case.srhhydro" or "case_SIF.dat". It can only be a srhhydro or SIF file.
    refinement : int
        Number of points to generate per cell (>= 1)

    Returns
    -------
    None
        Saves points to mesh_points.json
        Saves domain and boundary meshes to domain_{filename}.xdmf and boundaries_{filename}.xdmf
    """

    #check the file extension
    if not srhcontrol_file.endswith(".srhhydro") and not srhcontrol_file.endswith("_SIF.dat"):
        raise ValueError("The input file must be a srhhydro or SIF.dat file")

    #create the SRH_2D_Data object
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(srhcontrol_file)

    #get the SRH-2D mesh
    my_srh_2d_mesh = my_srh_2d_data.srhgeom_obj

    #process the internal points: Loop through all cells and generate points
    equation_points_dict = {}
    equation_points_list = []  # List to store points for VTK output
    pointID = 0

    print(f"Total number of cells in the mesh: {my_srh_2d_mesh.numOfElements}")

    for iCell in range(my_srh_2d_mesh.numOfElements):
        #get the cell's vertices
        cell_vertices = my_srh_2d_mesh.elementNodesList[iCell, :my_srh_2d_mesh.elementNodesCount[iCell]]
        num_vertices = my_srh_2d_mesh.elementNodesCount[iCell]
        
        # Get vertex coordinates
        vertex_coords = [my_srh_2d_mesh.nodeCoordinates[vid-1] for vid in cell_vertices]  # -1 because SRH-2D is 1-based
        
        if num_vertices == 3:  # Triangle
            p0, p1, p2 = vertex_coords
            
            # Generate sampling points using barycentric coordinates
            st_all = generate_random01_exclude_boundaries_with_center(
                centers=[1.0/3.0, 2.0/3.0],
                size=refinement
            )

            for st in st_all:
                point = point_on_triangle(p0, p1, p2, st[0], st[1])
                equation_points_list.append([point[0], point[1], point[2]])
                
                equation_points_dict[str(pointID)] = {
                    "x": point[0],
                    "y": point[1],
                    "z": point[2],
                    "spatial_dimensionality": 2
                }
                pointID += 1
                
        elif num_vertices == 4:  # Quad
            p0, p1, p2, p3 = vertex_coords
            
            # Generate sampling points using bilinear interpolation
            st_all = generate_random01_exclude_boundaries_with_center(
                centers=[0.5, 0.5],
                size=refinement
            )

            for st in st_all:
                s, t = st
                # Bilinear interpolation
                bottom_edge = (1-s)*p0 + s*p1
                top_edge = (1-s)*p3 + s*p2
                point = (1-t)*bottom_edge + t*top_edge
                
                equation_points_list.append([point[0], point[1], point[2]])
                
                equation_points_dict[str(pointID)] = {
                    "x": point[0],
                    "y": point[1],
                    "z": point[2],
                    "spatial_dimensionality": 2
                }
                pointID += 1
        else:
            raise ValueError(f"Unsupported cell type: {num_vertices}")

    # Write equation points for visualization
    equation_points = np.array(equation_points_list)
    meshio.write_points_cells("equation_points.vtk", equation_points, {}, {}, {}, binary=True)

    # Process boundary points
    boundary_points_dict = {}
    nodes = my_srh_2d_mesh.nodeCoordinates

    # Calculate total number of boundary points
    total_boundary_points = 0
    for bc_id, edge_ids in my_srh_2d_mesh.boundaryEdges.items():
        total_boundary_points += len(edge_ids) * refinement

    # Initialize arrays for boundary points
    bc_points = np.zeros((total_boundary_points, 3))
    bc_normal_vectors = np.zeros((total_boundary_points, 3))
    bc_IDs = np.zeros(total_boundary_points, dtype=np.int64)
    bc_represented_lengths = np.zeros(total_boundary_points, dtype=np.float64)

    pointID = 0
    #Loop over all boundaries
    for bc_id, edge_ids in my_srh_2d_mesh.boundaryEdges.items():
        for edge_id in edge_ids:
            # Get the nodes for this edge 
            edge_nodes = my_srh_2d_mesh.edges_r[abs(edge_id)]  #edge_id might be negative when the direction of the edge is the opposite of boundary
            node1, node2 = edge_nodes[0], edge_nodes[1]
                        
            if edge_id > 0:           #ensure the correct order of the nodes
                p1 = nodes[node1-1]  
                p2 = nodes[node2-1]  
            else:
                p1 = nodes[node2-1]  
                p2 = nodes[node1-1]  
            
            line_length = np.sqrt(np.sum((p2 - p1)**2))
            represented_length = line_length / refinement

            # Generate sampling points
            if refinement == 1:
                s_all = np.array([0.5])
            elif refinement == 2:
                s_all = np.array([0.333, 0.666])
            else:
                s_all = np.linspace(0.05, 0.95, refinement)

            for s in s_all:
                point = point_on_line(p1, p2, s)
                bc_IDs[pointID] = bc_id  # Boundary ID
                bc_points[pointID, :3] = point

                # Calculate and normalize normal vector
                normal = np.array([p2[1] - p1[1], -(p2[0] - p1[0]), 0.0])
                normal /= np.linalg.norm(normal[:2])
                bc_normal_vectors[pointID] = normal

                bc_represented_lengths[pointID] = represented_length
                pointID += 1

    # Write boundary points for visualization
    point_data = {
        'bc_ID': bc_IDs,
        'normal_vector': bc_normal_vectors,
        'represented_length': bc_represented_lengths
    }
    meshio.write_points_cells("boundary_points.vtk", bc_points, {}, point_data, {}, binary=True)

    # Organize points by boundary ID
    unique_BC_IDs = np.unique(bc_IDs)
    for bc_ID in unique_BC_IDs:
        mask = bc_IDs == bc_ID
        current_boundary_dict = {}
        
        for i, point in enumerate(bc_points[mask]):
            current_boundary_dict[str(i)] = {
                "x": point[0],
                "y": point[1],
                "z": 0.0,
                "normal_x": bc_normal_vectors[mask][i][0],
                "normal_y": bc_normal_vectors[mask][i][1],
                "normal_z": bc_normal_vectors[mask][i][2],
                "spatial_dimensionality": 2,
                "represented_length": bc_represented_lengths[mask][i]
            }

        boundary_points_dict[f"boundary_{bc_ID}"] = current_boundary_dict

    # Assemble all points
    all_points = {
        "training_points": {
            "equation_points": equation_points_dict,
            "boundary_points": boundary_points_dict
        }
    }

    # Save to JSON
    with open("mesh_points.json", 'w') as f:
        json.dump(all_points, f, indent=4, sort_keys=False)


    # Process the boundary conditions: create a json file for boundary conditions dictionary from srhcontrol_file depending whether it is a srhhydro or SIF file
    boundary_conditions_dict = {}

    if srhcontrol_file.endswith(".srhhydro"):
        # Read the srhhydro file
        with open(srhcontrol_file, 'r') as f:
            lines = f.readlines()



            

    
    
