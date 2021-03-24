"""
Modified from the following (accessed 3/21/2021; no copyright specified, public domain).
http://jrper.github.io/2016/08/11/GMSHToVTKTool.html

"""

import vtk
from vtk.util import numpy_support
import numpy
import struct

class gmsh2DMesh(object):
    """GMSH 2D Mesh class for 2D Hydraulic model meshes.

    gmsh2DMesh is desgined to:
    1. store the same 2D Hydraulic model mesh information (nodes, cells, boundaries, etc) using GMSH's MSH format.
       https://gmsh.info/doc/texinfo/gmsh.html#File-formats
    2. extrude 2D mesh to 3D (either with only one layer or multiple layers), to be used in 3D models, for example OpenFOAM.
    3. save both 2D and extruded 3D mesh to GMSH mesh format.

    gmsh2DMesh can be constructed from SRH-2D's SRH_2D_SRHGeom object (nodes, elements, and boundary NodeStrings).

    Notes:
        1. GMSH is 1-based (for IDs)
        2. "Dimension" in GMSH: 0 - vertex/point/node, 1 - line/edge/curve, 2 - surface/face, 3 - volume/cell
        3. Most important parts of MSH file for our purpose:
            a. $MeshFormat  ... $EndMeshFormat
            b. $PhysicalNames ... $EndPhysicalNames (for the names of boundaries and cell zones if any)
            c. $Entities ... $EndEntities (for the definition of boundaries and cell zones). "Entities" in GMSH are
               elementary model entities, which include points, curves, surfaces, and volumes. Here, points are defined
               in GMSH's geo file, not the nodes for mesh.

               In OpenFOAM's gmshToFoam tool, only surface and volumes are read (points and curves are ignored).
            d. $Nodes ... $EndNodes. All nodes in the mesh. Again, the nodes are different from points in GMSH.
            e. $Elements ... $EndElements. All elements in the mesh.

    Attributes:
    -----

    Methods:
    -----


    """

    #class variables for a supported subset of GMSH elementType (for pyHMT2D's purpose)
    #Currently, only triangles and quadrilaterals are allowed in 2D mesh (correspondingly only
    #hexahydron and prim are allowed in the extruded 3D mesh. This is ok in current version of SMS, but not
    #in HEC-RAS.
    MSHLINE = 1  #2-node line

    MSHTRI = 2   #3-node triangle
    MSHQUAD = 3  #4-node quadrilateral

    MSHHEX = 5   #8-node hexahedron
    MSHPRISM = 6 #6-node prism

    #class variables for GMSH's entities: point, curve, surface, volume
    MSH_ENTITY_NONE = -1
    MSH_ENTITY_POINT = 0
    MSH_ENTITY_CURVE = 1
    MSH_ENTITY_SURFACE = 2
    MSH_ENTITY_VOLUME = 3


    class gmshEntities:
        """ class for GMSH entities

        """

        def __init__(self, msh_entity_type, msh_entity_tag, msh_entity_boundbox,
                     msh_bounding_subentities_tags):
            self.entity_type = msh_entity_type    #int
            self.entity_tag = msh_entity_tag      #int
            self.entity_boundbox = msh_entity_boundbox  #list: [minX, minY, minZ, maxX, maxY, maxZ]
            self.entity_num_bounding_subentities = len(msh_bounding_subentities_tags) #numBoundingCurves (for surace) or numBoundngSurfaces (for volume)
            self.entity_bounding_subentities_tags = msh_bounding_subentities_tags #List: [tag1, tag2, ..., tagN]

    class gmshNodeEntityBlock:
        """ class for GMSH node entity block

        """

        def __init__(self, nodeEntityBlockDim, nodeEntityBlockTag, numNodesInBlock):
            self.entityDim = nodeEntityBlockDim
            self.entityTag = nodeEntityBlockTag
            self.numNodesInBlock = numNodesInBlock
            self.nodeCoordinates = [] #list of node coordinates in this block

    class gmshElementEntityBlock:
        """ class for GMSH element entity block

        """

        def __init__(self, elementEntityBlockDim, elementEntityBlockTag, elementType, numElementsInBlock):
            self.entityDim = elementEntityBlockDim
            self.entityTag = elementEntityBlockTag
            self.elementType = elementType
            self.numElementsInBlock = numElementsInBlock
            self.elementTags = {}  #dictionary to store {elementTag: [list of nodeTags]}


    def __init__(self):
        """Initialise Gmsh data structure"""
        self.version = "4.1" #version of MSH file format
        self.file_type = 0   #0: ASCII, 1:binary (currently only support ASCII)
        self.data_size = 8   #only used in binary file type (not relevant here)

        self.physical_names = {}  #dictionary to store physical names:
                                  # {physical_name_id: [dimension(ASCII int), physicalTag(ASCII int), "name"(127 characters max)]}
                                  # e.g. {1: [1, 3, "inlet"]}

        self.entities = {} #dictionary to store entities (only surface and volume for now). {entityID: gmshEntities object}

        self.numNodeEntityBlocks = 0 #number of node entity blocks
        self.numNodes = 0 #number of nodes

        self.nodes = {} #dictionary to store node entity blocks. {nodeEntityID: gmshNodeEntityBlock object}

        self.numElementEntityBlocks = 0 #number of element entity blocks
        self.numElements = 0 #number of elements

        self.elements = {} #dictionary to store element entity blocks. {elementEntityID: gmshElementEntityBlock object}


    def reset(self):
        """Reset gmsh2DMesh data"""
        self.physical_names = {}
        self.entities = {}
        self.numNodeEntityBlocks = 0
        self.numNodes = 0
        self.nodes = {}
        self.numElementEntityBlocks = 0
        self.numElements = 0
        self.elements = {}

    def build_from_srhgeom(self, srhgeom_obj):
        """ Build gmsh2DMesh from a SRH_2D_SRHGeom object

        Parameters
        ----------
        srhgeom_obj

        Returns
        -------

        """







    def read(self, mshfile=None):
        """Read a Gmsh .msh file.

        Reads Gmsh format 1.0 and 2.0 mesh files, storing the nodes and
        elements in the appropriate dicts.
        """

        if not mshfile:
            mshfile = open(self.filename, 'r')

        readmode = 0
        print('Reading %s' % mshfile.name)
        line = 'a'
        while line:
            line = mshfile.readline()
            line = line.strip()
            if line.startswith('$'):
                if line == '$NOD' or line == '$Nodes':
                    readmode = 1
                elif line == '$ELM':
                    readmode = 2
                elif line == '$Elements':
                    readmode = 3
                elif line == '$MeshFormat':
                    readmode = 4
                else:
                    readmode = 0
            elif readmode:
                columns = line.split()
                if readmode == 4:
                    if len(columns) == 3:
                        vno, ftype, dsize = (float(columns[0]),
                                             int(columns[1]),
                                             int(columns[2]))
                        print(('ASCII', 'Binary')[ftype] + ' format')
                    else:
                        endian = struct.unpack('i', columns[0])
                if readmode == 1:
                    # Version 1.0 or 2.0 Nodes
                    try:
                        if ftype == 0 and len(columns) == 4:
                            self.nodes[int(columns[0])] = map(float,
                                                              columns[1:])
                        elif ftype == 1:
                            nnods = int(columns[0])
                            for N in range(nnods):
                                data = mshfile.read(4 + 3 * dsize)
                                i, x, y, z = struct.unpack('=i3d', data)
                                self.nodes[i] = (x, y, z)
                            mshfile.read(1)
                    except ValueError:
                        print('Node format error: ' + line, ERROR)
                        readmode = 0
                elif ftype == 0 and readmode > 1 and len(columns) > 5:
                    # Version 1.0 or 2.0 Elements
                    try:
                        columns = map(int, columns)
                    except ValueError:
                        print('Element format error: ' + line, ERROR)
                        readmode = 0
                    else:
                        (id, type) = columns[0:2]
                        if readmode == 2:
                            # Version 1.0 Elements
                            tags = columns[2:4]
                            nodes = columns[5:]
                        else:
                            # Version 2.0 Elements
                            ntags = columns[2]
                            tags = columns[3:3 + ntags]
                            nodes = columns[3 + ntags:]
                        self.elements[id] = (type, tags, nodes)
                elif readmode == 3 and ftype == 1:
                    tdict = {1: 2, 2: 3, 3: 4, 4: 4, 5: 5, 6: 6, 7: 5, 8: 3, 9: 6, 10: 9, 11: 10}
                    try:
                        neles = int(columns[0])
                        k = 0
                        while k < neles:
                            etype, ntype, ntags = struct.unpack('=3i',
                                                                mshfile.read(3 * 4))
                            k += ntype
                            for j in range(ntype):
                                mysize = 1 + ntags + tdict[etype]
                                data = struct.unpack('=%di' % mysize,
                                                     mshfile.read(4 * mysize))
                                self.elements[data[0]] = (etype,
                                                          data[1:1 + ntags],
                                                          data[1 + ntags:])
                    except:
                        raise
                    mshfile.read(1)

        print('  %d Nodes' % len(self.nodes))
        print('  %d Elements' % len(self.elements))

        mshfile.close()

    def write_ascii(self, filename=None):
        """Dump the mesh out to a Gmsh 2.0 msh file."""

        if not filename:
            filename = self.filename

        mshfile = open(filename, 'w')

        print('$MeshFormat\n2.0 0 8\n$EndMeshFormat', file=mshfile)
        print('$Nodes\n%d' % len(self.nodes), file=mshfile)

        for node_id, coord in self.nodes.items():
            print(node_id, ' ', ' '.join([str(c) for c in coord]), sep="",
                  file=mshfile)

        print('$EndNodes', file=mshfile)
        print('$Elements\n%d' % len(self.elements), file=mshfile)

        for ele_id, elem in self.elements.items():
            (ele_type, tags, nodes) = elem
            print(ele_id, ' ', ele_type, ' ', len(tags), ' ',
                  ' '.join([str(c) for c in tags]), ' ',
                  ' '.join([str(c) for c in nodes]), sep="", file=mshfile)
        print('$EndElements', file=mshfile)

    def write_binary(self, filename=None):
        pass

    def write_geo(self, filename=None, use_ids=True, use_physicals=True):
        """Dump the mesh out to a Gmsh .geo geometry file."""

        string = ""

        for k, x in self.nodes.items():
            string += 'Point(%d) = {%f,%f,%f};\n' % (k, x[0], x[1], x[2])

        line_no = 1;
        line_ids = {}
        surface_ids = {}

        for k, l in self.elements.items():
            if l[0] == 1:
                if use_physicals:
                    line_ids.setdefault(l[1][-1], []).append(line_no)
                else:
                    line_ids.setdefault(l[1][1], []).append(line_no)
                line_ids.setdefault(l[1][-1], []).append(line_no)
                string += "Line(%d) = {%d,%d};\n" % (line_no, l[2][0], l[2][1])
                line_no += 1
            elif l[0] == 2:
                if use_physicals:
                    surface_ids.setdefault(l[1][-1], []).append(k)
                else:
                    surface_ids.setdefault(l[1][1], []).append(k)
                string += "Line(%d) = {%d,%d};\n" % (line_no,
                                                     l[2][0], l[2][1])
                string += "Line(%d) = {%d,%d};\n" % (line_no + 1,
                                                     l[2][1], l[2][2])
                string += "Line(%d) = {%d,%d};\n" % (line_no + 2,
                                                     l[2][2], l[2][0])
                string += "Line Loop(%d) = {%d,%d,%d};" % \
                          tuple([k] + range(line_no, line_no + 3))
                string += "Plane Surface(%d) = %d;" % (k, k)
                line_no += 3

        if use_ids:
            for k, v in line_ids.items():
                string += "Physical Line(%d) = {%s};\n" % (k,
                                                           ",".join(map(str, v)))
            for k, v in surface_ids.items():
                string += "Physical Surface(%d) = {%s};\n" % (k,
                                                              ",".join(map(str, v)))

        geofile = open(filename, 'w')
        geofile.write(string)
        geofile.close()

    def as_vtk(self, elementary_index=0):
        """Convert to a VTK unstructured grid object, ugrid."""

        etype = {1: vtk.VTK_LINE,
                 2: vtk.VTK_TRIANGLE}
        point_map = {};

        ugrid = vtk.vtkUnstructuredGrid()

        pts = vtk.vtkPoints()

        for i, v in enumerate(self.nodes.items()):
            (node_id, nodes) = v

            pts.InsertNextPoint(nodes)
            point_map[node_id] = i

        ugrid.SetPoints(pts)
        ugrid.Allocate(len(self.elements))

        physical_ids = vtk.vtkIntArray()
        physical_ids.SetNumberOfComponents(1)
        physical_ids.SetNumberOfTuples(len(self.elements))
        physical_ids.SetName("PhysicalIds")

        elementary_entities = vtk.vtkIntArray()
        elementary_entities.SetNumberOfComponents(1)
        elementary_entities.SetNumberOfTuples(len(self.elements))
        elementary_entities.SetName("ElementaryEntities")

        for i, v in enumerate(self.elements.items()):

            k, ele = v
            ids = vtk.vtkIdList()

            for node in ele[2]:
                ids.InsertNextId(point_map[node])

            ugrid.InsertNextCell(etype[ele[0]], ids)
            elementary_entities.SetValue(i, ele[1][elementary_index])
            physical_ids.SetValue(i, ele[1][-1])

        ugrid.GetCellData().AddArray(elementary_entities)
        ugrid.GetCellData().AddArray(physical_ids)

        return ugrid

    def from_vtk(self, ugrid):
        """Convert from a VTK unstructured grid object, ugrid."""

        etype = {vtk.VTK_LINE: 1,
                 vtk.VTK_TRIANGLE: 1}

        self.reset()

        for i in range(ugrid.GetNumberOfPoints()):
            self.nodes[i] = ugrid.GetPoint(i)

        for i in range(ugrid.GetNumberOfCells()):

            cell = ugrid.GetCell(i)
            tags = []
            ids = cell.GetPointIds()

            cd = ugrid.GetCellData()

            if cd.HasArray("ElementaryEntities"):
                tags.append(cd.GetArray("ElementaryEntities").GetValue(i))
            if cd.HasArray("PhysicalIds"):
                tags.append(cd.GetArray("PhysicalIds").GetValue(i))

            self.elements[i] = (etype[cell.GetCellType()],
                                tags,
                                [ids.GetId(_) for _
                                 in range(ids.GetNumberOfIds())])

        print('  %d Nodes' % len(self.nodes))
        print('  %d Elements' % len(self.elements))

