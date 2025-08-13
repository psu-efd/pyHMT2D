# Global variables

gVerbose = False
"""bool to define verbose or not"""

gMax_Nodes_per_Element = 8
"""maximum number of nodes for an element"""

gMax_Elements_per_Node = 10
"""maximum number of elements for a node"""

pyHMT2D_SCALAR = 1
"""scalar type in pyHMT2D"""

pyHMT2D_VECTOR = 2
"""vector type in pyHMT2D (with x and y components)"""

def setVerbose(bVerbose):
    """
    set the verbose level
    bVerbose: bool to define verbose or not
    """

    global gVerbose
    gVerbose = bVerbose

