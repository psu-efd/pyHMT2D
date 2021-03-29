from . import (
    HydraulicModel,
    SRH_2D,
    RAS_2D,
    Misc
)
from .__about__ import __version__
from .__common__ import *

__all__ = [
    "HydraulicModel",
    "SRH_2D",
    "RAS_2D",
    "Misc",
    "gMax_Nodes_per_Element",
    "gMax_Elements_per_Node",
    "__version__"
]

__version_info__ = (1, 0, 0)
__version__ = '.'.join(map(str, __version_info__))


