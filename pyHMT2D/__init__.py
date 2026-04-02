
from .__about__ import __version__
from .__common__ import *

__all__ = [
    "cli",
    "Hydraulic_Models_Data_Base",
    "Backwater_1D",
    "SRH_2D",
    "RAS_2D",
    "Misc",
    "gVerbose",
    "gMax_Nodes_per_Element",
    "gMax_Elements_per_Node",
    "__version__"
]


def __getattr__(name):
    """Lazy-load heavy submodules on first access.

    Replaces the previous eager ``from .Hydraulic_Models_Data import *``
    pattern so that importing any lightweight sub-package of pyHMT2D
    (e.g. ``pyHMT2D.AI_Tools``) does not pull in vtk or pywin32 until
    they are actually needed.
    """
    import importlib

    _SUBMODULE_MAP = {
        "SRH_2D":                     "pyHMT2D.Hydraulic_Models_Data.SRH_2D",
        "RAS_2D":                     "pyHMT2D.Hydraulic_Models_Data.RAS_2D",
        "Backwater_1D":               "pyHMT2D.Hydraulic_Models_Data.Backwater_1D",
        "Hydraulic_Models_Data_Base": "pyHMT2D.Hydraulic_Models_Data.Hydraulic_Models_Data_Base",
        "Misc":                       "pyHMT2D.Misc",
        "cli":                        "pyHMT2D.cli",
    }
    if name in _SUBMODULE_MAP:
        mod = importlib.import_module(_SUBMODULE_MAP[name])
        globals()[name] = mod  # cache so __getattr__ is not called again
        return mod

    raise AttributeError(f"module 'pyHMT2D' has no attribute {name!r}")
