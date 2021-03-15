"""
SRH_2D_Model:
"""

import win32com.client as win32
from pyHMT2D.__common__ import HydraulicModel

from .helpers import *
import sys

class SRH_2D_Model(HydraulicModel):
    """SRH_2D Model

    SRH_2D_Model controls the run of SRH_2D. A plan can be loaded and executed.

    Currently SRH_2D_Model has very limited capability to modify geometry, mesh, and flow data.

    Attributes:
        _faceless: {bool} -- whether show the SRH-2D window when it is running
        _srh_path: path to the SRH_2D program (e.g., srh_2d.exe depending on the version)
        _srh_pre_path: path to the SRH_2D preprocessing program (e.g., srh_2d_pre.exe depending on the version)

        _project_file_name: HEC-RAS project name (including the full path)
        _project: HEC_RAS_Project object corresponding to current project

    """
    def __init__(self, version, faceless=False):
        HydraulicModel.__init__(self, "SRH-2D", version)

        #whether run HEC-RAS without GUI interface showing up
        self._faceless = faceless