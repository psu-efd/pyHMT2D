import sys
import numpy as np
import vtk
from vtk.util import numpy_support as VN
from scipy import interpolate
from scipy import integrate
import os.path

from pyHMT2D.Hydraulic_Models_Data import HydraulicModel

from pyHMT2D.__common__ import *

class Backwater_1D_Model(HydraulicModel):
    """Simple 1D backwater curve model .


    Attributes:


    """

    def __init__(self):
        """SRH_2D_Model constructor

        Parameters
        ----------
        """

        HydraulicModel.__init__(self, "Backwater-1D", "1.0")

        #Backwater-1D case data (an object of Backater_1D_Data)
        self._backwater_1d_data = None

        #initialize the model
        self.init_model()


    def init_model(self):
        """Initialize Backwater-1D model

        Returns
        -------

        """

        print("Initializing Backwater-1D ...")

    def set_simulation_case(self, backwater_1d_data):
        """Set the simulation case to backwater_1D_data (if it has been created already)

        Parameters
        ----------
        backwater_1d_data : Backwater_1D_Data
            an object from class Backwater_1D_Data, which should be created before calling

        Returns
        -------

        """

        self._backwater_1d_data = backwater_1d_data

    def get_simulation_case(self):
        """Get the simulation case

        Returns
        -------

        """

        if not self._backwater_1d_data:
            raise Exception("Simulation case, i.e.,, the Backwater_1D_Data object, has not been set yet. Exiting ...")

        return self._backwater_1d_data

    def run_model(self):
        """Run the Backwater-1D model

        It will run the current Backwater-1D case.

        Returns
        -------

        """

        print("Running Backwater-1D ...")

        #check whether the case has been set up
        if not self._backwater_1d_data:
            raise Exception("Backwater-1D case has not been setup yet. Call set_simulation_case(...) first. Exiting ...")

        #solve the GVF ODE equation and save the results back into self._backwater_1d_data
        self._backwater_1d_data.normalDepth, \
        self._backwater_1d_data.criticalDepth, \
        self._backwater_1d_data.gridx, \
        self._backwater_1d_data.waterDepth, \
        self._backwater_1d_data.WSE, \
        self._backwater_1d_data.gridManningN = \
                    self._ocf_1D_backwater_curve(self._backwater_1d_data.slope,
                                                 self._backwater_1d_data.ManningNFunc,
                                                 self._backwater_1d_data.startx,
                                                 self._backwater_1d_data.startH,
                                                 self._backwater_1d_data.startZ,
                                                 self._backwater_1d_data.riverLength,
                                                 self._backwater_1d_data.nGrid,
                                                 self._backwater_1d_data.specificDischarge
                                                 )


    def _ocf_1D_backwater_curve(self, slope, ManningNFunc, startx, startH, startZ,
                                riverLength, nGrid, specificDischarge):
        """ Open channel flow: 1D backwater curve in a wide, rectangular channel

        Parameters
        ----------
        slope : float
            slope of channel
        ManningNFunc : interp1d
            Interpolation function for Manning's n
        startx : float
            starting x coordinate
        startH : float
            starting water depth H
        startZ : float
            starting bottom elevation
        riverLength : float
            length of river
        nGrid : int
            number of grid to be used
        specificDischarge : float
            specific discharge (discharge per unit width)

        Returns
        -------
        normalDepth : float
            normal depth (only makes sense for uniform Manning's n)
        criticalDepth: float
            critical depth (only makes sense for uniform Manning's n)
        xcoords : numpy array
            x coordinates of the backwater curve profile
        waterDepth : numpy array
            water depth along the backwater curve profile
        WSE : numpy array
            water surface elevation (waterDepth + bottom elevation)

        """
        # normal flow depth (only makes sense for uniform Manning's n value;
        # here we use the Manning's n at x = startx as a constant)
        const_ManningN = ManningNFunc(startx)
        normalDepth = (const_ManningN * specificDischarge / np.sqrt(slope)) ** (3.0 / 5.0)

        # critical flow depth
        criticalDepth = (specificDischarge ** 2 / 9.81) ** (1.0 / 3.0)

        print("Hn, Hc = ", normalDepth, criticalDepth)

        x = np.linspace(startx, startx + riverLength, nGrid)

        #calculate Manning's n at grid points
        ManningN = ManningNFunc(x)

        #units system factor
        Kn = 1.486 if self._backwater_1d_data.units == "EN" else 1.0

        waterDepth = integrate.odeint(self._F_H_backwater_curve, startH, x, args=(Kn, ManningNFunc, slope, specificDischarge))
        waterDepth = waterDepth[:, 0]  # convert the returned 2D array to a 1D array

        # negate the x-coordinate (go from downstream to upstream)
        # negX = -x

        # water surface elevation
        WSE = waterDepth + startZ + (x - startx) * slope

        return normalDepth, criticalDepth, x, waterDepth, WSE, ManningN

    def _F_H_backwater_curve(self, H, x, Kn, ManningNFunc, S, qw):
        """ F(H) function for the backwater curve equation (the right hand side)

        Parameters
        ----------
        H: water depth
        x: x coordinate
        Kn: units system factor (1.486 for EN, 1.0 for SI)
        ManningNFunc: Interpolation function for Manning's n
        S: channel slope
        qw: specific discharge

        Returns
        -------
        F(H)

        """

        return -(S - qw ** 2 * (ManningNFunc(x) / Kn) ** 2 / (max(H, 0)) ** (10.0 / 3)) / \
               (1 - qw ** 2 / 9.81 / (max(H, 0)) ** 3)