"""
Some misc. functions
"""
import numpy as np
from scipy import integrate

def ocf_1D_backwater_curve(slope, ManningN, startx, startH, startZ, riverLength, nGrid, specificDischarge):
    """ Open channel flow: 1D backwater curve in a wide, rectangular channel

    Parameters
    ----------
    slope: slope of channel
    ManningN: Manning's n
    startx: starting x coordinate
    startH: starting water depth H
    startZ: starting bottom elevation
    riverLength: length of river
    nGrid: number of grid to be used
    specificDischarge: specific discharge (discharge per unit width)

    Returns
    -------
    normalDepth: normal depth
    criticalDepth: critical depth
    xcoords: x coordinates of the backwater curve profile
    waterDepth: water depth along the backwater curve profile
    WSE: water surface elevation (waterDepth + bottom elevation)

    """
    # normal flow depth
    normalDepth = (ManningN * specificDischarge / np.sqrt(slope)) ** (3.0 / 5.0)

    #critical flow depth
    criticalDepth = (specificDischarge**2/9.81)**(1.0/3.0)

    print("Hn, Hc = ", normalDepth, criticalDepth)

    x = np.linspace(startx, startx+riverLength, nGrid)

    waterDepth = integrate.odeint(_F_H_backwater_curve, startH, x, args=(ManningN,slope,specificDischarge))
    waterDepth = waterDepth[:, 0]  # convert the returned 2D array to a 1D array

    #negate the x-coordinate (go from downstream to upstream)
    #negX = -x

    #water surface elevation
    WSE = waterDepth + startZ + (x-startx) * slope

    return normalDepth, criticalDepth, x, waterDepth, WSE


def _F_H_backwater_curve(H,x,n,S,qw):
    """ F(H) function for the backwater curve equation (the right hand side)

    Parameters
    ----------
    H: water depth
    x: x coordinate
    n: Manning's n
    S: channel slope
    qw: specific discharge

    Returns
    -------
    F(H)

    """

    return -(S-qw**2*n**2/(max(H,0))**(10.0/3))/(1-qw**2/9.81/(max(H,0))**3)