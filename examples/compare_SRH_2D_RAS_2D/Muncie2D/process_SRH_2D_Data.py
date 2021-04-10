import numpy as np
import os

import matplotlib.pyplot as plt

import pyHMT2D


def main():
    """ Testing SRH_2D_Data class

    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("SRH-2D/Muncie.srhhydro")

    #read SRH-2D result in XMDF format (*.h5)
    #wether the XMDF result is nodal or cell center
    bNodal = False

    my_srh_2d_data.readSRHXMDFFile("SRH-2D/Muncie_XMDFC.h5", bNodal)

    #export to VTK
    my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir='SRH-2D')


if __name__ == "__main__":
    main()

    print("All done!")