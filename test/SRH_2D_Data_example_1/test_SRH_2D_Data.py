import numpy as np
import os

import pyHMT2D

def main():
    """ Testing SRH_2D_data class

    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie2D.srhhydro", "Muncie2D.srhgeom", "Muncie2D.srhmat")

    # User specified SRH-2D result in SRH (point) or SRHC (cell center) format
    srhFileName = 'Muncie2D_SRHC2.dat'

    # whehter it is cell center data or at point (need to be checked by user)
    bCellData = True

    # Read SRH-2D result file
    resultVarNames, resultData = my_srh_2d_data.readSRHFile(srhFileName)
    # print(resultVarNames, resultData)

    # output SRH-2D result to VTK
    srhName = os.path.splitext(srhFileName)[0]
    vtkFileName = srhName + ".vtk"
    print("vtk file name = ", vtkFileName)

    my_srh_2d_data.outputVTK(vtkFileName, resultVarNames, resultData, bCellData)

    print("All done!")


if __name__ == "__main__":
    main()