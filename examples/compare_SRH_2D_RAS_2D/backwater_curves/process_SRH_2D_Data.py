import numpy as np
import os

import pyHMT2D

def main():
    """ Testing SRH_2D_Data class

    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("SRH-2D/backwater.srhhydro")

    #read SRH-2D result in XMDF format (*.h5)
    #wether the XMDF result is nodal or cell center (currently it is all nodal)
    bNodal = False

    my_srh_2d_data.readSRHXMDFFile("SRH-2D/backwater_curve_XMDFC.h5", bNodal)

    #outputXMDFDataToVTK(bNodal, timeStep=-1,lastTimeStep=False, dir=''):
    my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir="SRH-2D")

    # User specified SRH-2D result in SRH (point) or SRHC (cell center) format
    #srhFileName = 'backwater_curve_SRHC24.dat'

    # whehter it is cell center data or at point (need to be checked by user)
    #bCellData = True

    # Read SRH-2D result file
    #resultVarNames, resultData = my_srh_2d_data.readSRHFile(srhFileName)
    # print(resultVarNames, resultData)

    # output SRH-2D result to VTK
    #srhName = os.path.splitext(srhFileName)[0]
    #vtkFileName = srhName + ".vtk"
    #print("vtk file name = ", vtkFileName)

    #my_srh_2d_data.outputVTK(vtkFileName, resultVarNames, resultData, bCellData)

    print("All done!")


if __name__ == "__main__":
    main()