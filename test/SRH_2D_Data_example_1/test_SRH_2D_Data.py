import numpy as np
import os

import pyHMT2D

def use_SRHC_data():
    """ Testing SRH_2D_data class

    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie2D.srhhydro", "Muncie2D.srhgeom", "Muncie2D.srhmat")

    # User specified SRH-2D result in SRH (point) or SRHC (cell center) format
    srhFileName = 'Muncie2D_SRHC2.dat'

    # whehter the data is nodal (True) or at cell center (False) (need to be correctly set by user)
    bNodal = False

    # Read SRH-2D result file
    resultVarNames, resultData = my_srh_2d_data.readSRHFile(srhFileName)
    print(resultVarNames)
    print(resultData)
    print(resultData.shape)
    print(resultData[0])

    #copy the data to a numpy array
    resultDataArray = np.zeros((resultData.shape[0], len(resultData[0])), dtype="float32")
    for i in range(resultData.shape[0]):
        for j in range(len(resultData[0])-1):
            resultDataArray[i,j] = resultData[i][j]

    print(resultDataArray.shape)
    print("resultDataArray = \n", resultDataArray)

    # output SRH-2D result to VTK
    srhName = os.path.splitext(srhFileName)[0]
    vtkFileName = srhName + ".vtk"
    print("vtk file name = ", vtkFileName)

   # Note: in resultData, the 1st, second, third, and fourth columns are Point_ID,
    #      X_unit, Y_unit, and Bed_Elev_unit respectively. There is no need to save them
    #      in vtk.
    my_srh_2d_data.outputVTK(vtkFileName, resultVarNames[4:], resultDataArray[:,4:], bNodal)

def use_XMDF_data():
    """ Testing SRH_2D_data class

    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie2D.srhhydro", "Muncie2D.srhgeom", "Muncie2D.srhmat")

    # User specified SRH-2D result in XMDF format
    xmdfFileName = 'Muncie2D_XMDFC.h5'

    # whehter the data is nodal (True) or at cell center (False) (need to be correctly set by user)
    bNodal = False

    #read XMDF file
    my_srh_2d_data.readSRHXMDFFile(xmdfFileName, bNodal)

    # output SRH-2D result to VTK
    my_srh_2d_data.outputXMDFDataToVTK(bNodal)

def test_srhhydro_file():
    """Test the reading and writing of SRHHYDRO file

    Returns
    -------

    """
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie2D.srhhydro", "Muncie2D.srhgeom", "Muncie2D.srhmat")

    #test changing Manning's n
    my_srh_2d_data.srhhydro_obj.modify_ManningsN(1,0.03)

    my_srh_2d_data.srhhydro_obj.write_to_file("new.srhhydro")


if __name__ == "__main__":
    #use_SRHC_data()

    #use_XMDF_data()

    test_srhhydro_file()

    print("All done!")