import pyHMT2D

def test_RAS_2D_Data_construction(datafiles):
    """Test the construction of SRH_2D_Data object from HDF and terrain file

    Returns
    -------

    """

    path = str(datafiles)  # Convert from py.path object to path (str)
    print(path)

    #pytest run from the project level. So here the data path is test/data
    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("tests/data/Muncie2D.p01.hdf", "tests/data/TerrainMuncie_composite.tif")

    #print(my_ras_2d_data.TwoDAreaFace_FacePoints[0][0] )

    assert my_ras_2d_data.TwoDAreaFace_FacePoints[0][0][0] == 3
    assert my_ras_2d_data.TwoDAreaFace_FacePoints[0][0][1] == 2


