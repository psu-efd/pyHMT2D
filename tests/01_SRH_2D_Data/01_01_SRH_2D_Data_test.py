import pyHMT2D

def test_SRH_2D_Data():
    """Test the SRH_2D_Data object

    Returns
    -------

    """

    #pytest run from the project level. So here the data path is test/data
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("tests/data/Muncie.srhhydro")

    print(my_srh_2d_data.srhhydro_obj.srhhydro_content["SimTime"])
    #print(my_srh_2d_data.srhhydro_obj.srhhydro_content)

    assert len(my_srh_2d_data.srhhydro_obj.srhhydro_content["SimTime"]) == 3


if __name__ == "__main__":

    test_SRH_2D_Data()