import pyHMT2D

def demo_SRH_2D_Data():
    """ Demonstrate the use of the SRH_2D_Data class to process SRH-2D result data

    Returns
    -------

    """

    #create the SRH_2D_Data object
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie.srhhydro")

    #read SRH-2D result in XMDF format (*.h5)
    #wether the XMDF result is nodal or cell center
    bNodal = False

    my_srh_2d_data.readSRHXMDFFile("Muncie_XMDFC.h5", bNodal)

    #export to VTK
    my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True)

if __name__ == "__main__":

    demo_SRH_2D_Data()

    print("All done!")