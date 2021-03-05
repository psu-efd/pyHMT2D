import numpy as np

import pyHMT2D

# testing
def main():
    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("Muncie2DOnly_SI.p01.hdf", "subterrain_exported.tif")

    # print(my_ras2d.TwoDAreaFace_FacePoints[0])

    my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(timeStep=2)

    # my_ras_2d_data.exportSRHGEOMFile("Muncie2D")

    # my_ras_2d_data.exportSRHMATFile("Muncie2D")

    # my_ras_2d_data.exportBoundariesToVTK("Muncie2D_boundaries")

    #my_ras_2d_data.exportFaceProfilesToVTK("Muncie2D_faceprofiles")

    # dump all data to screen (debug)
    # my_ras_2d_data.dump_all_data()

    print("All done!")


if __name__ == "__main__":
    main()
