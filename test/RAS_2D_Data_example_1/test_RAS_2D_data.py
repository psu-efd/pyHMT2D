import numpy as np

import pyHMT2D

# testing
def main():
    my_ras2d = pyHMT2D.RAS_2D.RAS_2D_Data("Muncie2DOnly_SI.p01.hdf", "subterrain_exported.tif")

    # print(my_ras2d.TwoDAreaFace_FacePoints[0])

    my_ras2d.saveHEC_RAS2D_results_to_VTK(timeStep=2)

    # my_ras2d.exportSRHGEOMFile("Muncie2D")

    # my_ras2d.exportSRHMATFile("Muncie2D")

    # my_ras2d.exportBoundariesToVTK("Muncie2D_boundaries")

    #my_ras2d.exportFaceProfilesToVTK("Muncie2D_faceprofiles")

    # dump all data to screen (debug)
    # my_ras2d.dump_all_data()

    print("All done!")


if __name__ == "__main__":
    main()
