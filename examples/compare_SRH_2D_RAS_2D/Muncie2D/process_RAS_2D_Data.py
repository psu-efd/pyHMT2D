
import pyHMT2D

def main():

    #build the RAS_to_SRH_Converter object
    my_ras_to_srh_converter = pyHMT2D.Misc.RAS_to_SRH_Converter("HEC-RAS/Muncie2D.p01.hdf",
                                                                "HEC-RAS/Terrain/TerrainMuncie_composite.tif",
                                                                "SRH-2D/Muncie")

    #convert to SRH-2D (mesh and material) and save
    print("Convert HEC-RAS 2D mesh and material to SRH-2D.")
    my_ras_to_srh_converter.convert_to_SRH()

    #also convert HEC-RAS 2D result data to VTK (for comparision with SRH-2D)
    print("Convert HEC-RAS 2D results to VTK.")
    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("HEC-RAS/Muncie2D.p01.hdf",
                                         "HEC-RAS/Terrain/TerrainMuncie_composite.tif")

    my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True, dir="HEC-RAS")


if __name__ == "__main__":
    main()

    print("All done!")
