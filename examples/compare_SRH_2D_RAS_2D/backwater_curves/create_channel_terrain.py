"""
This script prepares the files needed for a simple backwater curve case. The case is 2D.

Steps:
1. prepare the terrain based on given slope and domain extent
2. prepare 2D channel domain geometry

For georeferencing, an arbitrary projected coordinate system, EPSG:32128, NAD83 Pennsylvania North
is used.
"""


from pyHMT2D.Misc import Terrain

def main():
    """

    Returns
    -------

    """

    #define parameters
    slope = 1e-5  #channel streamwise slope

    channel_lenx = 10*1e3 #channel streamwise length (meter)
    channel_leny = 1000     #channel width (meter)
    extra_len = 1000 #extra len in both directions to have some free room

    #pixel width and height (each pixel represents how many meters/feet in real world)
    #the smaller these numbers, the higher the resolution of the terrain and the larger
    #the GeoTiff file.
    pixel_width = 5.0
    pixel_height = 5.0

    #elevation at the top-left corner (origin)
    elevation_origin = 0.0

    #create the Terrain object
    my_terrain = Terrain("constant_slope_channel")

    #create the elevation data
    my_terrain.create_constant_slope_channel_elevation(slope,channel_lenx,channel_leny,
                                                       pixel_width,pixel_height,elevation_origin,
                                                       extra_len)

    #set georeferencing

    #raster's top-left corner georeferenced location (hypothetical; for fun it is the
    #middle of the Beaver Stadium at Penn State)
    geoTopLeft_x = 259101.0
    geoTopleft_y = 4521837.0

    #print out the coordinates of the four corners (to be copy and pasted into HEC-RAS geometry definition)
    #these corners have to be slightly shifted to be inside the bathymetry
    shift = 500.0
    print("{:16.9f} {:16.8f}".format(geoTopLeft_x+shift, geoTopleft_y-shift))
    print("{:16.9f} {:16.8f}".format(geoTopLeft_x+shift+channel_lenx, geoTopleft_y-shift))
    print("{:16.9f} {:16.8f}".format(geoTopLeft_x+shift+channel_lenx, geoTopleft_y-shift-channel_leny))
    print("{:16.9f} {:16.8f}".format(geoTopLeft_x+shift, geoTopleft_y-shift-channel_leny))
    print("{:16.9f} {:16.8f}".format(geoTopLeft_x+shift, geoTopleft_y-shift))

    print("Channel real length in x = ", channel_lenx)
    print("Channel real length in y = ", channel_leny)

    #set coordiante reference system's EPSG code
    EPSGCode = 32128 # NAD83 Pennsylvania North (unit = meter)
    my_terrain.set_georeference(geoTopLeft_x, geoTopleft_y, EPSGCode)

    #save the terrain to file
    geoTiffFileName = "constant_slope_channel.tif"

    my_terrain.save_terrain_to_file(geoTiffFileName,'GTiff')

    print("All done!")


if __name__ == "__main__":
    main()