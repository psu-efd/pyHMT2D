"""
A Python class for terrain data I/O, creation, and manipulation
"""

import numpy as np
import sys

from osgeo import gdal
from osgeo import osr
from osgeo import ogr

from pyHMT2D.Hydraulic_Models_Data import HydraulicData

class Terrain(HydraulicData):
    """A Python class for terrain data I/O, creation, and manipulation

    Typical work flow is as follows:

    1. create the Terrain object
    2. create the terrain (elevation, pixel_width/height): user can either call
       some pre-defined terrains such as constant slope, or create the terrain by themsleves
       and then call set_terrain(...), and set_pixel_size(...)
    3. set the georeferencing by calling set_georeference(...)
    4. save the terrain to file by calling save_terrain_to_file(...)

    Attributes
    ----------
    name : str
        name of the terrain
    elevation : numpy.ndarray
        elevation of the terrain (2D numpy array)
    geoTopLeft_x : float
        raster's top-left corner georeferenced x-coordinate
    geoTopLeft_y : float
        raster's top-left corner georeferenced y-coordinate
    pixel_width : float
        pixel width (in x), i.e, each pixel is how many meters/feet wide?
    pixel_height : float
        pixel height (in y), i.e, each pixel is how many meters/feet high?
    EPSGCode : int
        EPSG (European Petroleum Survey Group) code that defines the coordinate reference system
    geoTransform : list
        affine transform for the raster image from image pixel to real coordinates
    supportedGDALDrivers : list
        list of supported GDAL drivers (short names only)

    """

    def __init__(self, name):
        """Terrain class constructor

        Parameters
        ----------
        name : str
            name of the terrain

        """

        HydraulicData.__init__(self, "Terrain")

        self.name = name      # name of the terrain
        self.elevation = np.array([]) # elevation of the terrain, should be 2D numpy array

        self.geoTopLeft_x = 0.0 # raster's top-left corner georeferenced location
        self.geoTopleft_y = 0.0

        self.pixel_width = -1.0 # pixel width (in x), i.e, each pixel is how many meters/feet wide?
        self.pixel_height = -1.0 # pixel height (in y)

        self.EPSGCode = -1 # EPSG (European Petroleum Survey Group) code that defines the coordinate reference system

        self.geoTransform = [] # affine transform for the raster image from image pixel to real coordinates

        self.supportedGDALDrivers = [] #list of supported GDAL drivers (short names only)
        self.build_supported_GDAL_drivers_list()


    def get_terrain_name(self):
        """Get the terrain name

        Returns
        -------
        name : str
            name of the terrain

        """

        return self.name

    def get_elevation(self):
        """Get the elevation array

        Returns
        -------
        elevation : numpy.array
            elevation 2D array
        """

        return self.elevation

    def set_elevation(self, elevation):
        """Set the elevation array

        Parameters
        _______
        elevation : numpy.ndarray
            elevation 2D array

        Returns
        -------

        """

        self.elevation = elevation

    def set_pixel_size(self, pixel_width, pixel_height):
        """Set the pixel size in x and y directions

        Parameters
        ----------
        pixel_width : float
            the width of each pixel in real world
        pixel_height : float
            the height of each pixel in real world

        Returns
        -------

        """

        self.pixel_width = pixel_width
        self.pixel_height = pixel_height

    def build_supported_GDAL_drivers_list(self):
        """Build the list of supported GDAL drivers list

        Returns
        -------

        """

        for i in range(gdal.GetDriverCount()):
            drv = gdal.GetDriver(i)
            self.supportedGDALDrivers.append(drv.ShortName)

        print("Supported GDAL drivers are: ", self.supportedGDALDrivers)

    def create_constant_slope_channel_elevation(self, slope, channel_lenx, channel_leny, pixel_width,
                                              pixel_height, elevation_origin=0, extra_len=0):
        """Create a constant slope channel elevation

        The slope is in the x direction only. The slope in the y direction is zero.

        Parameters
        -----------
        slope : float
            slope in x
        channel_lenx : float
            channel length in x
        channel_leny : float
            channel length in y
        elevation_origin : float
            the elevation at the origin (top left; does not account for the extra fringe)
        extra_len : float
            optional extra fringe length added to the channel domain (to
            have some free room in developing 2D models in e.g., SMS or HEC-RAS.

        Returns
        -------

        """

        # raster image width and height in pixels
        self.pixel_width = pixel_width
        self.pixel_height = pixel_height

        #raster size (without extra length)
        nx = int(channel_lenx/pixel_width)
        ny = int(channel_leny/pixel_height)

        #extra length (fringe)
        extra_len_nx = int(extra_len/pixel_width)
        extra_len_ny = int(extra_len/pixel_height)
        self.elevation = np.zeros((ny+extra_len_ny*2,nx+extra_len_nx*2))

        #set the elevation
        for iy in range(0, ny + extra_len_ny * 2):
            for ix in range(0, nx + extra_len_nx * 2):
                self.elevation[iy, ix] = elevation_origin -slope * (ix-extra_len_nx) * pixel_width

    def set_georeference(self, geoTopLeft_x, geoTopLeft_y, EPSGCode):
        """Set the georeferencing information

        Parameters
        ----------
        geoTopLeft_x : float
            raster's top-left corner georeferenced x location
        geoTopLeft_y : float
            raster's top-left corner georeferenced y location
        EPSGCode : int
            EPSG (European Petroleum Survey Group) code that defines the coordinate reference system

        Returns
        -------

        """

        if self.pixel_width < 0 or self.pixel_height < 0:
            print("Either pixel_width or pixel_height is negative. Need to create the terrain first.")
            print("Exiting...")
            sys.exit()

        self.geoTopLeft_x = geoTopLeft_x
        self.geoTopleft_y = geoTopLeft_y

        self.geoTransform = [geoTopLeft_x, self.pixel_width, 0, self.geoTopleft_y, 0, -self.pixel_height]

        self.EPSGCode = EPSGCode

    def save_terrain_to_file(self, terrainFileName, geoDriverName = 'GTiff'):
        """save terrain to file, such as GeoTiff

        Parameters
        -----------
        terrainFileName : str
            file name for the saved terrain
        geoDriverName : str
            GDAL raster driver names, such as 'GTiff'

        """

        print("Save the terrain to file:", terrainFileName)

        if geoDriverName not in self.supportedGDALDrivers:
            print("The provided geoDriverName", geoDriverName, "is not supported.")
            print("Supported GDAL driver names are:", self.supportedGDALDrivers)
            print("Exiting...")
            sys.exit()

        if (not self.elevation.size) or (self.pixel_width < 0) or (self.pixel_height < 0) \
                or (self.EPSGCode < 0) or (len(self.geoTransform) == 0):
            print("Terrain data is not complete. Missing elevation, pixel size, EPSG code, "
                  "or geoTransform.")
            print("Exiting...")
            sys.exit()

        driver = gdal.GetDriverByName(geoDriverName)
        dst_filename = terrainFileName
        dst_ds = driver.Create(dst_filename, self.elevation.shape[1],
                               self.elevation.shape[0], 1, gdal.GDT_Float32)

        dst_ds.SetGeoTransform(self.geoTransform)

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.EPSGCode)

        dst_ds.SetProjection(srs.ExportToWkt())

        dst_ds.GetRasterBand(1).WriteArray(self.elevation)

        # Once we're done, close properly the dataset
        dst_ds = None


