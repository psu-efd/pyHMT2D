"""
A Python class for terrain data I/O, creation, and manipulation
"""

import math

import numpy as np
import sys

from osgeo import gdal
from osgeo import osr

from pyHMT2D.Hydraulic_Models_Data import HydraulicData
import random

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

    def add_bedform(self, channel_lenx, channel_leny, Lb, Hb, alpha_lee, a_stoss, rotation=0,perturbation=0):
        """Add bedform feature to the terrain

        Typical use scenario is firstly to create a constant slope channel and then add the bedform features.

        Parameters
        ----------
        channel_lenx : float
            channel length in x
        channel_leny : float
            channel length in y
        Lb : float
            bed form length
        Hb : float
            bed form height
        alpha_lee : float
            bed-form's lee side slope angle (in degrees)
        a_stoss : float
            bed-form's stoss size sine function amplitude
        rotation : float
            optional rotation angle of the domain in degrees (default is zero)
        perturbation : float
            optional perturbation added to the terrain (default is zero)

        Returns
        -------

        """

        #lee side length
        Llee = Hb/np.tan(math.radians(alpha_lee))

        #stoss side length
        Lstoss = Lb - Llee

        if Lstoss <=0:
            raise Exception("The calculated stoss lengh is negative. Check the bedform parameters.")

        print("Bedform lengths of stoss and lee sides are: ", Lstoss, Llee)

        #stoss side angle
        alpha_stoss = np.arctan(Hb/Lstoss)

        #raster size (without extra length)
        nx = int(channel_lenx/self.pixel_width)
        ny = int(channel_leny/self.pixel_height)

        # set the elevation
        for iy in range(self.elevation.shape[0]):
            for ix in range(self.elevation.shape[1]):

                # get current grid point's coordinate
                x = ix * self.pixel_width
                y = iy * self.pixel_height
                L = np.sqrt(x**2+y**2)

                x_new = x
                y_new = y

                # take care of optional rotation
                if np.abs(rotation) > 1e-3:
                    alpha = np.arctan(y/(x+1e-6))
                    alpha_prime = alpha + np.deg2rad(rotation)
                    x_new = L*np.cos(alpha_prime)
                    y_new = L*np.sin(alpha_prime)

                # get the location of current grid point within one bedform
                xprime = x_new % Lb

                # elevation due to the existence of bedform
                deltaZ = 0.0

                if xprime <= Lstoss:  #if current location is in the stoss side
                    deltaZ = xprime*np.tan(alpha_stoss) - a_stoss*np.sin(2*math.pi*xprime/Lstoss)
                else:                 #if the current location is in the lee side
                    deltaZ = Hb - Hb*(xprime-Lstoss)/(Lb-Lstoss)

                # take care of optional perturbation
                if np.abs(perturbation) > 1e-3:
                    deltaZ += (random.random() - 0.5)*2*perturbation

                self.elevation[iy, ix] += deltaZ

    def add_composite_channel(self, channel_lenx, channel_leny, L1, B, D, alpha):
        """Add composite channel to the terrain

        Typical use scenario is firstly to create a constant slope channel and then add the composite channel.

        Parameters
        ----------
        channel_lenx : float
            channel length in x
        channel_leny : float
            channel length in y
        L1 : float
            flood plain's width on one side (the other side can be calculated)
        B : float
            main channel bottom width
        D : float
            main channel depth
        alpha : float
            main channel's side slope angle (in degrees)

        Returns
        -------

        """

        #main channel slope's horizontal distance
        Lside = D*np.tan(math.radians(alpha))

        #raster size (without extra length)
        nx = int(channel_lenx/self.pixel_width)
        ny = int(channel_leny/self.pixel_height)

        # set the elevation
        for iy in range(self.elevation.shape[0]):
            for ix in range(self.elevation.shape[1]):

                # get current grid point's coordinate
                x = ix * self.pixel_width
                y = iy * self.pixel_height

                # elevation modification due to the existence of main channel
                deltaZ = 0.0

                if (y>L1 and y <= (L1+Lside)):
                    deltaZ = -(y-L1)/np.tan(math.radians(alpha))
                elif (y>(L1+Lside) and y<= (L1+Lside+B)):
                    deltaZ = -D
                elif (y>(L1+Lside+B) and y < (L1+Lside+B+Lside)):
                    deltaZ = -(L1+Lside+B+Lside - y)/np.tan(math.radians(alpha))
                else:
                    deltaZ = 0.0

                self.elevation[iy, ix] += deltaZ


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


