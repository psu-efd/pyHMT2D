"""
Some tools
"""

import vtk
import sys

import numpy as np

from osgeo import gdal
from osgeo import osr
from osgeo import ogr

def setNumpyArrayValueToNaN(array, value):
    # for numpy array value to NaN if it equals to value
    array[array==value]=np.nan

    return array



def printProgressBar(i,total,postText):
    """A simple function to print a progress bar (no need to use other libraries)

    Note: any print inbetween the calling of this function will create a new progress bar. So
    it is better not print anything in between the calling.

    Reference:
    https://stackoverflow.com/questions/3002085/python-to-print-out-status-bar-and-percentage

    Parameters
    ----------
    i
    total
    postText

    Returns
    -------

    """
    n_bar =20 #size of progress bar
    j= (i+1)/total
    j=max(j,0.01001)  #at least print something at the beginning
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%  {postText}")
    sys.stdout.flush()

def build_gdal_vrt(vrtFileName, sourceGeoTiffFileNameList):
    """build GDAL VRT to compose terrain from multiple GeoTiff files

    Note: GDAL's earlier versions (<2.3.3) can't read BigTIFF files.

    It is an equivalent of the GDAL command: gdalbuildvrt result.vrt *.tif
    https://gdal.org/programs/gdalbuildvrt.html

    Attributes
    -------
    vrtFileName: {string} -- result VRT file name
    sourceGeoTiffFileNameList: {list} -- list of source GeoTiff file names

    Returns
    -------

    """

    vrt_options = gdal.BuildVRTOptions(resampleAlg='cubic', addAlpha=True)

    #result_vrt = gdal.BuildVRT(vrtFileName, sourceGeoTiffFileNameList, options=vrt_options)
    result_vrt = gdal.BuildVRT(vrtFileName, sourceGeoTiffFileNameList)

    #write the vrt to file
    result_vrt = None
