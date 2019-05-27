# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 15:57:46 2017

@author: ciullo
"""

from osgeo import gdal, gdalconst
import numpy as np
import matplotlib.pyplot as plt

# from Stefano Bagli - stebubu

# -------------------------------------------------------------------------------
#   GDAL2Numpy
# -------------------------------------------------------------------------------


def GDAL2Numpy(pathname):
    dataset = gdal.Open(pathname, gdalconst.GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    geotransform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    wdata = band.ReadAsArray(0, 0, cols, rows).astype("float32")
    nodata = band.GetNoDataValue()
    return (wdata, geotransform, projection, nodata)
# -------------------------------------------------------------------------------
#   Numpy2GTiff
# -------------------------------------------------------------------------------


def Numpy2GTiff(arr, geotransform, projection, filename):
    if isinstance(arr, np.ndarray):
        rows, cols = arr.shape
        if rows > 0 and cols > 0:
            dtype = str(arr.dtype)
            if dtype in ["uint8"]:
                fmt = gdal.GDT_Byte
            elif dtype in ["uint16"]:
                fmt = gdal.GDT_UInt16
            elif dtype in ["uint32"]:
                fmt = gdal.GDT_UInt32
            elif dtype in ["float32"]:
                fmt = gdal.GDT_Float32
            elif dtype in ["float64"]:
                fmt = gdal.GDT_Float64
            else:
                fmt = gdal.GDT_Float64

            driver = gdal.GetDriverByName("GTiff")
            dataset = driver.Create(filename, cols, rows, 1, fmt)
            if (geotransform is not None):
                dataset.SetGeoTransform(geotransform)
            if (projection is not None):
                dataset.SetProjection(projection)
            dataset.GetRasterBand(1).WriteArray(arr)
            dataset = None
    return None

# from Joel Lawhead in Learning Geospatial Analysis with Python - Second
# Edition


def floodFill(c, r, mask):
    """
    Crawls a mask array containing
    only 1 and 0 values from the
    starting point (c=column,
    r=row - a.k.a. x, y) and returns
    an array with all 1 values
    connected to the starting cell.
    This algorithm performs a 4-way
    check non-recursively.
    """
    # cells already filled
    filled = set()
    # cells to fill
    fill = set()
    fill.add((c, r))
    width = mask.shape[1] - 1
    height = mask.shape[0] - 1
    # Our output inundation array
    flood = np.zeros_like(mask, dtype=np.int8)
    # Loop through and modify the cells which
    # need to be checked.

    while fill:
        # Grab a cell
        x, y = fill.pop()  # pop returns what is in the set and empties the set
        if y == height or x == width or x < 0 or y < 0:
            # Don't fill
            continue
        if mask[y][x] == 1:  # if you already filled it
            # Do fill
            flood[y][x] = 1
            filled.add((x, y))
            # Check neighbors for 1 values
            west = (x - 1, y)
            east = (x + 1, y)
            north = (x, y - 1)
            south = (x, y + 1)
            if west not in filled:
                fill.add(west)    # by stating so you keep doing the 'while'
            if east not in filled:
                fill.add(east)
            if north not in filled:
                fill.add(north)
            if south not in filled:
                fill.add(south)
    return flood


def plot_depthmap(dem, wlevel, coordinates):
        # dem: dem with terrain + dike ring with right heights
    x_ground, y_ground = coordinates
    a = np.where(dem < wlevel, 1, 0)
    # fld input is col,row,mask
    fld = floodFill(x_ground, y_ground, a)

    wd_map = (fld * wlevel - fld * dem) / 100.0  # in m

    plt.figure(figsize=(10, 10))
    plt.imshow(wd_map)
    plt.colorbar()
