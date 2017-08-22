#!/usr/bin/env python

"""
Utilities

The code is still under dubugging.


"""
import os
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
import numpy as np
from obspy import read
import operator

def read_xyz(file):
    """
    Read XYZ data from file.

    Usage: points, values = read_xyz(file)

    :param file: The file contains data in XYZ format. Each row contains three column Lon, Lat, Depth
    :return: points, values
    """

    with open(file, "r") as fp:
        lst = fp.readlines()

    points2 = []
    values2 = []
    for line in lst:
        row = line.split()
        try:
            lon = float(row[0])
            lat = float(row[1])
            dep = float(row[2])

            points2.append([lon, lat])
            values2.append(dep)
        except:
            continue

    points = np.array(points2)
    values = np.array(values2)

    return points, values

def read_xyzv(file):

    with open(file, "r") as fp:
        lst = fp.readlines()

    points = []
    values = []

    lata = []
    lona = []
    depa = []
    # vela = []

    mylist = []

    for i in range(0, len(lst)):
        row = lst[i].split()

        lon = float(row[0])
        lat = float(row[1])
        dep = float(row[2])
        vel = float(row[3])

        mylist.append([lon, lat, dep, vel])

    mylist.sort(key=operator.itemgetter(0, 1, 2))


    for line in mylist:
        lon = line[0]
        lat = line[1]
        dep = line[2]
        vel = line[3]
        points.append([lon, lat, dep])
        values.append(vel)

        if dep not in depa:
            depa.append(dep)

        if lat not in lata:
            lata.append(lat)

        if lon not in lona:
            lona.append(lon)

    # points = np.array(points)
    values = np.array(values)

    lata = np.array(lata)
    lona = np.array(lona)
    depa = np.array(depa)

    nlon = len(lona)
    nlat = len(lata)
    ndep = len(depa)

    if nlon * nlat * ndep != len(values):
        print "Not regular grids. nlon=%d, nlat=%d, ndep=%d, nvel=%d" % (nlon, nlat, ndep, len(values))
        os._exit(-1)

    # print lona.shape, lata.shape, depa.shape
    points = (lona, lata, depa)
    values = values.reshape(nlon, nlat, ndep)

    return points, values

def interp2(points, values, lon, lat, method="linear"):
    """
    Interpolate XYZ at any point.

    :param points: Data point coordinates. Can either be an array of shape (n, D), or a tuple of ndim arrays.
    :param values: Data values.
    :param lon:
    :param lat:
    :param method: linear, nearest, cubic
    :return:
    """
    zval  = griddata(points, values, (lon, lat), method=method)
    return zval

def interp3(points, values, lon, lat, dep):

    # construct interpolator
    my_interp3d = RegularGridInterpolator(points, values)
    pt = np.array([[lon, lat, dep]])
    val = my_interp3d(pt)
    return val[0]

def constuct_1D_data_from_3D(coors, file3d, path1d=".", depmin=5, depmax=300, depint=5, suffix=".vp"):
    """
    Construct 1-D data (e.g., 1-D velocity profile along depth) from a 3-D data cube


    Parameters
    ----------
    coors : list
        containing station id, station longitude, and station latitude
        [[AU.FORT.10, 110.0, -20.0], [CH.NONE, 120.1234, 25.1234]]

    file3d: string
        file name of 3D datacube in xyzv format

    path1d : string
        path to save constructed 1-D data
        default is "."

    depmin : float or int
        min depth for interpolation

    depmax : float or int
        max depth of interpolation

    depint : float or int
        depth interval of interpolation

    Notes
    -----
    #. depmin and depmax should be in the range of the gridfile.

    """

    try:
        os.makedirs(path1d)
    except:
        pass

    # read 3-D data
    points, values = read_xyzv(file=file3d)
    # print pnt[0]


    for line in coors:
        id = line[0]
        stlo = line[1]
        stla = line[2]
        print id

        fn = path1d + "/" + id + suffix
        fp = open(fn, "w")
        for depth in np.arange(depmin, depmax, depint):
            val = interp3(points=points, values=values, lon=stlo, lat=stla, dep=depth)
            fp.write("%s\t%f\t%f\t%f\t%f\n" % (id, stlo, stla, depth, val))

        fp.close()

    pass

def get_stid_sac(file):
    """
    Get station id from a sacfile.

    Station id equals to trace id without channel, AU.FORT.10 or G.CAN. this means LOC is none.

    :param file:
    :return:
    """

    tr = read(file, headonly=True)[0]
    header = {}
    sac = tr.stats.get("sac", {})
    header.update(sac)


    net = header["knetwk"]
    stn = header["kstnm"]
    loc = tr.stats.location

    stid = ".".join([net, stn, loc])

    stlo = header["stlo"]
    stla = header["stla"]

    return [stid, stlo, stla]




