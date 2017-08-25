#!/usr/bin/env python

"""
Utilities for Seismic Daylight Imaging

"""

import os
import util
from obspy import read
import gmt as gmt5
import numpy as np

def _read_topo(file):

    with open(file, "r") as fp:
        lst = fp.readlines()

    line = lst[0]
    row = line.split()
    return float(row[-1])

def write_discontinuity_into_header(sacfiles, datapath):

    """



    :param sacfiles:
    :param datapath:
    :return:


    Notes
    -----
    #. Moho depth saved in the header of user1, LAB in the header of user3, MLD in the header user5. if available from other sources


    """
    for file in sacfiles:

        coor = util.get_stid_sac(file=file)
        stid = coor[0]
        print stid

        mohofile = datapath + "/" + stid + ".moho"
        labfile = datapath + "/" + stid + ".lab"

        tr = read(file)[0]

        if os.path.exists(mohofile):
            moho_depth = _read_topo(mohofile)
            tr.stats.sac.user1 = moho_depth

        if os.path.exists(labfile):
            lab_depth = _read_topo(labfile)
            tr.stats.sac.user3 = lab_depth

        tr.write(filename=file, format="SAC")


def plot_reflectivity_velocity(reflectivity_files, datapath, figpath="figs"):

    try:
        os.makedirs(figpath)
    except:
        pass


    for file in reflectivity_files:

        tr = read(file)[0]
        trid = tr.id
        stid = util.stid_from_trid(trid=tr.id)

        header = {}
        sac = tr.stats.get("sac", {})
        header.update(sac)
        depmoho = header["user1"]
        deplab = header["user3"]

        acfile = file
        agcfile = file.replace(".sac.d", ".agc.d")
        ifnfile = file.replace(".sac.d", ".ifn.d")

        vpfile = datapath + "/" + stid + ".vp"
        vsfile = datapath + "/" + stid + ".vs"

        # print vpfile
        # print vsfile

        if not os.path.exists(vpfile):
            print vpfile, "not found. Exit."
            os._exit(0)

        if not os.path.exists(vsfile):
            print vsfile, "not found. Exit."
            os._exit(0)

        psfile = figpath + "/acvel_" + trid + ".ps"
        stnm = tr.stats.sac.kstnm
        print psfile

        gmt = gmt5.Gmt()

        gmt.comment("base")
        gmt.cmd("psbasemap", "-JX5c/-8c -R-1/4/0/200 -Bya50f10+l'Depth [km]' -BWsNe+t'%s' -K  -Bxcintfile+l'P Reflectivity' > %s" % (stnm, psfile))

        proc_sac(file=acfile, xyfile="ac.xy", idn=0.0, scale=1.0)
        gmt.cmd("psxy", "ac.xy -J -R -K -O -W1p >> %s" % (psfile))

        proc_sac(file=agcfile, xyfile="agc.xy", idn=1.5, scale=0.5)
        gmt.cmd("psxy", "agc.xy -J -R -K -O -W1p >> %s" % (psfile))

        proc_sac(file=ifnfile, xyfile="ifn.xy", idn=3, scale=1.0)
        gmt.cmd("psxy", "ifn.xy -J -R -K -O -W1p >> %s" % (psfile))

        # plot moho and depth
        gmt.shell("echo 5.0 %f | gmt psxy -J -R -K -O -SB3p -Gsandybrown -t50  >> %s" % (depmoho, psfile))
        gmt.shell("echo 5.0 %f | gmt psxy -J -R -K -O -SB3p -Gred -t50  >> %s" % (deplab, psfile))

        try:
            gmt.shell("echo 5.0 %f | gmt psxy -J -R -K -O -SB3p -Gcyan -t50  >> %s" % (tr.stats.sac.t5, psfile))
        except:
            pass

        # plot velocity
        gmt.cmd("psbasemap", "-JX4c/-8c -R3/9/0/200 -Bxa1+l'Velocity [km/s]' -Bya50f10+l'Depth [km]' -BwsNe+t'%s' -K -O -X5c >> %s" % (stnm, psfile))

        read_vel(file=vpfile, xyfile="vp.xy")
        gmt.cmd("psxy", "vp.xy -J -R -K -O -W1p >> %s" % (psfile))
        gmt.cmd("psxy", "ak135_vp.xy -J -R -K -O -W1p,- >> %s" % (psfile))

        read_vel(file=vsfile, xyfile="vs.xy")
        gmt.cmd("psxy", "vs.xy -J -R -K -O -W1p >> %s" % (psfile))
        gmt.cmd("psxy", "ak135_vs.xy -J -R -K -O -W1p,- >> %s" % (psfile))

        gmt.shell("echo 4.5 20 Vs | gmt pstext -J -R -K -O >> %s" % (psfile))
        gmt.shell("echo 8.0 20 Vp | gmt pstext -J -R -K -O >> %s" % (psfile))

        # plot moho and depth
        gmt.shell("echo 10.0 %f | gmt psxy -J -R -K -O -SB3p -Gsandybrown -t50  >> %s" % (depmoho, psfile))
        gmt.shell("echo 10.0 %f | gmt psxy -J -R -K -O -SB3p -Gred  -t50  >> %s" % (deplab, psfile))
        try:
            gmt.shell("echo 10.0 %f | gmt psxy -J -R -K -O -SB3p -Gcyan -t50  >> %s" % (tr.stats.sac.t5, psfile))
        except:
            pass

        gmt.comment("end")
        gmt.cmd("psxy", "-J -R -O -T >> %s" % psfile)
        gmt.cmd("psconvert", "-A -P -Tj %s" % psfile)
        gmt.cmd("psconvert", "-A -P -Tf %s" % psfile)

        gmt.execute()

    pass

def proc_sac(file, xyfile, idn, scale):

    tr = read(file)[0]
    z = np.arange(tr.stats.npts) * tr.stats.delta
    d = tr.data
    d = d/max(abs(d))*scale + idn

    fp = open(xyfile, "w")
    for i in range(len(z)):
        fp.write("%f\t%f\n" % (d[i], z[i]))
    fp.close()

def read_vel(file, xyfile):

    with open(file, "r") as fp:
        lst = fp.readlines()

    fp = open(xyfile, "w")
    for line in lst:
        row = line.split()
        dep = float(row[3])
        val = float(row[4])
        fp.write("%f\t%f\n" % (val, dep))
    fp.close()

def sac2xy(file, savepath, std=False, scale=10):
    """
    convert sac file to xy for gmt plot
    :param file:
    :return:
    """

    tr =  read(file)[0]

    npts = tr.stats.npts
    delta =  tr.stats.delta

    # tr.data[:201] = 0.0
    # x = tr.data
    # std = np.std(x)
    #
    # print std, np.amax(x)
    # a = np.where(abs(x)>scale*std)
    # x[a] = np.sign(x[a])*scale*std
    #
    # y = x/max(abs(x))
    tr.normalize()
    y = tr.data

    x = np.arange(0, npts)*delta

    print file
    bn = os.path.basename(file) + ".xy"
    fn = "/".join([savepath, bn])

    fp = open(fn, "w")
    for i in range(len(x)):
        # fp.write("%f\t%f\n" % (x[i], y[i]))
        fp.write("%f\t%f\n" % (y[i], x[i]))
    fp.close()
