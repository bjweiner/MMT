#!python
#
#  Set up a grid of points in RA, Dec, transform it to tangent projection
#  or focal plane coords using routines from desimeter, then either
#  transform it back at another time, to get offsets in RA-Dec space,
#  or compute a second transform from RA, Dec to get offsets in focal
#  plane space
#

# See for ex desimeter/fiberassign.py for some examples of using these
# routines

# desimeter uses a 1/tan approximation for the amount of refraction(elev)
# but I believe that should be close enough (better than 10% fractional) at
# elev > 15 deg

# this or similar environment ?
# source activate spectro2 
# source activate astroconda

import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import (EarthLocation, SkyCoord, ICRS, FK5, AltAz,
                                 Angle, Longitude)
import astropy.units as u
from astropy.time import Time

# hack to let it find desimeter
import sys
sys.path.insert(0, '/Users/bjw/desi/desihub/desimeter/py')
# sys.path.insert(0, '/Users/bjw/desi/desihub')
# these didn't work due to dependencies
#  from desimeter.transform.tan2fp.raytracefit import tan2fp, fp2tan
#      File "/Users/bjw/desi/desihub/desimeter/py/desimeter/transform/tan2fp/raytracefit.py", line 5, in <module>
#      from importlib.resources import files as resource_files
#      ImportError: cannot import name 'files' from 'importlib.resources' (/Users/bjw/software/miniconda3/envs/specdev/lib/python3.8/importlib/resources.py)
# it's complaining about trying to import 'files' from importlib.resources
# perhaps some other thing under desihub ?  But importlib
# should be a python built-in, importlib.resources must be
# in a newer version. It did not work in env astroconda w/ python 3.5
# where importlib doesn't have resources.py, did find
# importlib.resources in env specdev w/python 3.8 but didn't find 'files'.
# look in ../desihub/desimeter/py/desimeter/transform/tan2fp/raytracefit.py
# it's using this in get_raytracefit() and is a real dependency,
# that tan2fp, fp2tan need to work, but it's basically a glorified
# way to get a path to the data/raytrace-tan2fp.json file.

# I don't think I'm actually using tan2fp or fp2tan yet.

from desimeter.transform.radec2tan import radec2tan, tan2radec, hadec2xy
from desimeter.trig import sincosd
# from desimeter.transform.tan2fp.raytracefit import tan2fp, fp2tan

# import these or copy them ?
from area_refraction import construct_grid, get_uv_from_separ, get_times, get_coords
from area_refraction import parangle
from area_refraction import arrow_plot, plot_exposure_offsets

degtorad = np.pi / 180.0
radtodeg = 180.0 / np.pi
VERBOSE = 1
plotnumber = 1

hopkins = EarthLocation.of_site('mmt')
kittpeak = EarthLocation.of_site('kpno')

pressure_kpno = 80. * 1000. * u.Pa
wl = 0.7 * u.micron
maxradius = 1.6

# Return a float of mjd in days and an LST in degrees
# for a location (within the AltAz object) and time.
# the altaz_frame already has the obstime
def time_to_mjd_lst(altaz_frame):
    lon = altaz_frame.location.lon
    mjd = (altaz_frame.obstime).to_value('mjd')
    lmst = altaz_frame.obstime.sidereal_time('mean', longitude=lon)
    return mjd, lmst.deg

# Somewhat like calc_xy_grid in area_refraction.py
# but using desimeter radec2tan, returns tangent plane coordinates
# which are (from the radec2tan comments):
#    x_tan = sin(theta)*cos(phi)
#    y_tan = sin(theta)*sin(phi)
#    where theta,phi are polar coordinates.
#    theta=0 for along the telescope pointing.
#    phi=0 for a star with the same Dec as the telescope
#    pointing but a larger HA (or smaller RA).
# So the x and y are the components of the radial distance from
# the tel axis, and for small theta away from the telescope axis,
# the distance is approx the angle, since theta ~ sin(theta).
# So x and y ought to be in radian units;
# The sign of x should be the same as HA and the opposite of RA,
# so +x means smaller RA; this will require attention
def calc_tan_grid(index_cen, field_grid, altaz_frame):
    ra_tel_deg = field_grid[index_cen].ra.deg
    dec_tel_deg = field_grid[index_cen].dec.deg
    mjd, lst_deg = time_to_mjd_lst(altaz_frame)
    # hourang, parang = parangle(field_grid[index_cen], altaz_frame)
    # I don't know how to treat the hexapod rotation - does it totally
    # follow the parangle ? Since this is an equatorial telescope, I think it
    # should be moving only by a small amount. Connie said they have some
    # pre-computed rate at which to rotate it, to compensate for
    # other misalignments possibly
    hex_ang = 0.0
    ra_deg = field_grid.ra.deg
    dec_deg = field_grid.dec.deg
    # Inputs to radec2tan should be in degrees
    x_tan, y_tan = radec2tan(ra_deg, dec_deg, ra_tel_deg, dec_tel_deg, mjd, lst_deg, hex_ang)
    return x_tan, y_tan

# take a grid of tan angles back to ra, dec
# You only need to pass in the center SkyCoord for the tel ra, dec, not the whole grid
def calc_radec_from_tan(x_tan, y_tan, fieldcen, altaz_frame):
    ra_tel = fieldcen.ra.deg
    dec_tel = fieldcen.dec.deg
    mjd, lst_deg = time_to_mjd_lst(altaz_frame)
    hex_ang = 0.0
    ra_grid, dec_grid = tan2radec(x_tan, y_tan, ra_tel, dec_tel, mjd, lst_deg, hex_ang)
    return ra_grid, dec_grid
    
# This is like calc_exposure_offsets() in area_refraction.py
# To work like area_refraction.py, this is supposed to return two sets
# at times t1, t2,  of x,y offsets in arcsec from field center.

def calc_offsets_desimeter(fieldcen, tstart, tend):
    if VERBOSE > 0:
        print("Using pressure, wl: ", pressure_kpno, wl)
    altaz_frame_kpno_start = AltAz(obstime = tstart, location = kittpeak, pressure = pressure_kpno, obswl = wl)
    altaz_frame_kpno_end   = AltAz(obstime = tend, location = kittpeak, pressure = pressure_kpno, obswl = wl)
    # construct a grid of SkyCoord
    field_grid, nra, ndec = construct_grid(fieldcen.ra, fieldcen.dec, radius=90.0, step=15.0)
    index_cen = int(nra/2) * ndec + int(ndec/2)
    # index_cen = int(nra * ndec /2)
    # index_cen1, xsep1, ysep1, az_start, alt_start, sinth1, costh1 = calc_xy_grid(field_grid, altaz_frame_kpno_start)
    # I think these are same thing
    index_cen1 = index_cen
    xtan1, ytan1 = calc_tan_grid(index_cen, field_grid, altaz_frame_kpno_start)
    xtan2, ytan2 = calc_tan_grid(index_cen, field_grid, altaz_frame_kpno_end)
    
    # here we have to compute the separation, but maybe
    # xsep1 = xtan1 and so on, so try that
    # But the xtan1 and so on aren't like astropy Angles, they are numbers in radians
    # will save some hassle and allow plot_exposure_offsets to work
    # without changes, if we convert the sep's to Angles
    # also the sign of xtan is opposite the sign of RA, so try negating
    # it here
    xsep1 = - Angle(xtan1, u.radian)
    ysep1 = Angle(ytan1, u.radian)
    xsep2 = - Angle(xtan2, u.radian)
    ysep2 = Angle(ytan2, u.radian)
    if VERBOSE > 0:
        print('tan grid first element: ',xtan1[0], ytan1[0], xtan2[0], ytan2[0])
    
    # compute sep_rms - same as in area_refraction.py
    xsepdiff = xsep2 - xsep1
    ysepdiff = ysep2 - ysep1
    if VERBOSE > 1:
        print ("x y grid")
        for i in range(nra*ndec):
            print(i, xsep1[i], ysep1[i], xsepdiff[i], ysepdiff[i])
    # Since these weren't astropy angles, needed to convert to arcsec
    # sepdiff_arcsec_sq = (xsepdiff**2 + ysepdiff**2) / (radtodeg * 3600.)**2
    sepdiff_arcsec_sq = xsepdiff.arcsec**2 + ysepdiff.arcsec**2
    sepdiff_arcsec = np.sqrt( sepdiff_arcsec_sq )
    # Limit the stats to points inside DESI 1.6 deg radius
    # posradius is no longer necessary, I was using it as a check
    # posradius = np.sqrt((xsep1**2 + ysep1**2)*radtodeg**2)
    # inradius = np.where(posradius < maxradius)
    # if using astropy angles , can do this:
    posradius = np.sqrt((xsep1**2 + ysep1**2))
    inradius = np.where(xsep1.deg**2 + ysep1.deg**2 < maxradius**2)
    if VERBOSE > 0:
        print(maxradius, 'tan radius in deg?', posradius)
    sepdiff_sec_rms = np.sqrt(np.mean(sepdiff_arcsec_sq[inradius] ))
    sepdiff_sec_median = np.median( sepdiff_arcsec[inradius] )
    print("RMS offset from tstart - tend, arcsec: {:8.4f}".format(sepdiff_sec_rms))
    print("median offset from tstart - tend, arcsec: {:8.4f}".format(sepdiff_sec_median))
    sep_rms = sepdiff_sec_rms
    
    # get az and alt.
    field_cen_altaz_p_start = field_grid[index_cen].transform_to(altaz_frame_kpno_start)
    field_cen_altaz_p_end   = field_grid[index_cen].transform_to(altaz_frame_kpno_end)
    az_start  = field_cen_altaz_p_start.az
    alt_start = field_cen_altaz_p_start.alt
    az  = az_start
    alt = alt_start
    dec = fieldcen.dec
    
    # compute angle_list - this is (sintheta and costheta of the rot_ang
    #   where rot_ang = -parang) at start, and at end
    # So we need to calc the parang.  This is just used for plotting the
    # up-direction
    hourang1, parang1 = parangle(field_cen_altaz_p_start, altaz_frame_kpno_start)
    hourang2, parang2 = parangle(field_cen_altaz_p_end, altaz_frame_kpno_end)
    angle_list = [ np.sin(-parang1), np.cos(-parang1), np.sin(-parang2), np.cos(-parang2) ]
    
    return nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, hourang1, dec, angle_list, sep_rms


##########

'''
# need an equivalent of the plot offsets
# The plotting routine was pretty involved so these are copied
# from area_refraction.py, and if the calc routine returns the
# same arguments, we can re-use the plot routines.

# Draw arrows at locations using annotate. Here the middle of the arrow
# is at the x,y point. The length is u,v * scale
# So if x,y plot axes are in deg, and u,v are in arcsec, and scale=0.5, then
# a 1 arcsec shift is plotted as an arrow 0.5 axis units (0.5 deg) long.
# Suppress arrrows outside some x-y radius, eg DESI field is circular 1.6 deg.
def arrow_plot(x, y, u, v, scale=1, plotradius=1.6, color='black'):
    npts = len(x)
    for i in range(npts):
        if (x[i]**2 + y[i]**2 < plotradius**2):
            dx = u[i] * scale / 2.0
            dy = v[i] * scale / 2.0
            plt.annotate('', xytext=(x[i]-dx,y[i]-dy), xy=(x[i]+dx,y[i]+dy), arrowprops=dict(color=color, width=0.5, headlength=5, headwidth=3) )
    return

def plot_exposure_offsets(nra, ndec, index_cen, xsep1, ysep1, xsep2, ysep2, az, alt, hourang, dec, exptime, angle_list, sep_rms):
    global plotnumber
    xsepdiff = xsep2 - xsep1
    ysepdiff = ysep2 - ysep1
    # Make an arrow scaled to sepdiff somehow and plot at the location of xsep1, ysep1
    # Cast the angles to deg or arcsec units

    ylimscale = 1.2
    xlimscale = ylimscale * 1.2
    # You can set x and y lims equal to get a square plot, or x larger to get
    # more room for annotations. set_aspect('equal') should keep the circle circular
    xlim1 = xlimscale * np.min(xsep1.deg)
    xlim2 = xlimscale * np.max(xsep1.deg)
    ylim1 = ylimscale * np.min(ysep1.deg)
    ylim2 = ylimscale * np.max(ysep1.deg)
    plt.close()
    plt.clf()
    fig, ax = plt.subplots(1,1)
    plt.axis([xlim2, xlim1, ylim1, ylim2])
    plt.xlabel("Field pos, x (deg arc RA)")
    plt.ylabel("Field pos, y (deg Dec)")
    # ax.set_aspect('equal', 'box')
    ax.set_aspect('equal')
    # toplabel = str(altaz_frame.az) + str(altaz_frame.alt) + str(exptime)
    toplabel = 'az alt {:6.1f} {:6.1f} , exptime {:5.1f} min'.format(az, alt, exptime)
    plt.text(0, ylim2*1.03, toplabel, horizontalalignment='center')
    hadeclabel = 'HA {:5.2f} hr, Dec {:6.1f}'.format(hourang.hour, dec.deg)
    xhadec_label = xlim1 * 0.92
    yhadec_label = ylim2 * 0.88
    plt.text(xhadec_label, yhadec_label, hadeclabel, horizontalalignment='right')

    # plt.scatter is useful to change marker sizes, but I also want to change
    # orientation
    # exaggerate by making 1 arcsec diff be a line about 6 arcmin long?
    # nupscale = 360.0
    # npts = len(xsep1)
    # There is probably a better way to do this
    # for i in range(npts):
        # draw a line from xsep1,ysep1 to xsep2,ysep2 upscaling the length
    # Use annotate to draw an arrow plot where I have control over the scale,
    # and a scale arrow. ascale is a scale to make the arrows fit, so multiply
    # both the arrows and the label by ascale.
    ascale = 0.6
    xscale_label = xlim2 * 0.9
    # Depending whether you print label left or right of scale arrow
    # xscale_label = xlim2 * 0.4
    yscale_label = ylim1 * 0.95
    arrow_plot(xsep1.deg, ysep1.deg, xsepdiff.arcsec, ysepdiff.arcsec, scale=ascale, plotradius=maxradius, color='blue')
    # the 0.3 here is 0.3 arcsec, it goes here and in the label text.
    arrowlength = 0.3*ascale
    plt.annotate('', xytext=(xscale_label,yscale_label), xy=(xscale_label-arrowlength,yscale_label), arrowprops=dict(color='black', width=0.5, headlength=5, headwidth=3) )
    # plt.text(xscale_label+0.02,yscale_label,'0.3" apparent motion',horizontalalignment='right',verticalalignment='center')
    plt.text(xscale_label-arrowlength-0.04,yscale_label,'0.3" apparent motion',horizontalalignment='left',verticalalignment='center')

    # sinth and costh are the angle between the v-axis (up) and the y-axis (Dec)
    # I want to draw an arrow in the v direction, meaning du=0 and dv=length
    # use the transform of 0,v to x,y as in derotate_grid.
    # Previously the arrow direction needed to have x flipped, prob due to
    # mix-up because u (az) increases to right and x (dec) increases to left;
    # after I fixed that in derotate_grid, the up-arrow seems to be ok.
    sinth1, costh1, sinth2, costh2 = angle_list
    vlength = 0.25
    xuparrow1 = -vlength * sinth1
    yuparrow1 = vlength * costh1
    xuparrow2 = -vlength * sinth2
    yuparrow2 = vlength * costh2
    xup_label = xlim1 * 0.8
    yup_label = ylim1 * 0.8
    plt.annotate('', xytext=(xup_label, yup_label), xy=(xup_label+xuparrow1, yup_label+yuparrow1), arrowprops=dict(color='blue', width=0.5, headlength=5, headwidth=3) )
    plt.annotate('', xytext=(xup_label, yup_label), xy=(xup_label+xuparrow2, yup_label+yuparrow2), arrowprops=dict(color='red', width=0.25, headlength=5, headwidth=3) )
    plt.text(xup_label, yup_label, 'up', horizontalalignment='center', verticalalignment='center')
    xrms_label = xlim2 * 0.92
    yrms_label = ylim2 * 0.88
    plt.text(xrms_label, yrms_label, 'rms {:5.2f}'.format(sep_rms), horizontalalignment='left')

    # I think this can be done with pyplot.quiver
    # could use scale keyword, I want to turn off autoscale to make plots comparable,
    # but having trouble controlling the scale / units
    # putting x,y in deg, u,v (=xsep,ysep) in arcsec, setting scale=0.1 and scale_units='x'
    # should mean that a sep of 1 arcsec has the length of 0.1 deg on the x-axis.
    # quiverkey should plot an arrow to show the scale length but not sure how or units yet
    #q = ax.quiver(xsep1.deg, ysep1.deg, xsepdiff.arcsec, ysepdiff.arcsec, pivot='middle', scale_units='x', scale=0.1)
    # ax.quiverkey(q, X=0.3, Y=1.1, U=1, label='Arrow key, length=1', labelpos='E')
    plotname = 'refrac_map_{:04d}.pdf'.format(plotnumber)
    plt.savefig(plotname)
    plotnumber = plotnumber + 1
    plt.show(block=False)
    
    return

'''

def print_plot_arguments(nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, hourang, dec, exptime, angle_list, sep_rms):
    print('nra, ndec, index_cen1: ', nra, ndec, index_cen1)
    print('xsep1 array ', xsep1)
    print('ysep1 array ', ysep1)
    print('xsep2 array ', xsep2)
    print('ysep2 array ', ysep2)
    print('az, alt, hourang, dec, exptime: ', az, alt, hourang, dec, exptime)
    print('angle_list : ', angle_list)
    print('sep_rms : ', sep_rms)
    

##########

# main uses calc_exposure_offsets and plot_exposure_offsets
# update for desimeter-based routines
# if we return the exact same arguments, can use the same plot_exposure_offsets


def main():
    while True:
        fieldcen = get_coords()
        if fieldcen == -1:
            break
        tstart, tend, exptime = get_times()
        # plt.close()
        # nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, angle_list, sep_rms = calc_exposure_offsets(fieldcen, tstart, tend)
        nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, hourang, dec, angle_list, sep_rms = calc_offsets_desimeter(fieldcen, tstart, tend)
        if VERBOSE > 1:
            print_plot_arguments(nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, hourang, dec, exptime, angle_list, sep_rms)
        plot_exposure_offsets(nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, hourang, dec, exptime, angle_list, sep_rms)
    print("quitting")
    return
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()



