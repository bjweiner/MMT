#!python
#
# Differential refraction changes over some area around a point on sky
#

#  The Bennett and Saemundsson formulae used here are from the wikipedia
#  article and are for visual so about 500-550 nm.
#  See R.C. Stone 1996 http://adsabs.harvard.edu/abs/1996PASP..108.1051S
#  for a much more detailed treatment.

#  And see J. Meeus, Astronomical Algorithms chap 13-16

#  see https://docs.astropy.org/en/stable/coordinates/example_gallery_plot_obs_planning.html
# for an example of converting RA, Dec to AltAz at a time, etc
# although they don't use refraction

# this or similar environment
# source activate spectro2 

import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import (EarthLocation, SkyCoord, ICRS, AltAz,
                                 Angle, Longitude)
import astropy.units as u
from astropy.time import Time

degtorad = np.pi / 180.0
VERBOSE = 1
plotnumber = 1

# Bennett formula for refrac as func of apparent altitude, returns in radians
# with a in deg: R in min = cot(a + 7.31/(a + 4.4))
def refrac_apparent(elev_rad):
    degtorad = np.pi / 180.0
    elev_deg = elev_rad / degtorad
    r_arcmin = 1.0 / np.tan( (elev_deg + 7.31 / (elev_deg + 4.4)) * degtorad)
    r_rad = r_arcmin / 60.0 * degtorad
    return r_rad

# Saemundsson formula for refrac as func of true altitude, returns in radians
# with a = atrue in deg: R in min = 1.02 * cot(a + 10.3/(a + 5.11))
def refrac_true(elev_rad):
    degtorad = np.pi / 180.0
    elev_deg = elev_rad / degtorad
    r_arcmin = 1.02 / np.tan( (elev_deg + 10.3 / (elev_deg + 5.11)) * degtorad)
    r_rad = r_arcmin / 60.0 * np.pi / 180.0
    return r_rad

##### Stuff that belonged to the old plotting refrac(elevation) script

# Default pressure = 101.0 kPa, temp = 10 deg C
# 1 atm = 101.3 kPa.  At 6500-8000 ft, expect ~80-75 kPa.
pressure = 101.0
# pressure = 101.0 * 1000. * u.Pa
temp = 10.0

# refrac as function of elevation

elevation = np.arange(1, 90, 0.2) 
elevrad = elevation * degtorad
refrac_rad = refrac_true(elevrad)

refrac = (pressure / 101.0) * 283.0 / (273.0 + temp) * refrac_rad
refrac_sec = refrac / degtorad  * 3600.0

#####

# old
fieldsize_deg = 2.0
fieldsize_min = fieldsize_deg * 60.0
fieldsize_rad = fieldsize_deg * degtorad

# -----

# Set up grid around the center - in RA, Dec 
# I cannot figure out how to add offset numbers to a SkyCoord to make another SkyCoord
#  later: you could do this with the method .spherical_offsets_by() 
# but you can do math with Angles
def construct_grid(fieldcen_ra, fieldcen_dec, radius=60.0, step=10.0):
    # radius = 60.0
    # step = 10.0
    # Larger step to make a smaller array for printing/testijng
    # step = 30.0
    ra_offset_vec  = Angle(np.arange(-radius, radius+0.01, step) * u.arcmin)
    dec_offset_vec = Angle(np.arange(-radius, radius+0.01, step) * u.arcmin)
    nra  = len(ra_offset_vec)
    ndec = len(dec_offset_vec)

    # with np.tile and np.repeat, make long vectors of offsets representing a grid, 
    # where ra_offset is -60 n times, then -55 n times, and so on, while dec_offset 
    # changes faster (like a least significant bit)
    # although these are Angles, I don't think they can automatically know the cos(dec)
    # factor, so I need to convert arcmin to min of time (I will then take it out
    # again later)
    dec_offset_array = np.tile(dec_offset_vec, nra)
    field_grid_dec = fieldcen_dec + dec_offset_array

    ra_offset_array  = np.repeat(ra_offset_vec, ndec) / np.cos(field_grid_dec.radian)
    field_grid_ra  = fieldcen_ra + ra_offset_array 

    # Grid of coords around the center in RA, Dec
    field_grid = SkyCoord(field_grid_ra, field_grid_dec, frame='icrs')
    return field_grid, nra, ndec

# use the method .spherical_offsets_by(), passing in a SkyCoord rather than
# ra and dec of the center. Can you pass it vectors? But this method only
# exists in astropy > v x.x
def construct_grid_offsets(fieldcen, radius=60.0, step=10.0):
    ra_offset_vec  = Angle(np.arange(-radius, radius+0.01, step) * u.arcmin)
    dec_offset_vec = Angle(np.arange(-radius, radius+0.01, step) * u.arcmin)
    nra  = len(ra_offset_vec)
    ndec = len(dec_offset_vec)
    dec_offset_array = np.tile(dec_offset_vec, nra)
    # I think here I'm passing in an offset angle and don't need to divide by cos(dec)
    ra_offset_array  = np.repeat(ra_offset_vec, ndec)

    field_grid = fieldcen.spherical_offsets_by(ra_offset_array, dec_offset_array)
    
    return field_grid, nra, ndec

# u and v are offsets from center in arc angles, so convert from az back to arc units
# a circle of azimuth is shorter than arc units by cos(alt).
# This is just using az and alt and not doing a reprojection of spherical coordinates
# onto the center point (like a tan reprojection), which could be significant.
# Both changes in az to arc conversion and orthogonality of axes likely matter.
# However, I did think this would be closer than it tests ...
def get_uv_from_altaz(grid_altaz, index_cen):
    # u = (grid_altaz.az - grid_altaz[index_cen].az) * np.cos(grid_altaz.alt * degtorad)
    u1 = (grid_altaz.az - grid_altaz[index_cen].az) * np.cos(grid_altaz.alt.radian)
    v1 = grid_altaz.alt - grid_altaz[index_cen].alt
    return u1, v1

# Meeus's formulae for separation and PA: angular sep d
#  cos d = sin dec1 sin dec2 + cos dec1 cos dec2 cos (ra1-ra2)    Meeus eq 17.1
# he says it doesn't work near d = 0 or 180 deg, because of taking the arccos
# position angle P of point 1 measured from point 2 is:
#  tan P = sin(ra1 - ra2) / (cos dec2 tan dec1 - sin dec2 cos(ra1-ra2) )   Meeus p.116
# This has not been giving me sensible results.

def get_uv_from_meeus(grid_altaz, index_cen):
    dec1 = grid_altaz.alt
    dec2 = grid_altaz[index_cen].alt
    ra1  = grid_altaz.az
    ra2  = grid_altaz[index_cen].az
    radiff = ra1 - ra2
    cos_sep = np.sin(dec1.radian) * np.sin(dec2.radian) + np.cos(dec1.radian) * np.cos(dec2.radian) * np.cos(radiff.radian)
    sep = np.arccos(cos_sep)  # sep is in radians?
    tan_pa = np.sin(radiff.radian) / (np.cos(dec2.radian) * np.tan(dec1.radian) - np.sin(dec2.radian) * np.cos(radiff.radian))
    pa  = np.arctan(tan_pa)   # or use atan2 ?
    # PA is measured from the v axis through -u, so need to use PA+90 to find u and v
    # This will not work at sep = 0.0 because PA will blow up?
    if (sep > 1.0e-6):
        u_rad = sep * np.cos(pa + np.pi/4.0)
        v_rad = sep * np.sin(pa + np.pi/4.0)
    else:
        u_rad = 0.0
        v_rad = 0.0
    # make into Angles
    u_ang = Angle(u_rad, u.radian)
    v_ang = Angle(v_rad, u.radian)
    return u_ang, v_ang

# try using astropy as if altaz were ra dec?

# I am not sure if you can pass a vector of coords to the spherical offsets method
# but presumably you can find the offsets from the grid to center, and then
# change the sign.  What I don't know is if this method will actually give the
# offsets in apparent angle, which is what I need. It may outsmart me by
# using the true coordinates that belong to each point.
# Testing seems to show that if you call spherical_offsets_to() on AltAz objects
# it returns the apparent offsets, which I think are in arc units,
# but for spherical_offsets_to() on icrs objects, it returns the true offsets.
# This is sensible behavior.

# Note this original version is computing the offset from the grid point
# to the center and then flipping the sign. I tested and this is _not_
# always equal to computing the offset from center to grid point, presumably
# due to different center of projection of coordinates (offsetting along
# different great-circles)
def get_uv_from_separ_old(grid_altaz, index_cen):
    center_pt = grid_altaz[index_cen]
    diff1, diff2 = grid_altaz.spherical_offsets_to(center_pt)
    u1 = -diff1
    v1 = -diff2
    return u1, v1

# I wasn't sure if you could pass a vector of SkyCoord as argument
# to spherical_offsets_to().  I think you can; vs trying to make
# the vector by a for loop appending to an array was giving me
# hell with astropy unit errors because of the angle units.
# I think casting the angles to radians with .radian might solve that.
# This version computes the offset from center to each grid point.
def get_uv_from_separ(grid_altaz, index_cen):
    center_pt = grid_altaz[index_cen]
    u1, v1 = center_pt.spherical_offsets_to(grid_altaz)
    return u1, v1

# Or maybe I could use astropy skyoffset_frame() for something here.

# old
def set_observatory():
    hopkins = EarthLocation.of_site('mmt')
    kpno = EarthLocation.of_site('kpno')
    # Old numbers for ephem package
    # hopkins.lon = '-110.883'
    # hopkins.lat = '31.688'
    # hopkins.elevation = 2608
    # bok.lon = '-111:36:01.6'
    # bok.lat = '31:57:46.5'
    # bok.elevation = 2071
    return kpno

# -----
# Set the location

hopkins = EarthLocation.of_site('mmt')
kittpeak = EarthLocation.of_site('kpno')

# Test fields and times. Use a field that is fairly low at transit
# Dec -25 transits at about airmass 1.85, altitude 32.9 deg
# If precession is a concern, maybe try two diff fields since the motion
# of the north pole is toward about RA ~ 3 ??

testra = "13:00:00"
testdec = "-25:00:00"
fieldcen_ra = Angle("13h00m00s")
fieldcen_dec = Angle("-25d00m00s")
fieldcen = SkyCoord(fieldcen_ra, fieldcen_dec, frame='icrs')
# fieldcen = SkyCoord(testra, testdec, frame='icrs', unit=(u.hourangle, u.deg))
# For this field and KPNO, I used skycalc to find when it transits
# on ~ 2000 Apr 2 and 2025 Apr 2, in UTC
time1str = "2000-04-02 07:43:00"
time2str = "2025-04-02 07:44:00"
time2000 = Time(time1str)
time2025 = Time(time2str)

time_obs = time2025

# If not hardcoding, input some UT and some RA, Dec

def get_times():
    time_str = input("Time string for exposure start, eg: 2025-04-02 07:30:00: ")
    if time_str == '':
        print("using 2025-04-02 07:30:00")
        time_str = "2025-04-02 07:30:00"
    try:
        time_start = Time(time_str)
    except Exception as e:
        print("Error, using 2025-04-02 07:30:00")
        time_str = "2025-04-02 07:30:00"
        time_start = Time(time_str)
    length_str = input("Exposure length in min: ")
    explength = float(length_str)
    time_end = time_start + explength * u.min
    return time_start, time_end, explength

def get_coords():
    # coord_str = input("Field ra dec, eg: 13h00m00s +10d00m00s : ")
    coord_str = input("Field ra dec, eg: 13h00m00s +10d00m00s [-1 -1 to quit]: ")
    coords = coord_str.split()
    if len(coords) < 2:
        return -1
    rastr = coords[0]
    decstr = coords[1]
    if rastr == "-1" or decstr == "-1":
        return -1
    else:
        try:
            fieldcen = SkyCoord(Angle(rastr), Angle(decstr), frame="icrs")
        except:
            print("Error, using 13h00m00s -10d00m00s")
            rastr = "13h00m00s"
            decstr = "-10d00m00s"
            fieldcen = SkyCoord(Angle(rastr), Angle(decstr), frame="icrs")
        return fieldcen
            
        
# Calculate hour angle, (apparent?) azimuth, elevation

# if you don't specify a pressure or set pressure = 0, you should get
# no refraction.  transform_to() was rather slow - maybe just the first time
# if it downloads an iers file
pressure_kpno = 80. * 1000. * u.Pa
wl = 0.7 * u.micron
maxradius = 1.6

altaz_frame_p0 = AltAz(obstime = time_obs, location = kittpeak)
altaz_frame_pkpno = AltAz(obstime = time_obs, location = kittpeak, pressure = pressure_kpno, obswl = wl)
altaz_frame_pkpno2000 = AltAz(obstime = time2000, location = kittpeak, pressure = pressure_kpno, obswl = wl)
fieldcen_altaz_p0 = fieldcen.transform_to(altaz_frame_p0)
fieldcen_altaz_pkpno = fieldcen.transform_to(altaz_frame_pkpno)

# print(fieldcen_altaz_p0)
# print(fieldcen_altaz_pkpno)
# The important parts are: no refrac az, alt (179.96314566, 32.89824134)
# with refrac az, alt (179.96314566, 32.91837616), it defaulted to 0.0 deg C and obswl=1 micron
# refrac makes a difference of 0.020135 deg = 72.5 arcsec, in Alt, no diff in Az
# with obswl = 0.7 micron, az, alt (179.96314566, 32.91849711), which is 0.44 arcsec
# different from 1 micron (as one might expect from the usual size of the atmospheric
# dispersion effect)

def run_initial_tests(ra, dec, altaz_frame_p0, altaz_frame_p):
    altaz_frame_pkpno = altaz_frame_p

    # Does separation returned in AltAz depend on refraction?
    # These altaz frames include an obstime in Apr 2025, and pressure = 0 or 80 kPa
    test1a = SkyCoord('180.0d', '20.0d', frame=altaz_frame_p0)
    test1b = SkyCoord('180.0d', '22.0d', frame=altaz_frame_p0)
    test2a = SkyCoord('180.0d', '20.0d', frame=altaz_frame_pkpno)
    test2b = SkyCoord('180.0d', '22.0d', frame=altaz_frame_pkpno)
    test_u_p0, test_v_p0 = test1a.spherical_offsets_to(test1b)
    test_u_pkpno, test_v_pkpno = test2a.spherical_offsets_to(test2b)
    print("p = 0:", test_u_p0, test_v_p0)
    print("p = kpno:", test_u_pkpno, test_v_pkpno)
    # There is a microscopic difference, which suggests that the sep returned this way 
    # is in apparent coords (az, alt) and not in true coords
    # p = 0: 0d00m00s 2d00m00s
    # p = kpno: 0d00m00s 2d00m00.0004s
    test3a = test1a.transform_to('icrs')
    test3b = test1b.transform_to('icrs')
    test4a = test2a.transform_to('icrs')
    test4b = test2b.transform_to('icrs')
    test_u_p0, test_v_p0 = test3a.spherical_offsets_to(test3b)
    test_u_pkpno, test_v_pkpno = test4a.spherical_offsets_to(test4b)
    print("p = 0:", test_u_p0, test_v_p0)
    print("p = kpno:", test_u_pkpno, test_v_pkpno)
    # When calculated on coords in ICRS, the separations are a little
    # different from 0,2 deg. The RA diff should be due to precession and the
    # Dec diff mostly refraction.
    # p = 0: 0d00m05.9487s 2d00m00.1871s
    # p = kpno: 0d00m05.9624s 2d00m12.7898s
    test5a = SkyCoord('180.0d', '28.5d', frame=altaz_frame_pkpno)
    test5b = SkyCoord('180.0d', '31.5d', frame=altaz_frame_pkpno)
    # test5a = SkyCoord('180.0d', '28.5d', frame=altaz_frame_pkpno2000)
    # test5b = SkyCoord('180.0d', '31.5d', frame=altaz_frame_pkpno2000)
    test6a = test5a.transform_to('icrs')
    test6b = test5b.transform_to('icrs')
    test_u_pkpno, test_v_pkpno = test5a.spherical_offsets_to(test5b)
    print("p = kpno, altaz:", test_u_pkpno, test_v_pkpno)
    test_u_pkpno, test_v_pkpno = test6a.spherical_offsets_to(test6b)
    print("p = kpno, icrs:", test_u_pkpno, test_v_pkpno)
    # Here I get the amount of differential refraction across 3 deg when at
    # elevation 30 deg due south.  I should really use obstime 2000 to eliminate
    # precession; meanwhile I get that the 3 degree field is compressed from
    # ICRS to AltAz by 10 arcsec. That is 0.9 part in 1000. If that's
    # the relevant fraction, and we sit on a field for 15 minutes = 15 min of time
    # theta ~ 3.75 deg ~ 0.065 radian, then can I expect points to move off their
    # locations by about 0.065 / 1000 ? As the direction of field compression
    # changes as the field rotates?  I'm not sure this logic works, but if so,
    # a point at the edge of the field, 5400 arcsec from center, would move by
    # ~ 0.35 arcsec. Which is suspiciously close to the amount of motion that
    # dither tests are seeing?  Later I worked this out more carefully, see
    # Readme.offset_script, and the fractional change in distance from center
    # is about ~ theta * (dy-dx), where theta is the change in par angle and
    # dy and dx are the altitude and azimuth compression, about 880 and 220 ppm.
    # Trying with the 2000 obstime, it nearly eliminates the RA offset so that
    # was precession, and the Dec offset is basically unchanged so that is pure
    # differential refraction and the above conclusions hold.
    # Note that was just about the change in direction of field compression,
    # and not the change in compression due to field rise/set.
    # for the 2025 obstime:
    # p = kpno, altaz: 0d00m00s 3d00m00.0001s
    # p = kpno, icrs: 0d00m07.943s 3d00m10.0668s
    # for the Apr 2000 obstime:
    # p = kpno, altaz: 0d00m00s 3d00m00.0001s
    # p = kpno, icrs: 0d00m00.7478s 3d00m10.0664s

def run_grid_tests(ra, dec, altaz_frame_p0, altaz_frame_p):
    fieldcen_ra = ra
    fieldcen_dec = dec
    altaz_frame_pkpno = altaz_frame_p
    field_grid, nra, ndec = construct_grid(fieldcen_ra, fieldcen_dec)
    # Make a grid with fewer points for simple printing-out etc
    # field_grid, nra, ndec = construct_grid(fieldcen_ra, fieldcen_dec, radius=60.0, step=30.0)

    # Transformed to alt az with pressure=0 / no refrac, and with refrac 
    field_grid_altaz_p0 = field_grid.transform_to(altaz_frame_p0)
    field_grid_altaz_pkpno = field_grid.transform_to(altaz_frame_pkpno)

    # How to determine difference in offset in apparent position over the field
    # between two times?
    
    # the parallactic angle is the sky PA from N to E that points to the zenith, 
    # so the up axis in az-alt.  So if I calc. offsets from center in u and v that are
    # aligned with az-alt, then rotate to axes N and W, I'm de-rotating by the
    # parallactic angle. Imagine if paral angle is 30, then the alt axis is 30 deg
    # ccw from the Dec axis. So if u,v are az,alt, and x,y are RA,Dec, then
    #   u =  x cos theta + y sin theta
    #   v = -x sin theta + y cos theta
    #   x =  u cos theta - v sin theta
    #   y =  u sin theta + v cos theta

    # it seems like SkyCoord can't tell me the parallactic angle and I would
    # need to use astroplan. But I should also be able to get this angle from
    # the angle in az-alt system between points in the grid on a line of const RA.

    # compute coords in focal plane, where u and v are aligned with az-alt but
    # in deg of arc, and x and y are aligned with RA and Dec but in deg of arc
    # and apparent, so not real RA and Dec

    # center of grid, and northernmost point in the vector that ran through
    # the center.  The dec should change faster in the grid than the ra
    ndechalf = int(ndec/2)
    index_cen = int(nra/2) * ndec + int(ndec/2)
    index_north = int(nra/2) * ndec + ndec-1
    index_south = int(nra/2) * ndec
    ngrid = nra * ndec

    # print (ngrid, index_south, index_cen, index_north)
    # 169 78 84 90

    # These routines should be returning apparent separations in u, v which
    # are degrees of arc referenced to the grid center.  Which is not quite the
    # same as separation in az, alt.

    u1, v1 = get_uv_from_altaz(field_grid_altaz_pkpno, index_cen)
    print (u1)
    print (v1)
    u2, v2 = get_uv_from_meeus(field_grid_altaz_pkpno, index_cen)
    u3, v3 = get_uv_from_separ(field_grid_altaz_pkpno, index_cen)

    indexes = [0, index_cen, ngrid-1]
    print(indexes)
    for i in indexes:
        print(i, "u values: ",u1[i], u2[i], u3[i])
        print(i, "v values: ",v1[i], v2[i], v3[i])

    #   0 u values:  0d59m58.7701s 0.0247212rad 0d59m19.2368s
    #   0 v values:  -1d00m32.2828s 4.62185e-05rad -1d00m51.6164s
    #  84 u values:  0d00m00s nanrad 0d00m00s
    #  84 v values:  0d00m00s nanrad -0d00m00s
    # 168 u values:  -0d59m58.7959s 0.0246215rad -1d00m40.4948s
    # 168 v values:  0d59m22.5168s -5.36019e-05rad 0d59m01.1697s

    # The 3 methods do not agree yet. The Meeus values are way off after converting to deg.
    # The altaz values are off from separ by typ ~ 40 sec in u, 20 sec in v. Because
    # I just computed the diff in az and alt coordinates rather than a correct
    # projection, the altaz method is probably just wrong by that amount and I trust
    # separ more.  Just looking at separ, the last column, the 0th cell separ
    # is less+ in u by 41" and more in v by 52" (compared to 1.0 deg).
    # The center cell is at 0,0 by design. The 168th cell is more- in u by 41" and
    # less in v by 58".
    # This may be affected by rotation because the field isn't exactly on the meridian,
    # and not all due to compression / expansion. Better to compare along the center
    # row/col of the grid, center to sides, rather than center to corners; see below.

    # Print columns in grid: left, middle, right col; bottom, middle, top cell of each col
    indexes = [0, 0+ndechalf, ndec-1, index_cen-ndechalf, index_cen, index_cen+ndechalf, ngrid-ndec, ngrid-1-ndechalf, ngrid-1]
    print(indexes)
    # [0, 6, 12, 78, 84, 90, 156, 162, 168]
    for i in indexes:
        print(i, "u values: ",u1[i], u2[i], u3[i])
        print(i, "v values: ",v1[i], v2[i], v3[i])

    #   0 u values:  0d59m58.7701s 0.0247212rad 0d59m19.2368s
    #   0 v values:  -1d00m32.2828s 4.62185e-05rad -1d00m51.6164s
    #   6 u values:  0d59m59.1429s 0.0123895rad 0d59m58.6678s
    #   6 v values:  -0d00m35.3002s -0.0122868rad -0d00m55.6183s
    #  12 u values:  0d59m59.51s -5.84666e-05rad 1d00m41.2085s
    #  12 v values:  0d59m21.8037s 0.0246215rad 0d59m00.4483s
    #  78 u values:  -0d00m00.3613s 0.0123301rad -0d00m00.3574s
    #  78 v values:  -0d59m57.0782s 0.0123326rad -0d59m57.0782s
    #  84 u values:  0d00m00s nanrad 0d00m00s
    #  84 v values:  0d00m00s nanrad -0d00m00s
    #  90 u values:  0d00m00.3555s 0.0123306rad 0d00m00.3597s
    #  90 v values:  0d59m57.228s 0.0123331rad 0d59m57.228s
    # 156 u values:  -0d59m59.4901s 4.12541e-05rad -0d59m19.9563s
    # 156 v values:  -1d00m31.5627s 0.0247212rad -1d00m50.9042s
    # 162 u values:  -0d59m59.14s -0.0122892rad -0d59m58.673s
    # 162 v values:  -0d00m34.5838s 0.012387rad -0d00m54.9019s
    # 168 u values:  -0d59m58.7959s 0.0246215rad -1d00m40.4948s
    # 168 v values:  0d59m22.5168s -5.36019e-05rad 0d59m01.1697s

    # Just looking at the last column, u and v from separ(), use to disentangle
    # compression and rotation etc.  Comparing along the row or column through the
    # center is better than comparing the corners of the grid.
    # Comparing along row i= 6, 84, 162: u is compressed by 1.33 arcsec per 1 deg,
    #   where a compression of about 0.22/1000 = 0.8" is expected
    #   v is pulled down by 55 arcsec on each side wrt center, with ~ 0.3" asymmetry 
    #   over 1 deg. the v pull down is perhaps from the edges of the field being compressed
    #   differently from the center? very un-intuitive. Small asym due to rotation.
    # Comparing u along row i= 0, 78, 156: u is inward by 40-41" per 1 deg in bottom row
    #   v along row i= 0, 78, 156: pull-down of ~54" in v at edges.
    # Comparing u along row i= 12, 90, 168: u is outward by 40-41" per 1 deg in top row
    #   v along row i= 12, 90, 160: pull-down of ~59" in v at edges.
    # along col i= 78, 84, 90: u is off by -/+0.36" per 1 deg at bot/top prob due to rotation
    #   v is compressed by 2.92" over 1 deg at lower edge, by 2.77" over 1 deg at upper.
    #   This v compression of ~0.8/1000 is about what I expect.
    # along col i= 0, 6, 12: as we go up, u expands + outward by about 40" over 1 deg
    #   v is compressed by about 4" over 1 deg
    # along col i= 156, 162, 168: as we go up, u expands - outward by about 40" over 1 deg
    #   v is compressed by about 4" over 1 deg

# Call this routine with an altaz_frame that will already have a time and
# pressure defined.  The field_grid is in ra, dec and comes from construct_grid()
def calc_altaz_grid(field_grid, altaz_frame_p):

    # field_grid, nra, ndec = construct_grid(ra, dec)
    # Make a grid with fewer points for simple printing-out etc
    # field_grid, nra, ndec = construct_grid(ra, dec, radius=60.0, step=30.0)

    # Transformed to alt az with pressure=0 / no refrac, and with refrac 
    field_grid_altaz_p = field_grid.transform_to(altaz_frame_p)
    return field_grid_altaz_p

# Call this routine with an altaz_frame that will already have a time and
# pressure defined.  The field_grid is in ra, dec and comes from construct_grid()
# u, v are separations from the grid center, in apparent, and should be
# in arc units, but along az-alt axes. I think.
def calc_uv_grid(field_grid, altaz_frame_p):

    # Transformed to alt az with pressure=0 / no refrac, and with refrac 
    field_grid_altaz_p = field_grid.transform_to(altaz_frame_p)
    # Assume that the grid is centered on the middle element
    ngrid = len(field_grid)
    index_cen = int(ngrid/2)
    az = field_grid_altaz_p[index_cen].az
    alt = field_grid_altaz_p[index_cen].alt
    print("az, alt: ", az, alt)
    # u_sep, v_sep = get_uv_from_altaz(field_grid_altaz_p, index_cen)
    u_sep, v_sep = get_uv_from_separ(field_grid_altaz_p, index_cen)
    return index_cen, u_sep, v_sep, az, alt

# index_cen is where the center is, n_off is number of steps
# away to use to find the angle between the systems, int(ndec/2) or less.
# Rotate u, v onto offsets x,y that are aligned with the ra-dec axes.
# I'm computing the rotation between apparent u,v and ICRS ra,dec,
# which is the parallactic angle in 2000, not epoch-of-date. It's not
# clear that this matters as long as I rotate start and end onto the
# same system.
# Note that: if you're looking south from eg KPNO at a field on the
# equator, azimuth increases westward (to your right), but RA increases
# eastward (to your left). Alt and Dec are both up.
# If you look north, azimuth and RA both increase to right, but now Alt is
# up and Dec is down. The southward opposition of axes may have been
# causing a sign error in the uv-xy rotation earlier.
# I am not sure if S and N require treating as separate cases. So far,
# I put in the S formulae and looks like it also draws North fields correctly.

# I drew little axis diagrams for uv and xy rotated by theta, and looking south:
#    u = - x costheta - y sintheta
#    v = - x sintheta + y costheta
#    x = - u costheta - v sintheta
#    y = - u sintheta + v costheta
# Looking north and over-the-celestial pole, the directions of x and y are reversed, so
#    u = + x costheta + y sintheta
#    v = + x sintheta - y costheta
#    x = + u costheta + v sintheta
#    y = + u sintheta - v costheta
def derotate_grid(index_cen, n_off, ra_grid, dec_grid, u1, v1):
    index_south = index_cen - n_off
    index_north = index_cen + n_off
    ra1  = ra_grid[index_cen]
    dec1 = dec_grid[index_cen]
    ra2  = ra_grid[index_south]
    dec2 = dec_grid[index_south]
    ra3  = ra_grid[index_north]
    dec3 = dec_grid[index_north]
    utmp1 = u1[index_cen]
    vtmp1 = v1[index_cen]
    utmp2 = u1[index_south]
    vtmp2 = v1[index_south]
    utmp3 = u1[index_north]
    vtmp3 = v1[index_north]
    dra21  = ra2 - ra1
    ddec21 = dec2 - dec1
    dutmp21 = utmp2 - utmp1
    dvtmp21 = vtmp2 - vtmp1
    # if we aligned with the grid colum, the dra = dx term is 0.  
    sintheta21 = -dutmp21 / ddec21
    costheta21 = dvtmp21 / ddec21
    #
    dra31  = ra3 - ra1
    ddec31 = dec3 - dec1
    dutmp31 = utmp3 - utmp1
    dvtmp31 = vtmp3 - vtmp1
    # if the dra term is 0.  Does this make sintheta an Angle or dimensionless?
    sintheta31 = -dutmp31 / ddec31
    costheta31 = dvtmp31 / ddec31
    if VERBOSE > 0:
        print("sintheta, costheta from south grid step: {:9.6f} {:9.6f}".format(sintheta21, costheta21))
        print("sintheta, costheta from north grid step: {:9.6f} {:9.6f}".format(sintheta31, costheta31))
        # print("sintheta, costheta from south grid step: ", sintheta21, costheta21)
        # print("sintheta, costheta from north grid step: ", sintheta31, costheta31)
    sintheta = (sintheta21 + sintheta31) / 2.0
    costheta = (costheta21 + costheta31) / 2.0
    if VERBOSE > 0: 
        print("sintheta, costheta: {:9.6f} {:9.6f}".format(sintheta, costheta))
        # print("sintheta, costheta: ", sintheta, costheta)
    # print(u1)
    # print(v1)
    # I sometimes got astropy unit hell here, it may be necessary to cast some
    # of the angles to a unit, like: sintheta.radian, if sintheta is an Angle, which
    # it really should not be.
    x_sep = -u1 * costheta - v1 * sintheta
    y_sep = -u1 * sintheta + v1 * costheta
    return index_cen, x_sep, y_sep, sintheta, costheta


# Here we calculate the separations in apparent units, and then
# derotate to put them on the ra-dec axes
def calc_xy_grid(field_grid, altaz_frame_p):

    # assumes that the grid is centered on the middle element
    index_cen, u_sep, v_sep, az_cen, alt_cen = calc_uv_grid(field_grid, altaz_frame_p)
    # I do need to know the num of elements of a grid column here (ndec)
    # although I can just use an offset of 1 or 2 element if the calc is precise enough
    # However if there's a lot of field rotation that may not be enough
    # n_off = int(ndec/2)
    n_offmax = int(np.sqrt(len(field_grid)) / 2)
    n_off = np.min([2, n_offmax])
    index_south = index_cen - n_off
    index_north = index_cen + n_off

    index_tmp, x_sep, y_sep, sintheta, costheta = derotate_grid(index_cen, n_off, field_grid.ra, field_grid.dec, u_sep, v_sep)
    return index_cen, x_sep, y_sep, az_cen, alt_cen, sintheta, costheta

def calc_exposure_offsets(fieldcen, tstart, tend):
    # fieldcen is a SkyCoord
    # fieldcen = get_coords()
    # tstart, tend, exptime = get_times()

    if VERBOSE > 0:
        print("Using pressure, wl: ", pressure_kpno, wl)

    altaz_frame_kpno_start = AltAz(obstime = tstart, location = kittpeak, pressure = pressure_kpno, obswl = wl)
    altaz_frame_kpno_end   = AltAz(obstime = tend, location = kittpeak, pressure = pressure_kpno, obswl = wl)
    
    field_grid, nra, ndec = construct_grid(fieldcen.ra, fieldcen.dec, radius=90.0, step=15.0)
    # field_grid, nra, ndec = construct_grid_offsets(fieldcen, radius=60.0, step=30.0)
    
    # The two sets of uv grids should differ by the field rotation from start
    # to end so comparing them directly isn't that useful, although one could
    # imagine fitting for the optimum rotation, or calculating the parallactic
    # angle and de-rotating by that
    # index_cen1, usep1, vsep1 = calc_uv_grid(field_grid, altaz_frame_kpno_start)
    # index_cen1, usep2, vsep2 = calc_uv_grid(field_grid, altaz_frame_kpno_end)
    # print ("u v grid")
    # for i in range(nra*ndec):
    #    print(i, usep1[i], vsep1[i])
    
    # However, the xy grids should be derotated to the declination axis at field
    # center, so comparable to each other
    index_cen1, xsep1, ysep1, az_start, alt_start, sinth1, costh1 = calc_xy_grid(field_grid, altaz_frame_kpno_start)
    index_cen2, xsep2, ysep2, az_end, alt_end, sinth2, costh2 = calc_xy_grid(field_grid, altaz_frame_kpno_end)
    print("Start pos and exptime: {:s} {:s} {:8.2f} min".format(az_start, alt_start, (tend-tstart).sec/60))
    if (alt_start.deg < 10 or alt_end.deg < 10):
        print ('Warning, altitude < 10 deg, results unreliable.',az_start, alt_start, az_end, alt_end)

    xsepdiff = xsep2 - xsep1
    ysepdiff = ysep2 - ysep1
    if VERBOSE > 1:
        print ("x y grid")
        for i in range(nra*ndec):
            print(i, xsep1[i], ysep1[i], xsepdiff[i], ysepdiff[i])
    sepdiff_arcsec_sq = xsepdiff.arcsec**2 + ysepdiff.arcsec**2
    sepdiff_arcsec = np.sqrt( sepdiff_arcsec_sq )
    # Limit the stats to points inside DESI 1.6 deg radius
    inradius = np.where(xsep1.deg**2 + ysep1.deg**2 < maxradius**2)
    sepdiff_sec_rms = np.sqrt(np.mean(sepdiff_arcsec_sq[inradius] ))
    sepdiff_sec_median = np.median( sepdiff_arcsec[inradius] )
    print("RMS offset from tstart - tend, arcsec: {:8.4f}".format(sepdiff_sec_rms))
    print("median offset from tstart - tend, arcsec: {:8.4f}".format(sepdiff_sec_median))

    # return the sin, cos from the start to plot an arrow. 
    # The end might matter too if there is a lot of field rotation.
    angle_list = [sinth1, costh1, sinth2, costh2]
    return nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az_start, alt_start, angle_list, sepdiff_sec_rms

# Draw arrows at locations using annotate. Here the middle of the arrow
# is at the x,y point. The length is u,v * scale
# So if x,y plot axes are in deg, and u,v are in arcsec, and scale=0.5, then
# a 1 arcsec shift is plotted as an arrow 0.5 axis units (0.5 deg) long.
# Suppress arrrows outside some radius, eg DESI field is circular.
def arrow_plot(x, y, u, v, scale=1, plotradius=1.6, color='black'):
    npts = len(x)
    for i in range(npts):
        if (x[i]**2 + y[i]**2 < plotradius**2):
            dx = u[i] * scale / 2.0
            dy = v[i] * scale / 2.0
            plt.annotate('', xytext=(x[i]-dx,y[i]-dy), xy=(x[i]+dx,y[i]+dy), arrowprops=dict(color=color, width=0.5, headlength=5, headwidth=3) )
    return

def plot_exposure_offsets(nra, ndec, index_cen, xsep1, ysep1, xsep2, ysep2, az, alt, exptime, angle_list, sep_rms):
    global plotnumber
    xsepdiff = xsep2 - xsep1
    ysepdiff = ysep2 - ysep1
    # Make an arrow scaled to sepdiff somehow and plot at the location of xsep1, ysep1
    # Cast the angles to deg or arcsec units

    limscale = 1.2
    xlim1 = limscale * np.min(xsep1.deg)
    xlim2 = limscale * np.max(xsep1.deg)
    ylim1 = limscale * np.min(ysep1.deg)
    ylim2 = limscale * np.max(ysep1.deg)
    plt.close()
    plt.clf()
    fig, ax = plt.subplots()
    plt.axis([xlim2, xlim1, ylim1, ylim2])
    plt.xlabel("Field pos, x (deg arc RA)")
    plt.ylabel("Field pos, y (deg Dec)")
#    toplabel = str(altaz_frame.az) + str(altaz_frame.alt) + str(exptime)
    toplabel = 'az alt {:6.1f} {:6.1f} , exptime {:5.1f} min'.format(az, alt, exptime)
    plt.text(0, ylim2*1.03, toplabel, horizontalalignment='center')

    # plt.scatter is useful to change marker sizes, but I also want to change
    # orientation
    # exaggerate by making 1 arcsec diff be a line about 6 arcmin long?
    # nupscale = 360.0
    # npts = len(xsep1)
    # There is probably a better way to do this
    # for i in range(npts):
        # draw a line from xsep1,ysep1 to xsep2,ysep2 upscaling the length
    # Use annotate to draw an arrow plot where I have control over the scale,
    # and a scale arrow
    ascale = 0.6
    xscale_label = xlim2 * 0.9
    # Depending whether you print label left or right of scale arrow
    # xscale_label = xlim2 * 0.4
    yscale_label = ylim1 * 0.95
    arrow_plot(xsep1.deg, ysep1.deg, xsepdiff.arcsec, ysepdiff.arcsec, scale=ascale, plotradius=maxradius, color='blue')
    arrowlength = 0.3*ascale
    plt.annotate('', xytext=(xscale_label,yscale_label), xy=(xscale_label-arrowlength,yscale_label), arrowprops=dict(color='black', width=0.5, headlength=5, headwidth=3) )
    # plt.text(xscale_label+0.02,yscale_label,'0.3" apparent motion',horizontalalignment='right',verticalalignment='center')
    plt.text(xscale_label-arrowlength-0.02,yscale_label,'0.3" apparent motion',horizontalalignment='left',verticalalignment='center')

    # sinth and costh are the angle between the v-axis (up) and the y-axis (Dec)
    # I want to draw an arrow in the v direction, meaning du=0 and dv=length
    # use the transform of 0,v to x,y as in derotate_grid.
    # For some reason, the arrow direction needed to have x flipped, prob due to
    # mix up because u (az) increases to right and x (dec) increases to left, 
    # after I fixed that in derotate_grid the up-arrow seems to be ok.
    sinth1, costh1, sinth2, costh2 = angle_list
    vlength = 0.25
    xuparrow1 = -vlength * sinth1
    # xuparrow1 = -xuparrow1
    yuparrow1 = vlength * costh1
    xuparrow2 = -vlength * sinth2
    # xuparrow2 = -xuparrow2
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

##

# fieldcen = get_coords()
# tstart, tend, exptime = get_times()

# nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2 = calc_exposure_offsets(fieldcen, tstart, tend)

def main():
    while True:
        fieldcen = get_coords()
        if fieldcen == -1:
            break
        tstart, tend, exptime = get_times()
        # plt.close()
        nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, angle_list, sep_rms = calc_exposure_offsets(fieldcen, tstart, tend)
        plot_exposure_offsets(nra, ndec, index_cen1, xsep1, ysep1, xsep2, ysep2, az, alt, exptime, angle_list, sep_rms)
    print("quitting")
    return
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()


###  This is old stuff and notes
'''

############
        
# I think what's happening is that the atmo refraction is directed up along the local
# vertical, which is the same as the up/down in dec/alt in the center column of the grid,
# but offset in angle on the left/right sides of the the grid. So it pulls the
# grid out of square, pulling the L/R edges down (less up than the center col),
# pulling the bottom edge inward, and top edge outward.  This is pretty hard to visualize.

#####
# print(u1)
# print(u1[78],u1[84],u1[90])
# -0d00m00.3613s 0d00m00s 0d00m00.3555s
# print(u1[6],u1[84],u1[162])
# 0d59m59.1429s 0d00m00s -0d59m59.14s
# print(v1)
# print(v[0],v[6],v[12])
# -1d00m32.2828s -0d00m35.3002s 0d59m21.8037s
# print(v[78],v[84],v[90])
# -0d59m57.0782s 0d00m00s 0d59m57.228s
# u is very close to the evenly spaced grid from -1 to 1 deg, mostly within ~1 arcsec.
# I think it's at slightly smaller offset <1deg, because the refraction pushes objects 
# to higher elevation where the arc diameter of the circle of azimuth is smaller.
# (Meeus mentions this in context of apparent horiz diam of Sun, bottom of p.107.)
# the v differences are larger, in the center they are near-symmetrical and the
# grid is compressed by about 3 arcsec over 1 degree; this I understand.
# at the extremes of u (leftmost/rightmost) v is about -35 arcsec to south,
# compared to the original spacing, and I have a little trouble understanding
# why it's to south.  It might be because I ignored the tan reprojection onto the
# center, instead just using az/alt as the coordinates.

# astropy can compute separations, which is what I want, but not sure it can
# compute _apparent_ separations.  However, I could treat apparent az and alt
# as if they are an ra and dec, and compute separations.  Also, Meeus Chap 17
# gives formulae for relative separation and PA between two points given ra and dec.

    

# Rotation between the az-alt and ra-dec systems - But do I want apparent ra/dec
# not true ra/dec?

ra1  = field_grid_ra[index_cen]
dec1 = field_grid_dec[index_cen]
ra2  = field_grid_ra[index_south]
dec2 = field_grid_dec[index_south]
ra3  = field_grid_ra[index_north]
dec3 = field_grid_dec[index_north]
# or
ra1  = field_grid[index_cen].ra
dec1 = field_grid[index_cen].dec
ra2  = field_grid[index_south].ra
dec2 = field_grid[index_south].dec
ra3  = field_grid[index_north].ra
dec3 = field_grid[index_north].dec
# az1  = field_grid_altaz_pkpno[index_cen].az
# alt1 = field_grid_altaz_pkpno[index_cen].alt
# az2  = field_grid_altaz_pkpno[index_north].az
# alt2 = field_grid_altaz_pkpno[index_north].alt
# print(ra1, dec1, az1, alt1)
# 13h00m00s -25d00m00s 179d57m47.3244s 32d55m06.5896s
# print(ra2, dec2, az2, alt2)
# 13h00m00s -24d00m00s 179d57m47.7528s 33d55m03.8176s
utmp1 = u3[index_cen]
vtmp1 = v3[index_cen]
utmp2 = u3[index_south]
vtmp2 = v3[index_south]
utmp3 = u3[index_north]
vtmp3 = v3[index_north]

# Now I can get the rotation angle between these systems
dra21  = ra2 - ra1
ddec21 = dec2 - dec1
# daz  = az2 - az1
# dalt = alt2 - alt1
# print(dra, ddec, daz, dalt)
# 0h00m00s 1d00m00s 0d00m00.4284s 0d59m57.228s
# for the example time since it's nearly transiting the paral angle is near 0.
# dra is 0 by construction (in true ra/dec)
# sintheta = daz / ddec
# costheta = dalt / ddec
dutmp21 = utmp2 - utmp1
dvtmp21 = vtmp2 - vtmp1
sintheta21 = dutmp21 / ddec21
costheta21 = dvtmp21 / ddec21
print(sintheta21, costheta21)
# 9.928230432755806e-05 0.9991883856290287
dra31  = ra3 - ra1
ddec31 = dec3 - dec1
dutmp31 = utmp3 - utmp1
dvtmp31 = vtmp3 - vtmp1
sintheta31 = dutmp31 / ddec31
costheta31 = dvtmp31 / ddec31
print(sintheta31, costheta31)
# 9.991154972948647e-05 0.9992300054426211
# Although the sign of the diffs is reversed because south-cen vs north-cen,
# the diffs in uv and radec are both reversed, so the ratio has same sign,
# but the angle of rotation inferred is slightly different depending on
# north-cen vs south-cen, because the grid is distorted. For these values of
# sintheta 9.928e-5 and 9.991e-5, theta~=sintheta, these are ~ 20.5" and 20.6".
# The difference is 0.06e-5 in sintheta, so at a radius of 1 degree, the
# diff in rotation makes an offset of 0.002 arcsec, so not worth worrying about.
# We do need to keep track of the rotation (between u-v and x-y), it's only
# small here because the test was nearly on the meridian.

# Doing it in az and alt was wrong  because az is in longitude not degrees of arc,
# and we are at elev +30 or so. So if the edge of the grid is 1 deg of arc,
# it's about 1.15 deg of azimuth, and mixing the units may do weird things.

# daz_grid  = field_grid_altaz_pkpno.az - field_grid_altaz_pkpno[index_cen].az
# dalt_grid = field_grid_altaz_pkpno.alt - field_grid_altaz_pkpno[index_cen].alt
# dra_grid  = daz_grid * costheta - dalt_grid * sintheta
# ddec_grid = daz_grid * sintheta + dalt_grid * costheta
# u, v, x, y are already centered at the field center

dx_grid = u1 * costheta - v1 * sintheta
dy_grid = u1 * sintheta + v1 * costheta

# dra_grid and ddec_grid weren't really ra and dec, they were
# angular displacements along the directions of the ra and dec axes, but
# they are apparent with all the refraction change included. Also dra is
# not compensated for the cos(dec) compression, and we're at dec -25.
# So I changed those to x and y, or dx_grid and dy_grid to represent
# that they are offsets from field center

# Now I could do the same calculation at two different times for same field
# and see how much dra_grid and ddec_grid changed.  These are offsets relative
# to the center, where I'm assuming the guiders will drive the field to
# be at what they think is the correct field center and position angle.

# Calc apparent coords over grid ?

# Apply refraction?  Calc size of the effect

# In DESI, dither experiments seem to see offsets that vary from top to bottom of
# the field, but they're not straight up and down, and they're not completely aligned
# with the x-y axes, the line of zero offset was tilted by some degrees
# (I think this was in the az-alt plane; in any case it was done when the field was
# around transit and RA/Dec is nearly az/alt)
# That makes it more complex than a mere miscalculation of refraction, which should
# vary only along the alt axis.  Also when the dithers are pointed near the zenith,
# the offsets are small, which makes it seem like refraction is involved.
# What if precession is involved somehow?  Like the calc of refraction was done
# for an az-alt calculated from an unprecessed RA-Dec when it should have been
# precessed, and then the refraction offsets were applied to the actual az-alt ?
# With astropy skycoord I can't make it do a "wrong" calculation, but I could
# compute the offsets over the field caused by refraction for the correct position
# and for an incorrect position, and see if those differ by some pattern that
# looks like the observed issue.

# But shouldn't the main effect of precession be to offset the field by a few
# min of time and of dec? How much diff should that make in the differential
# corrections across the field, and wouldn't it be less than integrating on the
# field for >15 minutes or however long it takes?  Is it possible that the
# pattern they see is just a consequence of field rotation during the exposure,
# and perhaps computing the positions at the start of the exposure rather
# than at the midpoint?



# What happens if we calculate it for some other position, like
# unprecessed RA, Dec, and then look at how much that changes the
# positions?  We could look at the positions as a function of
# offset from the field center, assuming that the guiders will
# always get the center onto the desired center, and we need to know
# the difference in relative position. Probably need to use
# az-alt and parallactic angle to transform into some X-Y focal plane
# and compare that, due to field rotation



# old stuff?

'''
