#!python
#
#  Differential atmospheric refraction across some field size,
#  like a slitmask diameter
#
#  The Bennett and Saemundsson formulae used here are from the wikipedia
#  article and are for visual so about 500-550 nm.
#  See R.C. Stone 1996 http://adsabs.harvard.edu/abs/1996PASP..108.1051S
#  for a much more detailed treatment.
#
#  For par angles etc see also Tom's pages at
#  http://cholla.mmto.org/astronomy/
#  http://cholla.mmto.org/astronomy/altaz.html
#
#  And see J. Meeus, Astronomical Algorithms chap 13-16
#
#  https://www2.keck.hawaii.edu/inst/common/parallactic.html
#  Filippenko 1982 plots of par angle vs hour angle
#  http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1982PASP...94..715F

import numpy as np
import matplotlib.pyplot as plt

degtorad = np.pi / 180.0

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

azimuth = np.arange(360)
elevation = np.arange(1,90) 
elevrad = elevation * degtorad 

refrac_rad = refrac_true(elevrad)

# Default pressure = 101.0 kPa, temp = 10 deg C
pressure = 101.0
temp = 10.0
refrac = (pressure / 101.0) * 283.0 / (273.0 + temp) * refrac_rad

refrac_sec = refrac / degtorad  * 3600.0

# Use 0.25 deg for Binospec and 1.0 deg for Hectospec, or other
# truly wide field instruments (DESI = 1.6 deg?)
# don't forget to change the differential refrac y axis label and plot limit
# fieldsize_deg = 0.25
fieldsize_deg = 1.0
fieldsize_min = fieldsize_deg * 60.0
fieldsize_rad = fieldsize_deg * degtorad

refrac_diff_rad = refrac_true(elevrad) - refrac_true(elevrad + fieldsize_rad)
refrac_diff = (pressure / 101.0) * 283.0 / (273.0 + temp) * refrac_diff_rad
refrac_diff_sec = refrac_diff/ degtorad * 3600.0

plt.plot(elevation,refrac_sec/60.0,'g-')
fig = plt.axis([0,90,0,20])
fig = plt.xlabel('Elevation, deg')
fig = plt.ylabel('Refraction, arcmin')
# plt.savefig('elev_refrac_min.pdf')
plt.savefig('elev_refrac_min.png')

plt.clf()
plt.plot(elevation,refrac_diff_sec,'b-')
# fig = plt.axis([0,90,0,3])
fig = plt.axis([0,90,0,6])
fig = plt.xlabel('Elevation, deg')
# fig = plt.ylabel('Differential refraction across 15\', arcsec')
fig = plt.ylabel('Differential refraction across 60\', arcsec')
# plt.savefig('elev_diff_refrac.pdf')
plt.savefig('elev_diff_refrac.png')

ii = np.where(np.abs(elevation % 5) < 0.01)
print "elevation = ",elevation[ii]
print "Diff refrac = ", refrac_diff_sec[ii]

# From Tom Trebisky's webpage. Mt Hopkins is at +7:23:32 hr W, +31:41.3 deg N
lat_deg = 31 + 41.3/60
lat_rad = lat_deg * degtorad
ha_deg = np.arange(-120,120,1) 
ha_hr = ha_deg / 15.0
ha_rad = ha_hr * 15 * degtorad
dec_array_deg = np.arange(-40,90,10)
linestyle_array = ['k-','m-','b-','c-','g-','y-','r-','k-','m-','b-','c-','g-','y-','r-']
linestyle_array = linestyle_array[0:len(dec_array_deg)]

# Calculate az, alt, parang, in degrees, for a scalar latitude and
# declination and a vector of hour angle
def calc_az_alt_pa_for_dec(lat_deg, ha_hr, dec_deg):
    degtorad = np.pi / 180.0
    lat_rad = lat_deg * degtorad
    ha_rad = ha_hr * 15 * degtorad
    dec_rad = dec_deg * degtorad

    use_trebisky = True

    if use_trebisky:
        # Trebisky webpage
        sinh = np.sin(lat_rad) * np.sin(dec_rad) + np.cos(lat_rad) * np.cos(dec_rad) * np.cos(ha_rad)
        sinz = -np.cos(dec_rad) * np.sin(ha_rad)
        cosz = ( np.sin(dec_rad) - np.sin(lat_rad)*sinh) / np.cos(lat_rad)
        sin_alt = sinh
        # Trebisky webpage has this backward, tan_az should be sinz/cosz
        # tan_az = cosz / sinz
        alt_rad = np.arcsin(sin_alt)
        az_rad = np.arctan2(sinz, cosz)
        sinb = -np.sin(az_rad) * np.cos(lat_rad) / np.cos(dec_rad)
        cosb = -np.cos(ha_rad) * np.cos(az_rad) - np.sin(ha_rad) * np.sin(az_rad) * np.sin(lat_rad)
        # Same tan problem?
        tan_pa = cosb / sinb
        # pa_rad = np.arctan(tan_pa)
        # This one should be correct?
        # pa_rad = np.arctan2(cosb, sinb)
        # This one should be wrong?
        pa_rad = np.arctan2(sinb, cosb)
    else:
        # Using Meeus
        tan_az = np.sin(ha_rad) / (np.cos(ha_rad) * np.sin(lat_rad) - np.tan(dec_rad) * np.cos(lat_rad))
        sin_h = np.sin(lat_rad) * np.sin(dec_rad) + np.cos(lat_rad) * np.cos(dec_rad) * np.cos(ha_rad)
        sin_alt = sin_h
        alt_rad = np.arcsin(sin_alt)
        # az_rad = np.arctan(tan_az)
        sinz = -np.cos(dec_rad) * np.sin(ha_rad)
        cosz = ( np.sin(dec_rad) - np.sin(lat_rad)*sin_h) / np.cos(lat_rad)
        az_rad = np.arctan2(sinz,cosz)
        denom = (np.tan(lat_rad) * np.cos(dec_rad) - np.sin(dec_rad) * np.cos(ha_rad))
        tan_pa = np.sin(ha_rad) / denom
        pa_rad = np.arctan(tan_pa)
        # Fixing the arctan jumping from -90 to +90
        test1 = (ha_rad > 0) & (dec_rad > lat_rad) & (pa_rad > 0)
        pa_rad[test1] = pa_rad[test1] - np.pi
        test2 = (ha_rad < 0) & (dec_rad > lat_rad) & (pa_rad < 0)
        pa_rad[test2] = pa_rad[test2] + np.pi

    # Make pa go from 0 to 360, vs -180 to 180
    testneg = pa_rad < 0
    pa_rad[testneg] = pa_rad[testneg] +  np.pi
            
    alt_deg = alt_rad / degtorad
    az_deg = az_rad / degtorad
    pa_deg = pa_rad / degtorad
    return az_deg, alt_deg, pa_deg

plt.clf()
# fig = plt.axis([-8,8,-180,180])
fig = plt.axis([-8,8,0,180])
# fig = plt.axis([-8,8,-90,90])
fig = plt.xlabel("Hour angle, hours")
fig = plt.ylabel("Parallactic angle, deg")
for dec, sty in zip(dec_array_deg, linestyle_array):
    az_deg, alt_deg, pa_deg = calc_az_alt_pa_for_dec(lat_deg, ha_hr, dec)
    ifup = (alt_deg > 0)
    ifdown = (alt_deg <= 0)
    # can make parang go from -90 to 90 instead of 0 to 180 ?
    # This does not work well
    # testpa = pa_deg > 90.0
    # pa_deg[testpa] = pa_deg[testpa] - 180.0
    # plt.plot(ha_hr[ifup], pa_deg[ifup], 'b-')
    plt.plot(ha_hr[ifup], pa_deg[ifup], sty)
    # plt.plot(ha_hr[ifdown], pa_deg[ifdown], 'c.')
    plt.text(ha_hr[100],pa_deg[100]+5,str(dec))
plt.savefig('hourangle_parangle_bydec.pdf')


plt.clf()
fig = plt.axis([-8,8,-180,180])
fig = plt.xlabel("Hour angle, hours")
fig = plt.ylabel("Azimuth angle, deg")
for dec, sty in zip(dec_array_deg, linestyle_array):
    az_deg, alt_deg, pa_deg = calc_az_alt_pa_for_dec(lat_deg, ha_hr, dec)
    ifup = (alt_deg > 0)
    ifdown = (alt_deg <= 0)
    plt.plot(ha_hr[ifup], az_deg[ifup], sty)
    # plt.plot(ha_hr[ifdown], az_deg[ifdown], 'c.')
    plt.text(ha_hr[100],az_deg[100]+5,str(dec))
plt.savefig('hourangle_azimuth_bydec.pdf')

plt.clf()
fig = plt.axis([-8,8,-30,90])
fig = plt.xlabel("Hour angle, hours")
fig = plt.ylabel("Altitude angle, deg")
for dec, sty in zip(dec_array_deg, linestyle_array):
    az_deg, alt_deg, pa_deg = calc_az_alt_pa_for_dec(lat_deg, ha_hr, dec)
    ifup = (alt_deg > 0)
    ifdown = (alt_deg <= 0)
    plt.plot(ha_hr[ifup], alt_deg[ifup], sty)
    plt.plot(ha_hr[ifdown], alt_deg[ifdown], 'c.')
    plt.text(ha_hr[100],alt_deg[100]+5,str(dec))
plt.savefig('hourangle_altitude_bydec.pdf')

 
    
