#!python
#
#  For checking whether observations at a given PA (like Binospec masks)
#  will have a problem with the rotator limits
#  prompt for an observation date
#  Read a table of Name, RA, Dec, PA  columns (space-delimited)
#  make plots of airmass, parallactic angle, and rotator angle
#  Calling:
#    python rotator_angle.py <masks_name_radec_pa.dat> <instname>
#  where masks_name_radec_pa.dat is the name of your file with data
# and instname is 'binospec' or 'mmirs' (default binospec)
# If you omit the arguments, then it will prompt for the filename,
# but then you can't specify instrument.

# Ben Weiner Feb -April 2018


# requires pyephem

import ephem
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Read RA, Dec, PA columns from file as strings, return as lists.
# RA and Dec must be something that pyephem parses, so hh:mm:ss dd:mm:ss works,
# but decimal RA in degrees doesn't. However decimal RA in hours and
# decimal Dec in degrees do work.
def read_coords(fname):
    names = []
    ra = []
    dec = []
    pa = []
    f = open(fname,'r')
    for line in f:
        if (line.strip() != '' and line[0] != '#'):
            name1, ra1, dec1, pa1 = line.split()
            names.append(name1)
            ra.append(ra1)
            dec.append(dec1)
            pa.append(pa1)
    f.close()
    return names, ra, dec, pa

# Hard-code the observatory location for the moment.
# just set to Hopkins for now
#
def set_observatory():
    hopkins = ephem.Observer()
    hopkins.lon = '-110.883'
    hopkins.lat = '31.688'
    hopkins.elevation = 2608

    bigelow = ephem.Observer()
    bigelow.lon = '-110:44:04.3'
    bigelow.lat = '32:24:59.3'
    bigelow.elevation = 2510

    bok = ephem.Observer()
    bok.lon = '-111:36:01.6'
    bok.lat = '31:57:46.5'
    bok.elevation = 2071

    lemmon = ephem.Observer()
    lemmon.lon = '-110.7889'
    lemmon.lat = '32.44257'
    lemmon.elevation = 2805

    house = ephem.Observer()
    house.lon = '-110:55.274'
    house.lat = '32:13.846'
    house.elevation = 753

#    return bigelow
#    return house
#    return lemmon
    return hopkins
#    return bok

def get_a_time():
    tnow = ephem.now()
    localhms = raw_input('Enter UT date to start compute at: yyyy mm dd hh mm ss [now]: ')
    if len(localhms) > 1:
        localtime = localhms.split()
        # print localtime
        trequest_tuple = (int(localtime[0]), int(localtime[1]), int(localtime[2]), int(localtime[3]), int(localtime[4]), float(localtime[5]))
        trequest = ephem.Date(trequest_tuple)
    else:
        trequest = tnow
    # print trequest
    return trequest

# name, ra and dec are strings, pa is a string with degrees
def set_body(name, ra, dec, pa):
    body = ephem.FixedBody()
    body.name = name
    body._ra = ra
    body._dec = dec
    body._epoch = '2000/1/1 12:00:00'
    # Something weird happens here that causes positive PA to be
    # off by ~ +0.3 or -0.3 deg, and negative PA to be off by ~ +3.0
    # deg, closer to 0.  Maybe pyephem bug?
    pa_pos = float(pa)
    if pa_pos < 0.0:
        pa_pos = pa_pos + 360.0
    body._pa = 3.1416 / 180.0 * pa_pos
    return body

def compute_pos_at_time(body, obs_location, time):
    obs_time = obs_location
    obs_time.date = time
    body.compute(obs_time)
    return body

# compute parallactic and rotator angles and normalize the rot angle
# back into -180 to +180
def compute_rotangle_at_time(body, obs_location, time):
    obs_time = obs_location
    obs_time.date = time
    body.compute(obs_time)
    parang = body.parallactic_angle()
    rotang = ephem.degrees(parang - body._pa)
    rotang_norm = rotang.znorm
    if body.alt > 0.0:
        airmass = 1.0 / np.sin(body.alt)
    else:
        airmass = 1000.0
    return parang, rotang_norm, airmass

    
# pyephem does compute parallactic angles, see the def of FixedBody
# See atmo_refraction.py for the formulae for PA
# def compute_parangle(body, obs_location, time):

# Do a loop over time through the night, computing angles

def compute_angle_timeloop(body, observatory, timestart, duration=12.0):
    stephour = 0.05
    nstep = int (duration / stephour)
    hourtime = np.zeros(nstep)
    parang = np.zeros(nstep)
    rotang = np.zeros(nstep)
    airmass = np.zeros(nstep)
    for i in range(nstep):
        time1 = ephem.Date(timestart + i * stephour * ephem.hour)
        parang1, rotang1, airmass1 = compute_rotangle_at_time(body, observatory, time1)
        y, m, d, h, m, s = time1.tuple()
        hourtime[i] = h + m/60.0 + s/3600.0
        parang[i] = parang1
        rotang[i] = rotang1
        airmass[i] = airmass1
    return hourtime, parang, rotang, airmass

def plot_angles(body, observatory, timestart, rotlimits, plotnum=0, pdf_file=''):
    # compute angles as function of time
    # hourtime, parang, rotang, airmass = compute_angle_timeloop(body, observatory, timestart)
    # for a longer plot in winter
    hourtime, parang, rotang, airmass = compute_angle_timeloop(body, observatory, timestart, duration=14.0)
    # date for x label
    ymd_string = str(timestart).split()[0]
    # plot title with mask, ra, dec, pa. body._pa is just a float not an Angle for some reason.
    panorm = 180.0/3.1416 * body._pa
    if panorm > 180.0:
        panorm = panorm-360.0
    print body.name, panorm
    titlestring = '%s, %5s %6s, PA=%5.0f' % (body.name, str(body._ra)[0:5], str(body._dec)[0:6], panorm )
    # x plot limits
    hourmin = hourtime[0]
    hourmax = hourtime[len(hourtime)-1]
    ifup = (airmass < 10)
    # find 12 degree twilight times
    sun = ephem.Sun()
    observatory.horizon = '-12:00'
    observatory.date = timestart
    twi1 = observatory.next_setting(sun)
    twi2 = observatory.next_rising(sun)
    y,m,d,h,m,s = twi1.tuple()
    twi1_hour = h + m/60.0 + s/3600.0
    y,m,d,h,m,s = twi2.tuple()
    twi2_hour = h + m/60.0 + s/3600.0
    # make plot
    plt.clf()
    fig, ax1 = plt.subplots()
    fig = ax1.axis([hourmin,hourmax,-180,180])
    #xstring = body.name + ', UT time, hours'
    xstring = ymd_string + ', UT time, hours'
    fig = ax1.set_xlabel(xstring)
    fig = ax1.set_ylabel('angle, deg')
    # plot angles as fn of time
    ax1.plot(hourtime[ifup], 180.0/3.1416*parang[ifup], 'k:')
    ax1.plot(hourtime[ifup], 180.0/3.1416*rotang[ifup], 'b-')
    #ax1.text(hourmin+0.5,155.0,body.name,color='k')
    #plt.title(str(body.name))
    plt.title(titlestring)
    ax1.text(hourmin+0.5,140.0,'parallactic angle',color='k')
    # Horizontal lines for rotator limits
    ax1.plot([hourmin, hourmax], [rotlimits[0],rotlimits[0]], 'r-')
    ax1.plot([hourmin, hourmax], [rotlimits[1],rotlimits[1]], 'r-')
    ax1.text(hourmin+0.5,120.0,'rotator angle',color='b')
    ax1.text(hourmin+0.5,-165.0,'rotator limits',color='r')
    # Vertical lines for twilight
    ax1.plot([twi1_hour,twi1_hour], [-180,180], color='c')
    ax1.plot([twi2_hour,twi2_hour], [-180,180], color='c')
    ax1.text(twi1_hour,-140.0,'12 deg twilight',color='c')
    # Second y axis for plotting airmass
    ax2 = ax1.twinx()
    fig = ax2.axis([hourmin,hourmax,0,3])
    fig = ax2.set_ylabel('airmass')
    ax2.plot(hourtime[ifup], airmass[ifup], 'g-')
    ax2.text(hourmax-0.5,2.6,'airmass',color='g',horizontalalignment='right')
    if pdf_file == '':
        # Individual pdf files per object
        # plt.show()
        fname = 'rotangle_plot_%d.pdf' % plotnum
        plt.savefig(fname)
    else:
        # All in one pdf file
        pdf_file.savefig()
        plt.close()
    return hourtime, parang, rotang, airmass

# Loop through sets of coordinates making plots
# Rotator limits are hard coded; now switchable on instrument name

def plot_angle_loop(names,ra,dec,pa,instname):
    observatory = set_observatory()
    tstart = get_a_time()
    rotlimits_raw = [-174, +171]
    # Binospec has a -8 deg rotator offset that shifts the limits to -166,+179
    # MMIRS is +7,; limits are tabulated as -174,+165 in doc database,
    # so it doesn't work to apply the offset to the same raw limits.
    binospec_limits = [-166, +179]
    mmirs_limits = [-174, +165]
    if instname.lower() == 'binospec':
        rotlimits = binospec_limits
    elif instname.lower() == 'mmirs':
        rotlimits = mmirs_limits
    else:
        rotlimits = binospec_limits
    pdfname = 'rotangle_plots.pdf'
    pdffile = PdfPages(pdfname)
    nobj = len(ra)
    for i in range(nobj):
        body = set_body(names[i], ra[i], dec[i], pa[i])
        hourtime, parang, rotang, airmass = plot_angles(body, observatory, tstart, rotlimits, plotnum=i, pdf_file=pdffile)
    pdffile.close()
    print 'Wrote plots to file ',pdfname
    return


def main():
    if len(sys.argv) >= 2:
        fname = sys.argv[1]
    else:
        fname = raw_input('Enter filename with name, ra, dec, pa (space delimited, hh:mm:ss): ')
    if len(sys.argv) >= 3:
        instname = sys.argv[2]
    else:
        instname = 'binospec'
    names, ra, dec, pa = read_coords(fname)
    plot_angle_loop(names, ra, dec, pa, instname)
    return

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()


    
                  
