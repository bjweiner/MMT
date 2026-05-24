#!python
#
#  Extract some star spectra from individual exposures
#  from Binospec, use the sci_img_*_proc.fits as input
#
#  The input file format is:
#    objid
#    date
#    input_image_name.fits   output_extracted_1dspec_name.fits
#    ...
#  with as many rows of input_image output_spec as you want for a
#  given target.  The objid and date are just used to annotate the plot
#  and make the output plot filename.

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

if len(sys.argv) >= 2:
    fname = sys.argv[1]
else:
    fname = raw_input('Enter filename with obj, date, list of files: ')

f = open(fname, 'r')
line = f.readline()
objid = line.split()[0]
line = f.readline()
datestr = line.split()[0]

imnames = []
outnames = []
for line in f:
    if (line.strip() != '' and line[0] != '#'):
        imname1, outname1 = line.split()
        imnames.append(imname1)
        outnames.append(outname1)
f.close()

# Hard set the extraction location for now

# normal location of a single object longslit side B
ycen = 2081
# This was for a serendip object
# ycen = 2575

# to prompt for a y location to extract
yloc_str = raw_input('Enter y pixel location to extract: ')
yloc = int(yloc_str)
ycen = yloc

yrad = 10
ysky = ycen - 2 * yrad
nrows = 2*yrad

# Very bogus wavelength calibration, cwl is wavelength at x~2048
# I use a bogus calibration because we are extracting from raw data
# for G270 / 6500
# cwl = 6500.0
# dwdpix = 1.3
# for G600 / 6500
cwl = 6500.0
dwdpix = 0.6
# for G600 / 8500
# cwl = 8500.0
# dwdpix = 0.6

# to prompt for a wave calib
wcal_str = raw_input('Enter central wavelength and dw/dpix: ')
wcal_info = wcal_str.split()
cwl = float(wcal_info[0])
dwdpix = float(wcal_info[1])

plotsty = ['k-', 'm-', 'b-', 'c-', 'g-', 'r-']
nsty = len(plotsty)
# For stars we likely want to subtract the sky, but to extract
# from a flat field we don't want the subtracted version
# (and we could have averaged the extraction rather than sum)

fig, ax = plt.subplots()

# loop through spectra plotting each
nspec = len(imnames)
for i in range(nspec):
    hdulist1 = fits.open(imnames[i])
    # get a header from 1st image
    if i == 0:
        # For some reason the value of WSEEING may be NA in extension 2
        # header, so use extension 1
        hdr2d_1 = hdulist1[ 1 ].header
        airmass = hdr2d_1['AIRMASS']
        catpa = hdr2d_1['CATPA']
        posang = hdr2d_1['POSANG']
        parang = hdr2d_1['PA']
        rotatorang = hdr2d_1['ROT']
        wseeing = hdr2d_1['WSEEING']
        offsetang = posang - parang
        # if offsetang < -360.0:`
        #    offsetang = offsetang + 360.0
        #if offsetang > +360.0:`
        #    offsetang = offsetang - 360.0
        #if offsetang < -180.0
        #    offsetang = offsetang + 180.0
        #if offsetang > +180.0:
        #    offsetang = offsetang - 180.0
        #tmp1 = int(offsetang / 180)
        #offsetang = offsetang - tmp1 * 180.0
        # Report offset between slit pa and parang in 0 to 180 deg
        offsetang = offsetang % 180
        # This offset should be ~equal to the rotator angle, modulo 180
    # get the data of side B, 2nd fits extension
    indx = 2
    data2d_1 = hdulist1[ indx ].data
    # extract spectrum of a section, summing along y-axis
    objspec = np.sum( data2d_1[ycen-yrad:ycen+yrad,:], axis=0)
    skyspec = np.sum( data2d_1[ysky-yrad:ysky+yrad,:], axis=0)
    objsubspec = objspec - skyspec
    nx = len(objspec)
    xpix = np.arange(nx)
    wavelen = cwl + (xpix - nx/2.0) * dwdpix
    hdulist1.close()
    # make an astropy table
    table_spec = Table([xpix, wavelen, objspec, skyspec, objsubspec], names=('x', 'wavelen', 'objspec', 'skyspec', 'objsubspec'))
    # write it as a fits table
    table_spec.write( outnames[i], format='fits', overwrite=True)
    plt.plot(wavelen,objspec,plotsty[ i % nsty ])

# save the plot with annotation

titlestr = objid + ' ' + datestr
plotname = objid + '_' + datestr + '.png'

print (airmass, offsetang, rotatorang, wseeing)

airmasslabel = 'airmass %4.2f' % (airmass)
offsetanglabel = 'delta-parang %6.1f' % (offsetang)
rotatoranglabel = 'rot angle %6.1f' % (rotatorang)
if (wseeing > 0.0) and (wseeing <5.0):
    seeinglabel = 'WFS seeing %4.2f' % (wseeing)
else:
    seeinglabel = 'WFS seeing NA'
# plt.text(7500, 6e4, airmasslabel) 
# plt.text(7500, 5e4, offsetanglabel) 
# Best to plot as a fraction of the x and y limits, forget how
# ax.text(0.7, 0.8, airmasslabel, transform=ax.transAxes)
# ax.text(0.7, 0.7, offsetanglabel, transform=ax.transAxes)
plt.figtext(0.67, 0.8, airmasslabel)
# plt.figtext(0.67, 0.75, offsetanglabel)
plt.figtext(0.67, 0.75, rotatoranglabel)
plt.figtext(0.67, 0.7, seeinglabel)

plt.xlabel('wavelength, approx')
plt.ylabel('extracted sum, counts')
plt.title(titlestr)

# to set plot limits for the flatfield plot
# plt.axis([3800,7000,0,150e3])

plt.savefig(plotname)



