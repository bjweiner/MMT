#!python
# 
# read some spectra written by extract_spectra.py
# and combine them

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
# for median smoothing
import scipy.signal as scisig

# How many pixels to smooth by
nsmooth = 11
# Plot the object spectrum, or object -sky
ifplotskysub = True

def mycombine(x):
    result = np.nanmedian(x)
    # result = np.median(x)
    # result = np.nanmean(x)
    # result = np.mean(x)
    return result

if len(sys.argv) >= 2:
    fname = sys.argv[1]
else:
    fname = input('Enter filename with list of files: ')

outname = input('Enter combined-output fits filename: ')

f = open(fname, 'r')

specnames = []
for line in f:
    if (line.strip() != '' and line[0] != '#'):
        sname1 = line.split()
        # print(sname1, sname1[0])
        # split() returns a list so you need to get the first element
        specnames.append(sname1[0])
f.close()

print(specnames)

nspec = len(specnames)

# get the wavelength array of first spectrum
hdulist1 = fits.open(specnames[0])
data1 = hdulist1[1].data
wave_ref = data1['wavelen']
hdulist1.close()

nx = len(wave_ref)

sky_2d = np.zeros((nx, nspec))
obj_2d = np.zeros((nx, nspec))
objsub_2d = np.zeros((nx, nspec))

for i in range(nspec):
    hdulist1 = fits.open(specnames[i])
    data1 = hdulist1[1].data
    wave1 = data1['wavelen']
    skyspec1 = data1['skyspec']
    objspec1 = data1['objspec']
    objsubspec1 = data1['objsubspec']
    # do something about interpolating the fluxes onto wave_ref
    # store them somewhere
    sky_interp = np.interp(wave_ref, wave1, skyspec1)
    obj_interp = np.interp(wave_ref, wave1, objspec1)
    objsub_interp = np.interp(wave_ref, wave1, objsubspec1)
    sky_2d[:, i] = sky_interp
    obj_2d[:, i] = obj_interp
    objsub_2d[:, i] = objsub_interp
    hdulist1.close()
    
# do some kind of combine operation, or kernel smoothing
# with a kernel as fn of wavelength, etc

# if you wanted to do something like rescale the spectra
# before taking the median, could do that

# I'm not sure that combining the sky makes total sense and
# doing it separately means that objsub_combined is not
# equal to obj_combined - sky_combined, but may as well
# compute them all

sky_combined = np.zeros(nx)
obj_combined = np.zeros(nx)
objsub_combined = np.zeros(nx)

# using a for loop here is inefficient but ok for now

for j in range(nx):
    tmp1 = mycombine(sky_2d[j, :])
    sky_combined[j] = tmp1
    tmp1 = mycombine(obj_2d[j, :])
    obj_combined[j] = tmp1
    tmp1 = mycombine(objsub_2d[j, :])
    objsub_combined[j] = tmp1

# do something to write out the combined spectra - 
# same format as input

xpix = np.arange(nx)
# make an astropy table
table_spec = Table([xpix, wave_ref, obj_combined, sky_combined, objsub_combined], names=('x', 'wavelen', 'objspec', 'skyspec', 'objsubspec'))
# write it as a fits table
table_spec.write( outname, format='fits', overwrite=True)

# make a plot. Don't need plotsty if plotting only one
plotsty = ['k-', 'm-', 'b-', 'c-', 'g-', 'r-']
nsty = len(plotsty)

if ifplotskysub == True:
    spectoplot = objsub_combined
else:
    spectoplot = obj_combined
if nsmooth > 1:
    flux_smooth = scisig.medfilt(spectoplot, kernel_size=nsmooth)
else:
    flux_smooth = spectoplot
plt.plot(wave_ref, flux_smooth, plotsty[ 0 ])

# titlestr = objid + ' ' + datestr
# plotname = objid + '_' + datestr + '.png'
titlestr = outname
plotname = outname + '.png'

plt.xlabel('wavelength')
plt.ylabel('combined flux')
# plt.title(titlestr)

# to set fixed plot limits
# plt.axis([3800,8000,0,1000])

plt.savefig(plotname)



    
