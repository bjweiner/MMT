#!/usr/bin/env python

# Fit a line in a longslit spectrum like MMT Blue Channel
# in a series of images, reading an x-region to fit and a list
# of image names from the command line.  This allows quickly
# checking for wavelength shifts due to grating movement,
# flexure, etc.  - Ben Weiner, May 2022

# Requires numpy, astropy, and scipy to be installed.
# Run like:
#   python centroid_line.py list.nt1.object_5398 1500 1550
# with arguments: file w/list of FITS files, xmin, xmax
# The fit can go off to la-la land if the line is weak or there
# are multiple lines, so ignore fits where the stddev (line width)
# is very large, fit is offset from x-pixel of peak flux, etc.

import sys
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.stats import biweight
from scipy.signal import medfilt
from astropy.modeling import models, fitting

# Cheesy hard-coding x-regions for the moment to find a line
# for objects / night sky at cwl 4500 A
xlim1 = [1160, 1260]
# for objects / night sky at cwl 5500 A
xlim2 = [1370, 1430]
# for comps at cwl 4500 A
xlim3 = [790,860]
# for comps at cwl 5500 A
xlim4 = [1300,1380]
# Set which one to use
xlim = xlim4
# take average over rows in upper region of slit, away from target
# This assumes your images are about 256 rows; if you change the
# readout region 
ylim = [180, 240]

# allow passing in an x-region with command line arguments, as
# python centroid_line.py fname xlim1 xlim2
#
if len(sys.argv) >=4 :
    xlim = [ int(sys.argv[2]), int(sys.argv[3]) ]
    # print("Using xregion: ",xlim[0], xlim[1])

# open image file and retrieve data.  Assume data is in hdu #1 or 0?
# select x-section of data to analyze

def get_region(fname, xlim, ylim):
    xmin, xmax = xlim
    ymin, ymax = ylim
    hdulist = fits.open(fname)
    hdunum = 0
    pixdata = hdulist[hdunum].data
    # Python shapes the data with rows (y-col) first I guess
    ny, nx = pixdata.shape
    # print("nx, ny = ", nx, ny)
    if ((xmax > nx) or (ymax > ny)):
        print("Requested data region out of bounds, nx, ny:",nx,ny)
        xmax = min(xmax, nx)
        ymax = min(ymax, ny)
    pixregion = pixdata[ymin:ymax,xmin:xmax]
    hdulist.close()
    return pixregion


# take median in columns to filter out cosmic rays etc?
# Median doesn't work well here because the data is quantized integer

def filter_columns_fast(data_2d):
    # axis=0 should take the median/mean along the y-axis ?
    # print("data 2d shape ", data_2d.shape)
    data_1d = np.mean(data_2d, axis=0)
    # print("data 1d shape ", data_1d.shape)
    return data_1d

def filter_columns(data_2d):
    # for each column, compute some robust average like the biweight
    # print("data 2d shape ", data_2d.shape)
    ny, nx = data_2d.shape
    smdata = np.zeros(nx)
    for i in range(nx):
        coldata = data_2d[:,i].astype(float)
        bwt = biweight.biweight_location(coldata)
        smdata[i] = bwt        
    # print("data 1d shape ", smdata.shape)
    return smdata

# compute flux-weighted centroid along a row and return. This is
# too dependent on the input region to be useful, or it needs a
# first step to center the region on the line and a dc-level
# subtraction?
def row_centroid(data_1d):
    # subtract the DC level - this doesn't work well because
    # the denominator of sum(x*flux)/sum(flux) becomes zero
    # foo = np.mean(data_1d)
    # print("dc level = ",foo)
    # data_1d = data_1d - foo
    nx = data_1d.shape[0]
    x_1d = np.arange(nx)
    sum_flux = np.sum(data_1d)
    sum_xflux = np.sum(x_1d * data_1d)
    xcen = sum_xflux / sum_flux
    return xcen

# smooth row with a box median filter using scipy.signal.medfilt
def smooth_row(data_1d):
    filtsize = 3
    data_out = medfilt(data_1d, filtsize)
    return data_out

# When the data has been smoothed along each column, the location
# of the max appears to be good enough to find the line. 
# a light median filter along the row helps this too.
def row_max(data_1d):
    # nx = data_1d.shape[0]
    xcen = np.argmax(data_1d)
    return xcen

# Fit a gaussian to a line, passing in an initial guess for the
# location
def fit_gaussian(data_1d, x_loc=0.0):        
    nx = data_1d.shape[0]
    if x_loc < 0.01:
        x_loc = nx/2.0
    x_1d = np.arange(nx)
    # subtract the mean or median, although it may be better to have
    # the model fit a continuum with Const1D ?
    datamed = np.median(data_1d)
    # data_sub = data_1d - np.median(data_1d)
    gauss_init = models.Gaussian1D(amplitude=1.0, mean=x_loc, stddev=1.5) + models.Const1D(amplitude=datamed)
    fit_g = fitting.LevMarLSQFitter()
    gauss = fit_g(gauss_init, x_1d, data_1d)
    return gauss

# read list of input files from a file
def do_list_files(listname):
    # print a header line
    print("# filename     max_xpixel  fit_xpixel  amplitude  stddev")
    flist = open(listname,'r')
    for line in flist:
        fn = line.strip()
        data_2d = get_region(fn, xlim, ylim)
        data_1d = filter_columns(data_2d)
        smdata = smooth_row(data_1d)
        xmax_relative = row_max(smdata)
        gauss_fit = fit_gaussian(smdata, x_loc=xmax_relative)
        xbegin = xlim[0]
        xmax = xmax_relative + xbegin
        # mean_0 is for the 0th component of a compound model
        xfit = gauss_fit.mean_0.value + xbegin
        # print '%s %s %4d %s %7.2f %s %7.2f %s %5.2f' % (fn, " x pixel of flux max = ", xmax, "fit line loc = ", xfit, " ampl= ",gauss_fit.amplitude_0.value, " stddev= ",gauss_fit.stddev_0.value)
        # Leave out the text in each line; print a header at beginning
        print ('%s       %4d      %7.2f      %7.2f    %5.2f' % (fn, xmax, xfit, gauss_fit.amplitude_0.value, gauss_fit.stddev_0.value))
    flist.close()
    
def main():
    if len(sys.argv) >= 2:
        fname = sys.argv[1]
    else:
        fname = raw_input('Enter filename with list of spectrum images: ')
    do_list_files(fname)

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()


# for testing

# fn = fname1
# # fn = fname2
# data_2d = get_region(fn, xlim, ylim)
# data_1d = filter_columns(data_2d)
# smdata = smooth_row(data_1d)
# xcen_relative = row_centroid(smdata)
# xmax_relative = row_max(smdata)
# gauss_fit = fit_gaussian(smdata, x_loc=xmax_relative)
# xbegin = xlim[0]
# xmax = xmax_relative + xbegin
# print(xcen_relative,xmax_relative)
# print(fn, " xmax = ", xmax)
# print(gauss_fit.amplitude, gauss_fit.mean, gauss_fit.stddev)
# 
# flist1 = open('list.nt1.object_files','r')
