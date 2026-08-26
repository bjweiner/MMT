#!python

# BJW May 2026
# Borrowing routines from the not-completed combine_bino_spectra.py

# Typical 1d filenames  obj_counts_slits_extr.fits
#   obj_abs_slits_extr.fits
# Typical 2d filenames  obj_counts_slits_lin.fits
#   obj_abs_slits_lin.fits

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys
# for median smoothing
import scipy.signal as scisig

# Read mask info structures from ext 4 and 5 of the spec1d HDU
def get_mask_info(spec1dhdu):
    info_a = spec1dhdu[4].data
    info_b = spec1dhdu[5].data
    return info_a, info_b

# Open spec1d file
def open_bino_spec1d(spec1dname):
    spec1dhdu = fits.open(spec1dname)
    return spec1dhdu

# Open spec2d file
def open_bino_spec2d(spec2dname):
    spec2dhdu = fits.open(spec2dname)
    return spec2dhdu

# Read spec1d data    
def get_bino_spec1d(spec1dhdu):
    # These are 2-d images containing the 1-d spectra as rows
    flux1d_all = spec1dhdu[1].data
    error1d_all = spec1dhdu[2].data
    hdr1d_all = spec1dhdu[1].header
    # Fields like the object name, time, etc are in 0th extension
    hdr_main = spec1dhdu[0].header
    return hdr_main, hdr1d_all, flux1d_all, error1d_all

def close_bino_file(hdulist):
    hdulist.close()
    return 


# Read spec2d data
# def get_bino_spec2d(spec2dhdu):
#    # Get number of extensions somehow
#    nobj = something
#    for i = 1, nobj:
#        data2d = spec2dhdu[i].data

# plot object in row iobj, 0-based
def plot_object(hdr_main, hdr1d, flux1d, error1d, iobj):
    # test whether iobj is < n rows of image
    ny, nx = np.shape(flux1d)
    if (iobj >= ny):
        print("data has ",ny," rows but requested object ",iobj)
        return
    # name = hdr1d['CATID']
    name = hdr_main['OBJECT']
    crpix1 = hdr1d['CRPIX1']
    crval1 = hdr1d['CRVAL1']
    cd1_1 = hdr1d['CD1_1']
    wave_index = np.arange(nx)
    wave1d = (wave_index - (crpix1-1)) * cd1_1 + crval1

    # iobj = 1
    fluxobj = flux1d[iobj,:]
    errobj = error1d[iobj,:]
    # Could also use np.nanmax() to avoid nans
    isvalid = np.where( ~np.isnan(fluxobj) )
    fluxmin = np.min(fluxobj[isvalid])
    fluxmax = np.max(fluxobj[isvalid])
    wavemingood = wave1d[isvalid[0][0]]
    wavemaxgood = wave1d[isvalid[0][-1]]

    print("Valid flux from wavelengths: ", wavemingood, wavemaxgood)
    # To prompt for a plot wavelength range
    # wstring = input('Enter wavelength min max to plot: ')
    # waveminplot = float(wstring.split()[0])
    # wavemaxplot = float(wstring.split()[1])
    waveminplot = wavemingood
    wavemaxplot = wavemaxgood
    fminplot = fluxmin
    # fminplot = 0.0
    fmaxplot = fluxmax

    smoostr = input('Smooth flux by N pixels: ')
    nsmooth = int(smoostr)
    flux_smooth = scisig.medfilt(fluxobj[isvalid], kernel_size=nsmooth)
    err_smooth = errobj[isvalid] / np.sqrt(nsmooth)

    wlabel = 0.65*(wavemaxplot-waveminplot) + waveminplot
    flabel = 0.8*(fluxmax-fluxmin) + fluxmin

    plt.clf()
    plt.axis([waveminplot, wavemaxplot, fminplot, fmaxplot])
    # plt.plot(wave1d[isvalid], errobj[isvalid], 'r-')
    # plt.plot(wave1d[isvalid], fluxobj[isvalid], 'k-')
    plt.plot(wave1d[isvalid], err_smooth, 'r-')
    plt.plot(wave1d[isvalid], flux_smooth, 'k-')
    plt.xlabel('wavelength, nm')
    plt.ylabel('flux')
    plt.text(wlabel, flabel, name)
    # Could add iobj to name for spectrum plot
    plotname = name + '_spec1d_plot.pdf'
    plt.savefig(plotname)
    plt.show(block=False)
    return plotname

def plot_spec_from_file(fname, iobj):
    spec1dhdu = open_bino_spec1d(fname)
    hdr_main, hdr1d, flux1d, error1d = get_bino_spec1d(spec1dhdu)
    plotname = plot_object(hdr_main, hdr1d, flux1d, error1d, iobj)
    print('Saved plot in file ', plotname)
    spec1dhdu.close()
    return

def main():
    if len(sys.argv) >= 2:
      fname = sys.argv[1]
    else:
      fname = input('Filename of Binospec reduced 1-d spectra, eg obj_abs_slits_extr.fits: ').strip()
    if len(sys.argv) >= 3:    
      iobj = int(sys.argv[2])
    else:
      iobj = int(input('Spectrum row number to plot, 0-based; for Bino longslit on side B, enter 1: ').strip())
    plot_spec_from_file(fname, iobj)
    return       


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()

