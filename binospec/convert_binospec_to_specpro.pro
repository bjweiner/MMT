
; Convert files output by binospec_quickreduce to the
; format expected by specpro, a GUI for examining multi-slit data
; see http://specpro.caltech.edu/specpro_formats.html
; and http://specpro.caltech.edu/

; Arguments:
;   file1d: name of the file with extracted 1d spectra, eg
;   'obj_counts_slits_extr.fits'
;   file2d: name of the file with linearized 2d slit spectra, eg
;   'obj_counts_slits_lin.fits'
;   outdir (optional): name of a directory to write the output files
;   into. The directory must already exist.

; The output will be individual files like 
; spec1d.<maskname>.<slitnumber>.<targetnumber>.fits
; spec2d.<maskname>.<slitnumber>.<targetnumber>.fits
; info.<maskname>.<slitnumber>.<targetnumber>.dat

; I find that some modifications to the stock specpro code are
; needed to make it work better, perform redshift fitting, etc.
; These mods will have to be released separately.

; Normally one will run the convert code on the Binospec pipeline 
; output files:
; 1-d file, e.g. obj_counts_slits_extr.fits
; 2-d file, e.g. obj_counts_slits_lin.fits
; Could use not-linearized 2-d files obj_counts_slits.fits, but we
; would need a 2-d map of the wavelength at each pixel.

; outdir argument specifies where to write the data. Because this
; will write 3 x (number of slits) files, it can easily write
; up to ~400-500 files per mask.

; There is very little error-trapping in this code at the moment.

; Todo: allow passing in an image file or files with WCS to cutout and 
; resample postage stamps,
; allow passing in photometry and photo-z files

; Benjamin Weiner, 25 July - 2 August 2018

; Mar 2022. Changes to reconcile with CNAW version
;   loop_create_files: long64(sxpar(hdr2d,'SLITTARG'))
;   loop_create_files: pass hdr2d into write_spec1d_file
;   write_spec1d_file: put header into primary extension
;   write_spec2d_file: put header into primary extension
; Changes to reconcile with Steve Willner's
;   write_spec1d_file: add 1 to index when creating lambda array
;   write_spec2d_file: add 1 to index when creating lambda array

; Get mask header structures from FITS tables in 1d file. This
; assumes the mask structures are always in extensions 4 and 5 of
; the 1d file.
pro get_mask_headers, spec1dfile, header_a, header_b
  header_a = mrdfits(spec1dfile, 4)
  header_b = mrdfits(spec1dfile, 5)
end

; Get image header (FITS header format) from extension iext, ie slit
; iext, of the 2d spectra file
function get_spec2d_header, spec2dfile, iext
  data1d = mrdfits(spec2dfile, iext, hdr)
  return, hdr
end

; Create files, also print a table of slit, object name, RA, Dec

pro loop_create_files, spec1dfile, spec2dfile, maskstruct, maskid=maskid,outdir=outdir
  ; if keyword_set(outdir) then cd,outdir
  if keyword_set(maskid) then begin
     maskstr = string(maskid,format='(I)')
  endif else begin
     maskstr = string(maskstruct[0].mask_id,format='(I)')
  endelse
  maskstr = strcompress(maskstr, /remove_all)

  nobj = n_elements(maskstruct)
  flux1d_all = mrdfits(spec1dfile, 1, hdr1d_all)
  error1d_all = mrdfits(spec1dfile, 2)
  ; Note that FITS extensions are 1-based but other things indexed by
  ; IDL are 0-based
  ; File for list of slits
  if keyword_set(outdir) then begin
     foutname = outdir + '/' + 'slitlist.txt'
  endif else begin
     foutname = 'slitlist.txt'
  endelse
  fnum = 3
  openw, fnum, foutname
  printf,fnum,'# slit   object            target   RA       DEC'
  ; loop over slit extensions
  for i = 1, nobj do begin
     slitstr = strcompress(string(i,format='(I03)'), /remove_all)
     data2d = mrdfits(spec2dfile,i,hdr2d)
     targnum = long64(sxpar(hdr2d,'SLITTARG'))
     targstr = strcompress(string(targnum,format='(I)'), /remove_all)
     suffix = maskstr + "." + slitstr + "." + targstr + ".fits"
     suffixdat = maskstr + "." + slitstr + "." + targstr + ".dat"
     fname1d = 'spec1d.' + suffix
     fname2d = 'spec2d.' + suffix
     stampname = 'stamps.' + suffix
     photname = 'phot.' + suffixdat
     infoname = 'info.' + suffixdat
     write_spec1d_file, fname1d, flux1d_all, error1d_all, i, hdr1d_all, outdir=outdir
     ; CNAW passes the hdr2d header into the write_spec1d_file procedure
     ; write_spec1d_file, fname1d, flux1d_all, error1d_all, i, hdr2d, outdir=outdir
     write_spec2d_file, fname2d, data2d, hdr2d, outdir=outdir
     write_info_file, infoname, maskstruct, hdr2d, i, outdir=outdir, outcatfilenum=fnum
     ; Not implemented yet
     ; write_stamps_file, stampname, maskstruct, stamp_image, outdir=outdir
     ; write_phot_file, photname, maskstruct, phot_catalog, outdir=outdir
     
  endfor
  close, fnum

end


; Write a spec1d file with 1-d arrays of flux, ivar, lambda
; do these explicitly need to be double?
; iobj is the iobj'th slit, need to convert to 0-based
pro write_spec1d_file, fname, flux1d_all, error1d_all, iobj, hdr1d_all, outdir=outdir, wmin=wmin, wmax=wmax
  if not keyword_set(wmin) then wmin=3600.0
  if not keyword_set(wmax) then wmax=10500.0
  mindata = -200.0
  maxdata = 10000.0
  ; Don't let ivar be too small/large? Likely to fail in flux units
  minivar = 1.0e-10
  maxivar = 1.0e8
  ctype1 = sxpar(hdr1d_all,'ctype1')
  cunit1 = sxpar(hdr1d_all,'cunit1')
  crpix1 = sxpar(hdr1d_all,'crpix1')
  crval1 = sxpar(hdr1d_all,'crval1')
  cdelt1 = sxpar(hdr1d_all,'cdelt1')
  cd1_1  = sxpar(hdr1d_all,'cd1_1')
  flux = double(flux1d_all[*,iobj-1])
  error = double(error1d_all[*,iobj-1])
  ivar = 1.0 / error^2
  ; ibaderr = where(error lt 1.0e-8, icount)
  ; if icount gt 0 then ivar[ibaderr] = 0.d0
  ibad = where(flux lt mindata, icount)
  if icount gt 0 then flux[ibad] = mindata
  ibad = where(flux gt maxdata, icount)
  if icount gt 0 then flux[ibad] = maxdata
  ibad = where(ivar lt minivar, icount)
  if icount gt 0 then ivar[ibad] = minivar
  ibad = where(ivar gt maxivar, icount)
  if icount gt 0 then ivar[ibad] = maxivar

  npix = n_elements(flux)
  ; lambda = (dindgen(npix) - crpix1) * cdelt1 + crval1
  ; Fix the off-by-1 error noted by Steve Willner
  lambda = (dindgen(npix)+1 - crpix1) * cdelt1 + crval1
  ; convert from nm to angstroms
  if strcompress(cunit1, /remove_all) eq 'nm' or strcompress(cunit1, /remove_all) eq 'NM' then lambda = lambda*10.d0
  iuse = where(lambda gt wmin and lambda lt wmax, icount)
  if icount lt 1 then print, "Warning: no pixels found in wavelength ", wmin,wmax," for ",fname
  outstruct = {flux: flux[iuse], ivar: ivar[iuse], lambda: lambda[iuse]}
  ; added by CNAW to write a header into the primary extension
  im = ''
  sxaddpar,hdr1d_all,'NAXIS',0
  sxdelpar,hdr1d_all,'NAXIS1'
  sxdelpar,hdr1d_all,'NAXIS2'
  if keyword_set(outdir) then begin
     ; added by CNAW
     writefits, outdir + '/' + fname, im, hdr1d_all
     mwrfits, outstruct, outdir + '/' + fname, /create
  endif else begin
     mwrfits, outstruct, fname, /create
  endelse
end

; Write a spec2d file with 2-d arrays of flux, ivar, lambda
; we won't have an ivar so set to zero
pro write_spec2d_file, fname, data2d, hdr2d, outdir=outdir, wmin=wmin, wmax=wmax
  if not keyword_set(wmin) then wmin=3600.0
  if not keyword_set(wmax) then wmax=10500.0
  mindata = -200.0
  maxdata = 10000.0
  minivar = 1.0e-10
  maxivar = 1.0e8
  ctype1 = sxpar(hdr2d,'ctype1')
  cunit1 = sxpar(hdr2d,'cunit1')
  crpix1 = sxpar(hdr2d,'crpix1')
  crval1 = sxpar(hdr2d,'crval1')
  cdelt1 = sxpar(hdr2d,'cdelt1')
  cd1_1  = sxpar(hdr2d,'cd1_1')
  flux = data2d
  ivar = data2d-data2d + 0.1
  ibad = where(flux lt mindata, icount)
  if icount gt 0 then flux[ibad] = mindata
  ibad = where(flux gt maxdata, icount)
  if icount gt 0 then flux[ibad] = maxdata
  ; ibad = where(ivar lt minivar, icount)
  ; if icount gt 0 then ivar[ibad] = minivar
  ; ibad = where(ivar gt maxivar, icount)
  ; if icount gt 0 then ivar[ibad] = maxivar

  lambda = data2d-data2d
  ; create first a 1-d array of lambda then fill it into the 2-d array
  lsize = size(lambda)
  nx = lsize(1)
  ny = lsize(2)
  ; lambda1 = (dindgen(nx) - crpix1) * cdelt1 + crval1
  ; Fix the off-by-1 error noted by Steve Willner
  lambda1 = (dindgen(nx)+1 - crpix1) * cdelt1 + crval1
  ; convert from nm to angstroms
  if strcompress(cunit1, /remove_all) eq 'nm' or strcompress(cunit1, /remove_all) eq 'NM' then lambda1 = lambda1*10.d0
  iuse = where(lambda1 gt wmin and lambda1 lt wmax, icount)
  if icount lt 1 then print, "Warning: no pixels found in wavelength ", wmin,wmax," for ",fname
  for j = 0,ny-1 do begin
     lambda[*,j] = lambda1
  endfor
  lambda = lambda[iuse,*]
  flux = flux[iuse,*]
  ivar = ivar[iuse,*]
  outstruct = {flux: flux, ivar: ivar, lambda: lambda}
  if keyword_set(outdir) then begin 
     mwrfits, outstruct, outdir + '/' + fname, /create
  endif else begin
     mwrfits, outstruct, fname, /create
  endelse
end

; Use info from mask header to write the info file. Can we 
; translate the mask info into an extraction location?
; iobj is the iobj'th slit, need to convert to 0-based
; If some of these fields are missing, things break - notably zphot -
; so set them to zero
pro write_info_file, fname, maskstruct, hdr2d, iobj, outdir=outdir, outcatfilenum=outcatfilenum
  indobj = iobj-1
  nypix = sxpar(hdr2d,'NAXIS2')
  slitra = sxpar(hdr2d,'SLITRA')
  slitdec = sxpar(hdr2d,'SLITDEC')
  ; Slit dimensions are probably in mm in the focal plane. 
  ; ~ 0.158 or 0.167 mm/arcsec, 5.99 arcsec/mm according to bino_extract_1d_slits.pro
  ; slitheight = sxpar(hdr2d,'SLITHEIG')
  ; slitwidth = sxpar(hdr2d,'SLITWIDT')
  slitheight = maskstruct[indobj].height
  slitwidth = maskstruct[indobj].width
  slitheight_asec = slitheight * 5.99
  slitwidth_asec = slitwidth * 5.99
  offset = maskstruct[indobj].offset
  offset_asec = offset * 5.99
  ; These are guesses. Also ~ 0.24 arcsec / pixel on detector
  ; and these need to be in pixels.
  ; The position may actually be somewhere in maskstruct[indobj].bbox
  extractpos = nypix / 2.0 + offset_asec / 0.24
  extractwidth = 1.0 / 0.24
  fnum = 1
  if keyword_set(outdir) then begin 
     openw, fnum, outdir + '/' + fname
  endif else begin
     openw, fnum, fname
  endelse

  printf, fnum, "ID  ",maskstruct[indobj].object
  printf, fnum, "RA ",maskstruct[indobj].ra
  printf, fnum, "DEC ",maskstruct[indobj].dec
; need to understand height, width, offset fields in mask header
  printf, fnum, "extractpos ",extractpos
  printf, fnum, "extractwidth ",extractwidth
  printf, fnum, "slitRA ",slitra
  printf, fnum, "slitDEC ",slitdec
  printf, fnum, "slitlen ", slitheight_asec
  printf, fnum, "slitwid ", slitwidth_asec
; Assume theta is always PA of slit, mask_PA is PA of mask, could be
; different if slit is tilted
  printf, fnum, "slitPA ",maskstruct[indobj].theta
; Photo-z not implemented yet
  printf, fnum, "zphot ",0.0
  printf, fnum, "zpdf ",0.0
  printf, fnum, "zpdf-low ",0.0
  printf, fnum, "zpdf-up ",5.0
  close, fnum
  if keyword_set(outcatfilenum) then begin
     printf, outcatfilenum, iobj, maskstruct[indobj].object, maskstruct[indobj].target, maskstruct[indobj].ra, maskstruct[indobj].dec, format='(I03, 2x, A, 1x, I9, F11.5, F11.5)'
  endif
end

; Main function at end to force compiling all the subroutines

pro convert_binospec_to_specpro, spec1dfile, spec2dfile, outdir=outdir

  get_mask_headers, spec1dfile, maskheader_a, maskheader_b
  nobj_a = n_elements(maskheader_a) 
  nobj_b = n_elements(maskheader_b) 
  nobj = nobj_a + nobj_b
  print, "Found mask headers for ",nobj," objects"
  maskheader = [maskheader_a, maskheader_b]
  ; header1 = get_spec2d_header(spec2dfile,1)
  maskid = maskheader_a[0].mask_id
  loop_create_files, spec1dfile, spec2dfile, maskheader, maskid=maskid, outdir=outdir

end


  



