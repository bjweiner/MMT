; NAME:
;   zfindspec_sp
;
; PURPOSE:
;   Find possible redshift matches for a set of spectra using a set of
;   eigen-templates.
;
; CALLING SEQUENCE:
;   result = zfindspec_sp( eigenfile, event, info, $
;     [eigendir=, npoly=, zmin=, zmax=, zguess=, pwidth=, nfind=, width=, $
;     objflux=, objinv=,  _EXTRA= ]
;
; INPUTS:
;   eigenfile  - Input FITS file with an [NPIXSTAR,NSTAR] image with
;                either templates or eigenspectra.  If a wildcard appears
;                in the file name, then the file that appears last in a sort
;                is used.
;                The header keywords COEFF0, COEFF1 are used to specify
;                the wavelength mapping in log-10 Angstroms.
;   event      - From specpro
;   info       - From specpro
;
; OPTIONAL KEYWORDS:
;   eigendir   - Directory for EIGENFILE; default to $IDLSPEC2D/templates.
;   columns    - Column numbers of the eigenspectra image to use in the
;                PCA fit; default to all columns.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to none.
;   zmin       - Minimum redshift to consider; default to no lower bound.
;   zmax       - Maximum redshift to consider; default to no upper bound.
;   zguess     - Initial guess for redshift; search for a solution about
;                this value.  If specified with PWIDTH, then ZMIN and ZMAX
;                are ignoreed.
;   pwidth     - Search width in pixels about the intial guess redshift ZGUESS.
;                If specified with ZGUESS, then ZMIN and ZMAX are
;                ignored.
;   nfind      - Keyword for ZCOMPUTE().
;   width      - Keyword for ZCOMPUTE().
;         
;   _EXTRA     - Keywords for ZCOMPUTE(), such as PSPACE, DOPLOT, DEBUG.
;
; OUTPUTS:
;   result     - Structure with redshift-fit information.  Structure
;                elements are left blank if fewer than NFIND peaks are found.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   One can specify a search domain for the redshift with ZMIN and ZMAX, or
;   with ZGUESS and PWIDTH.  If none of those parameters are set, then all
;   possible redshifts that overlap the object and star (template) are tested.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   concat_dir()
;   djs_filepath()
;   fileandpath()
;   readfits()
;   splog
;   sxpar()
;   zcompute()
;
; INTERNAL SUPPORT ROUTINES:
;   sp1d_struct()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;   Revised June 27, 2002 by mcc:
;      Added additional parameter linear_lambda= so that a vector 
;   of wavelength values can be passed to routine. This way the 
;   wavelength info isn't passed via the header (hdr) parameter.
;   Commented out the portion where the hdr is referenced and 
;   added new lines which take input linear lambda vector and 
;   switch to log-lambda values. 
;
;  2011: Adapted for use in SpecPro by D. Masters
;   Sep 2011: Change to resolution matching between template and data,
;   by J. Trump.  Instead of simply binning data to match the template,
;   when the data are lower resolution the template is now binned to
;   match the data.
;
;-----------------------------------------------------------------------------
function zfindspec, eigenfile, event, info, eigendir=eigendir, $
                columns=columns, zmin=zmin, zmax=zmax, $
                zguess=zguess, pwidth=pwidth, nfind=nfind, width=width, $
                objflux=objflux, objivar=objivar, loglam=loglam, $
                wvmin=wvmin, wvmax=wvmax, _EXTRA=EXTRA

  common com_zfind, starflux, starloglam0, stardloglam, $
    nstars, npoints, starflux_corrected, thisfile

  specpro_get_state, event, info
  sp = *(info.spec1Dptr)
  idxzoomrange = where(sp.lambda ge info.xrange[0] and sp.lambda le info.xrange[1])

  ss1d = {spec:[sp.flux[idxzoomrange]], ivar:[sp.ivar[idxzoomrange]], lambda:[sp.lambda[idxzoomrange]]}
  dloglam = (alog10(max(ss1d.lambda))-alog10(min(ss1d.lambda))) / n_elements(ss1d.lambda)

  templateflux = readfits(eigenfile, shdr, /SILENT)

  ;get the template wavelength information from the header.
  ;fix by DCM for different format
  ;we need to convert linear lambda spectral template to log lambda.
  crval1 = sxpar(shdr,'CRVAL1')
  cd = sxpar(shdr,'CDELT1')
  crpix1 = sxpar(shdr,'CRPIX1')
  dim_temp=size(templateflux)
  xtemp = findgen(dim_temp(1))
  lambda = crval1+cd*(xtemp+1-crpix1)
  npoints = n_elements(templateflux)

  lambdamax = lambda[npoints-1]
  lambdamin = lambda[0]
  stardloglam = (alog10(lambdamax)-alog10(lambdamin)) / npoints

  ; BJW - debugging, print wavelength spacing of template and data
  ; There are often NaNs here even when spectrum is zoomed in. Why?
  ; print, "zfindspec: stardloglam, dloglam: ", stardloglam, dloglam
  ; test = where(finite(ss1d.spec,/nan), icount)
  ; print, "zfindspec: pre-bin ss1d flux NaN at ", icount," elements"

  ; Force resolution of template and data to match
  ; BJW: what this does is interpolate onto the lower-resolution 
  ; grid, but linear2log does not smooth, so one could actually miss fine
  ; features. Why not interpolate onto the *higher* resolution grid?
  ; This does work, but not any better? and takes a long time to run

; original
   if stardloglam gt dloglam then begin  ;data res better than template
; BJW reverse
;   if stardloglam lt dloglam then begin  ;data res better than template
      loglam = linear2log(ss1d, binsize=stardloglam, flux=objflux,ivar=objivar)
      starloglam = lin2log(lambda, stardloglam, flux=templateflux)
      stardloglam = starloglam[1]-starloglam[0] ;change by dcm (stardloglam changed slightly, this updates it)
      ; BJW: shouldn't stardloglam be recomputed _before_ rebinning ss1d?
   endif else begin ;template res better than data
       starloglam = lin2log(lambda, dloglam, flux=templateflux)
       stardloglam = starloglam[1]-starloglam[0]
       loglam = linear2log(ss1d, binsize=dloglam, flux=objflux, ivar=objivar)
   endelse
  starloglam0 = starloglam[0]

  ; BJW - debugging
  ; print, "zfindspec: stardloglam, dloglam: ", stardloglam, dloglam
  ; test = where(finite(ss1d.spec,/nan), icount)
  ; print, "zfindspec: post-log-bin ss1d flux NaN at ", icount," elements"
  
  ;change, dcm 10/2/11
  ;width = (n_elements(starloglam)/20)
  median_width = round(150/((10^(starloglam[n_elements(starloglam)-1])-10^(starloglam[0]))/n_elements(starloglam)))
  templatecont = djs_median(templateflux[*,0], width=median_width, boundary='reflect')
  templateflux[*,0] = templateflux[*,0] - double(templatecont)

 ; DETERMINE GRID SIZE (IN LOG LAMBDA) FOR THE OBJECT SPECTRUM AND
 ; DETERMINE THE WAVELENGTH RANGE FOR THE OBJECT.

  objloglam0 = loglam[0]
  objdloglam = stardloglam      ; loglam[1] - loglam[0]
  ;should be same as stardloglam, if code works

  ;smooth the object spectrum and remove continuum 
  lambdarange = max(sp.lambda[idxzoomrange])-min(sp.lambda[idxzoomrange])
  ;median_width = round(lambdarange / 10.) 

  ;change, dcm 10/2/11
  ;pick a pixel width such that we get 1000 angstroms in each bin
  median_width = round(1000/((10^(loglam[n_elements(loglam)-1])-10^(loglam[0]))/n_elements(loglam)))
  ;if this is not enough points (low res) then set to a minimum number
  ;if width lt 10 then width=10

  ; BJW - the djs_median generates errors, why?  Could the median_width 
  ; cause issues, esp. if <1000 A of data?
  ; print, 'zfindspec: smoothing continuum: Npix, width: ',n_elements(objflux), median_width
  objcont = djs_median(objflux, width=median_width, boundary='reflect')
  ; print, 'zfindspec: median of flux, continuum: ',median(objflux), median(objcont)
  objflux = objflux - double(objcont)
  test = where(finite(objflux,/nan), icount)
  ; print, "zfindspec: cont-sub flux NaN at ", icount," elements"

; check if we need to trim the input spectrum to just a subregion.
  if n_elements(wvmin) gt 0 then wvmin = wvmin[0] else wvmin = -1
  if n_elements(wvmax) gt 0 then wvmax = wvmax[0] else wvmax = -1
  if wvmax ge 0 and wvmin ge 0 then begin
      minpix = findpix(loglam, alog10(wvmin))
      maxpix = findpix(loglam, alog10(wvmax))
      loglam = loglam[minpix:maxpix]
      objflux = objflux[minpix:maxpix]
      objivar = objivar[minpix:maxpix]
      objloglam0 = loglam[0]
  endif

;;; CHECK IF THE zmin AND zmax ARGUMENTS WERE PASSED. IF SO, THEN
;;; CONVERT THE REDSHIFT VALUES INTO PIXEL VALUES pmin AND pmax.
;;; THIS IS ONLY TRUE IF objloglam0 = temploglam0?
  if n_elements(zmin) ne 0 then $
    pmin = FLOOR( ALOG10(1.0 + zmin) / objdloglam )
  if n_elements(zmax) ne 0 then $
    pmax = CEIL( ALOG10(1.0 + zmax) / objdloglam )

;;; CHECK IF A GUESS REDSHIFT zguess WAS PASSED ALONG WITH A PIXEL
;;; WINDOW pwidth. IF SO, THEN RESET pmin AND pmax ACCORDING TO THE
;;; GUESS VALUE AND THE WINDOW.
  if N_ELEMENTS(zguess) GT 0 AND KEYWORD_SET(pwidth) then begin
    if KEYWORD_SET(width) THEN width1 = width $
    else width1 = pwidth
      pmin = FLOOR( ALOG10(1.0 + zguess) / objdloglam - 0.5*(pwidth+1+width1))
      pmax = FLOOR( ALOG10(1.0 + zguess) / objdloglam + 0.5*(pwidth+1+width1))
  endif

  ;ALTERATION BY BJW, 8/21/03 
  maxp =long((objloglam0 + objdloglam*.99*n_elements(objflux) - starloglam0)/$
            objdloglam)  

  if maxp lt pmax then begin ;limit upper redshift range to have overlap
     pmax = maxp
     print, 'resetting pmax to: ', maxp 
  endif


  if abs(objdloglam - stardloglam) GT 0.05*objdloglam THEN $
    MESSAGE, 'Template and object lambda resolution do NOT match!'

   ;----------
   ;Compute the redshift difference between the first pixel of the object
   ;spectra and the template.
   poffset = (objloglam0 - starloglam0) / objdloglam

   ;divide template flux by rms to avoid numerical errors
   templateflux = double(templateflux)
   templateflux = templateflux / sqrt(mean(templateflux^2))
   templateflux = float(templateflux)

   if pmin ge pmax then begin
      print, 'Inappropriate z range for this autocorrelation with this template.'
      return, 0
   endif

   ; BJW - debugging
   ; print,"zfindspec: poffset, pmin, pmax: ",poffset,pmin,pmax
   ; try to make a separate plot window
   wset, 0
   device, decomposed=1
   plot, loglam, objflux, xtitle='log lambda',ytitle='zfindspec flux + error'
   ; ,xrange=[3.6,4.0]
   oplot, loglam, 1.0/sqrt(objivar)
   ; ,color='00ff00'

  ; BJW - debug
  test = where(finite(objflux,/nan), icount)
  ; print, "zfindspec: zcompute flux NaN at ", icount," elements"
  ; try filtering the NaNs by setting any point with NaN flux or ivar
  ; to have zero flux and large variance. This helps a lot
  ibad = where(finite(objflux,/nan) or finite(objivar,/nan), icount)
  if icount gt 0 then begin
     objflux[ibad] = 0.0
     objivar[ibad] = 1.0e-10
  endif

   ;call zcomputespec to find the redshift
   zans = zcompute(objflux, objivar, templateflux, poffset=poffset, $
                  pmin=pmin, pmax=pmax, nfind=nfind, width=width, $
                  plottitle=plottitle, _EXTRA=EXTRA)

   ;----------
   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.  Do not modify any errors that are less than zero, since those
   ; can be used as just warning flags from the fit.

   indx = where(zans.dof GT 0, npeak)
   if (npeak GT 0) then $
    zans[indx].z = 10.^(objdloglam * zans[indx].z) - 1.
   jndx = where(zans.dof GT 0 and zans.z_err GE 0)

   if (jndx[0] NE -1) then $
    zans[jndx].z_err = $
     alog(10d) * objdloglam * zans[jndx].z_err * (1 + zans[jndx].z)

   return, zans

end ;zfindspec
 
