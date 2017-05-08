'''
An NPR re-write of Aaron Meisner's IDL code that matches the
SDSS DR7Q and BOSS DR12Q to the NEOWISE(-R) data.

'''


def function my_galaxy_sample

; read in the decals dr3 specobj match catalog
  fname_decals = '/project/projectdirs/cosmo/data/legacysurvey/dr3/external/survey-dr3-specObj-dr13.fits'

  str = mrdfits(fname_decals, 1)

; read in sdss specobj catalog
  fname_sdss = '/project/projectdirs/cosmo/data/sdss/dr13/sdss/spectro/redux/specObj-dr13.fits'
  sdss = mrdfits(fname_sdss, 1)

; keep only rows of dr3 catalog that are matched to a class GALAXY 
; object

  w = where((str.objid NE -1) AND (strtrim(sdss.class,2) EQ 'GALAXY'))

  str = str[w]

; require at least four epochs of wise forced photometry, each with
; > 3 contributing exposures (in both bands ?)

  wise_lc_nobs = str.wise_lc_nobs

  wise_lc_nobs_w1 = reform(wise_lc_nobs[*, 0, *], 5, n_elements(str))
  wise_lc_nobs_w2 = reform(wise_lc_nobs[*, 1, *], 5, n_elements(str))

; require >3 contrib exposures in both bands
  n_good_epoch = total((wise_lc_nobs_w1 GT 3) AND (wise_lc_nobs_w2 GT 3), 1)

  wkeep = where(n_good_epoch GE 4) ; could try 3 at some point, might be fine

  str = str[wkeep]
  return, str

end

function my_qso_sample 

  fname = '/project/projectdirs/cosmo/data/legacysurvey/dr3/external/survey-dr3-DR12Q.fits'

; could cache this
  str = mrdfits(fname, 1)

; only objects that have a decals counterpart
  w = where(str.objid NE -1)
  str = str[w]

;;; hack -- merge in DR7Q objects

  fname_dr7 =  '/project/projectdirs/cosmo/data/legacysurvey/dr3/external/survey-dr3-DR7Q.fits'

  str_dr7 = mrdfits(fname_dr7, 1)
  w_dr7 = where(str_dr7.objid NE -1)
  str_dr7 = str_dr7[w_dr7]
  spherematch, str_dr7.ra, str_dr7.dec, str.ra, str.dec, 1.0d/3600.0d, m_dr7, $
      m_dr12
  matched = bytarr(n_elements(str_dr7))
  matched[m_dr7] = 1
  print, total(matched)
  str_dr7 = str_dr7[where(~matched)]
  str = struct_append(str, str_dr7)
;;;

; require at least four epochs of wise forced photometry, each with
; > 3 contributing exposures (in both bands ?)

  wise_lc_nobs = str.wise_lc_nobs

  wise_lc_nobs_w1 = reform(wise_lc_nobs[*, 0, *], 5, n_elements(str))
  wise_lc_nobs_w2 = reform(wise_lc_nobs[*, 1, *], 5, n_elements(str))

; require >3 contrib exposures in both bands
  n_good_epoch = total((wise_lc_nobs_w1 GT 3) AND (wise_lc_nobs_w2 GT 3), 1)

  wkeep = where(n_good_epoch GE 4) ; could try 3 at some point, might be fine

  str = str[wkeep]
  return, str
end

function qso_const_flux_chi2, x, p

  ; x values are irrelevant for the constant flux model, just exists
  ; to make mpfitfun happy

  return, (fltarr(n_elements(x)) + p[0])

end

function qso_best_const_flux, fluxes, sigmas, chi2=chi2, rchi2=rchi2

; should also calculate and return chi2/dof of best-fit const flux model

  guess = median(fluxes)
  result = mpfitfun('qso_const_flux_chi2', findgen(n_elements(fluxes)), $
                     fluxes, sigmas, [guess], /quiet, bestnorm=chi2)

  ndof = n_elements(fluxes) - 1 ; mean flux is the one free param
  rchi2 = chi2/float(ndof)

  return, result[0]
end

function qso_mean_snr, fluxes, sigmas, nbad=nbad

; nbad meant to be optional output

  ; get mean SNR of this qso averaged over the usable epochs
  if total(sigmas EQ 0) NE 0 then stop
  if total(~finite(sigmas)) NE 0 then stop

  snrs = fluxes/sigmas

  nbad = long(total(snrs LT (-3.0))) ; is this a good threshold ?

  return, mean(snrs)
end

function is_monotonic, fluxes, asc=asc, desc=desc
  ; is the lightcurve monotonically increasing/decreasing ?
  ; assumes the fluxes are input in order of ascending or descending MJD

  ; i believe that in my current setup they'll be input in order
  ; of ascending MJD

  ; asc and desc meant as optional outputs

  sind = sort(fluxes)

  asc = (total(sind EQ lindgen(n_elements(fluxes))) EQ n_elements(fluxes))
  desc = (total(sind EQ reverse(lindgen(n_elements(fluxes)))) EQ $
          n_elements(fluxes))

  monotonic = (asc OR desc)
  return, monotonic

end

function get_fluxes_sigmas, row

  n_w1 = (row.wise_lc_nobs)[*, 0]
  n_w2 = (row.wise_lc_nobs)[*, 1]

  wgood = where((n_w1 GT 3) AND (n_w2 GT 3), nwgood)
  if (nwgood LT 4) then stop

  fluxes_all = row.wise_lc_flux
  w1_fluxes_all = fluxes_all[*, 0]
  w2_fluxes_all = fluxes_all[*, 1]

  ivars_all = row.wise_lc_flux_ivar
  w1_ivars_all = ivars_all[*, 0]
  w2_ivars_all = ivars_all[*, 1]

  w1_fluxes = w1_fluxes_all[wgood]
  w2_fluxes = w2_fluxes_all[wgood]

  w1_ivars = w1_ivars_all[wgood]
  w2_ivars = w2_ivars_all[wgood]

; sanity check ivars here !!!
  if (total(w1_ivars EQ 0) NE 0) or (total(w2_ivars EQ 0) NE 0) then stop
  if (total(~finite(w1_ivars)) NE 0) or (total(~finite(w2_ivars)) NE 0) then $
      stop

  w1_sigmas = sqrt(1.0/w1_ivars)
  w2_sigmas = sqrt(1.0/w2_ivars)


  good_epoch_mask = bytarr(5)
  good_epoch_mask[wgood] = 1

  outstr = {w1_fluxes: w1_fluxes, $
            w2_fluxes: w2_fluxes, $
            w1_sigmas: w1_sigmas, $
            w2_sigmas: w2_sigmas, $
            good_epoch_mask : good_epoch_mask }

  return, outstr
end

function analyze_one_qso, row

; row is a structure obtained by extracting a single row from the
; quasar table

; figure out which epochs are present (and with the necessary n_exp)

  data = get_fluxes_sigmas(row)
  

  best_flux_w1 = qso_best_const_flux(data.w1_fluxes, data.w1_sigmas, $
      chi2=chi2_w1, rchi2=rchi2_w1)
  best_flux_w2 = qso_best_const_flux(data.w2_fluxes, data.w2_sigmas, $
      chi2=chi2_w2, rchi2=rchi2_w2)

; call qso_mean_snr -- what about cases w/ negative 'signal'
  mean_snr_w1 = qso_mean_snr(data.w1_fluxes, data.w1_sigmas, nbad=nbad_w1)
  mean_snr_w2 = qso_mean_snr(data.w2_fluxes, data.w2_sigmas, nbad=nbad_w2)

; compute pearson r for flux measurements in W1 versus W2
  pearson_r = correlate(data.w1_fluxes, data.w2_fluxes)
; return some sort of summary structure ?

; this could become problematic e.g. if all fluxes in some band are negative
  mag_range_w1 = 2.5*alog10(max(data[where(data.w1_fluxes GT 0)].w1_fluxes) / $
                            min(data[where(data.w1_fluxes GT 0)].w1_fluxes))
  mag_range_w2 = 2.5*alog10(max(data[where(data.w2_fluxes GT 0)].w2_fluxes) / $
                            min(data[where(data.w2_fluxes GT 0)].w2_fluxes))

  mono_w1 = is_monotonic(data.w1_fluxes, asc=asc_w1, desc=desc_w1)
  mono_w2 = is_monotonic(data.w2_fluxes, asc=asc_w2, desc=desc_w2)

  outstr = {mean_snr_w1 : mean_snr_w1,              $
            mean_snr_w2 : mean_snr_w2,              $
            pearson_r : pearson_r,                  $
            n_epoch : n_elements(data.w1_fluxes),   $
            nbad_w1 : nbad_w1,                      $
            nbad_w2 : nbad_w2,                      $
            mono_w1 : mono_w1,                      $
            mono_w2 : mono_w2,                      $
            best_flux_w1 : best_flux_w1,            $
            best_flux_w2 : best_flux_w2,            $
            chi2_w1 : chi2_w1,                      $
            chi2_w2 : chi2_w2,                      $
            rchi2_w1 : rchi2_w1,                    $
            rchi2_w2 : rchi2_w2,                    $
            flux_min_w1 : min(data.w1_fluxes),      $
            flux_max_w1 : max(data.w1_fluxes),      $
            flux_min_w2 : min(data.w2_fluxes),      $
            flux_max_w2 : max(data.w2_fluxes),      $
            mag_range_w1 : mag_range_w1,            $
            mag_range_w2 : mag_range_w2,            $
            good_epoch_mask : data.good_epoch_mask, $
            rising_w1 : asc_w1,                     $
            falling_w1 : desc_w1,                   $
            rising_w2 : asc_w2,                     $
            falling_w2 : desc_w2                     }

  return, outstr
end

function analyze_all_qso, cat=cat, galaxy=galaxy

; cat meant as optional output

; galaxy is kw arg, set to 1 to run SDSS galaxies rather than SDSS quasars

; loop over analyze_one_qso for each quasar, and store the outputs
; in a summary table

  if ~keyword_set(galaxy) then $
      cat = my_qso_sample() $
  else $
      cat = my_galaxy_sample()

  for i=0L, n_elements(cat)-1 do begin
      print, i+1, ' of  ' , n_elements(cat)
      str = analyze_one_qso(cat[i])
      if ~keyword_set(outstr) then outstr = replicate(str, n_elements(cat))
      outstr[i] = str
  endfor

  return, outstr
end

pro plot_one_row, row

  data = get_fluxes_sigmas(row)

  !p.multi = [0, 1, 2]
; w1 first
  y_range = (max(data.w1_fluxes) - min(data.w1_fluxes))
  ymin = ((min(data.w1_fluxes) - 0.05*y_range) < 0)
  ymax = max(data.w1_fluxes) + 0.05*y_range
  plot, lindgen(n_elements(data.w1_fluxes)), data.w1_fluxes, psym=1, $
      xrange=[-1, 5], yrange=[ymin, ymax]
  errplot, lindgen(n_elements(data.w1_fluxes)), data.w1_fluxes-data.w1_sigmas, $
      data.w1_fluxes+data.w1_sigmas

; and repeat for w2

  plot, lindgen(n_elements(data.w2_fluxes)), data.w2_fluxes, psym=1, $
      xrange=[-1, 5]
  errplot, lindgen(n_elements(data.w2_fluxes)), data.w2_fluxes-data.w2_sigmas, $
      data.w2_fluxes+data.w2_sigmas

end

; need wrapper that loops over the qso's
pro test_analyze

  str = my_qso_sample()

  blat = analyze_one_qso(str[0])
  
end

pro get_variable_sample, galaxy=galaxy

; get a small sample w/ correlation coeff > 0.9, monotonic 
; flux increase/decrease in both bands, and > 1 mag peak to peak variability
; in at least one of w1/w2

  summary = analyze_all_qso(cat=cat, galaxy=galaxy)

   ; cut to sample of interest, should be ~50 objects or so
  wgood = where(summary.mono_w1 AND summary.mono_w2 AND $
               (summary.pearson_r GT 0.90) AND (summary.nbad_w1 EQ 0) AND $
               (summary.nbad_w2 EQ 0) AND (summary.rchi2_w1 GE 7) AND $
               (summary.rchi2_w2 GE 7) AND ((summary.mag_range_w1 GE 0.5) OR $
               (summary.mag_range_w2 GE 0.5)), nwgood)

   summary = summary[wgood]
   cat = cat[wgood]
   ; determine coadd_id and pixel (x,y) for each object of interest

   x_coadd = fltarr(nwgood)
   y_coadd = fltarr(nwgood)
   coadd_id = strarr(nwgood)
   for i=0L, nwgood-1 do begin
       id = coadd_id_from_radec(cat[i].ra, cat[i].dec, xbest=xbest, ybest=ybest)
       coadd_id[i] = id
       x_coadd[i] = xbest
       y_coadd[i] = ybest
   endfor

   ; append this info to summary structure
   addstr = replicate({coadd_id: '', x_coadd: 0.0, y_coadd: 0.0}, nwgood)
   addstr.coadd_id = coadd_id
   addstr.x_coadd = x_coadd
   addstr.y_coadd = y_coadd

   summary = struct_addtags(summary, addstr)

   ; now check to see if any of the interesting objects are
   ; flagged by my bright star mask
   msk_vals = bytarr(nwgood)
   for i=0L, nwgood-1 do begin
       msk_name = $
           '/global/projecta/projectdirs/cosmo/work/wise/outputs/merge/' + $
           'neo2/fulldepth/' + strmid(summary[i].coadd_id, 0, 3) + '/' + $
           summary[i].coadd_id + '/*-msk.fits.gz'
       print, msk_name
       msk = readfits(msk_name)
       print,total(msk NE 0)/n_elements(msk)
       msk_val = msk[long(round(summary[i].x_coadd)), $
                     long(round(summary[i].y_coadd))]
       msk_vals[i] = msk_val
   endfor
   print, msk_vals

   summary = summary[where(msk_vals EQ 0)]
   cat = cat[where(msk_vals EQ 0)]

   ; read in official SDSS dr12q file

   ; spherematch against SDSS dr12q file

   ; write multi-extension output:
   ; ex = 1 -- my summary table
   ; ex = 2 -- SDSS dr12q 
   ; ex = 3 -- full decals measurements for each object from 
   ;           survey-dr3-DR12Q.fits

   outname = keyword_set(galaxy) ? 'dr3_wise_lc_sample-galaxy.fits' : $
                                   'dr3_wise_lc_sample.fits'
   if file_test(outname) then stop
   mwrfits, summary, outname
   mwrfits, cat, outname

   ; now get all of the SDSS information
   fname_spec = keyword_set(galaxy) ? '/project/projectdirs/cosmo/data/sdss/dr13/sdss/spectro/redux/specObj-dr13.fits' : '$SCRATCH/DR12Q.fits'

   spec = mrdfits(fname_spec, 1)
   spec_ra = keyword_set(galaxy) ? spec.plug_ra : spec.ra
   spec_dec = keyword_set(galaxy) ? spec.plug_dec : spec.dec
   spherematch, cat.ra, cat.dec, spec_ra, spec_dec, 1.0d/3600.0d, m_cat, $
       m_spec, dist

   for i=0L, n_elements(cat)-1 do begin
       sdss_info = spec[m_spec[where(m_cat EQ i)]]
       if ~keyword_set(sdss) then sdss = sdss_info else $
           sdss = struct_append(sdss, sdss_info)
   endfor

   mwrfits, sdss, outname
end

;dr3_wise_lc_sample.fits
pro gather_finder_cubes, fname=fname

  if ~keyword_set(fname) then fname = 'dr3_wise_lc_sample.fits'
  str = mrdfits(fname, 1)

  sidelen = 101 ; ?? not sure what the right size cutout should be
  half = sidelen/2

  basedir = '/project/projectdirs/cosmo/work/wise/outputs/merge'

; could do an indstart / nproc thing at some point
  for i=0L, n_elements(str)-1 do begin
      coadd_id = str[i].coadd_id
      x = str[i].x_coadd
      y = str[i].y_coadd
      ix = long(round(x))
      iy = long(round(y))

      ix_min = ix - half
      ix_max = ix + half
      iy_min = iy - half
      iy_max = iy + half

      ix_min_clip = (ix_min > 0)
      ix_max_clip = (ix_max < 2047)
      iy_min_clip = (iy_min > 0)
      iy_max_clip = (iy_max < 2047)

      ; initialize data cubes in both w1 and w2
      n_epoch = str[i].n_epoch ; should be either 4 or 5

      good_epochs = where(str[i].good_epoch_mask)

      cube_w1 = fltarr(sidelen, sidelen, n_epoch)
      cube_w2 = fltarr(sidelen, sidelen, n_epoch)

      ; now actually extract the cutouts
      for j=0L, n_epoch-1 do begin
          ; construct file names in w1 and w2
          this_epoch = good_epochs[j]
          fname_w1 = basedir + '/e' + string(this_epoch, format='(I03)') + $
              '/' + strmid(coadd_id, 0, 3) + '/' + coadd_id + '/unwise-' + $
              coadd_id + '-w1-img-u.fits'
          fname_w2 = repstr(fname_w1, '-w1-', '-w2-')

          ; want to read in the -n-u files so that i can throw out 
          ; regions with (n_u <= 2)
          fname_w1_n = repstr(fname_w1, '-img-u', '-n-u') + '.gz'
          fname_w2_n = repstr(fname_w1_n, '-w1-', '-w2-')

          w1 = readfits(fname_w1)
          w2 = readfits(fname_w2)
          w1_n = readfits(fname_w1_n)
          w2_n = readfits(fname_w2_n)
          w1 *= (w1_n GT 2)
          w2 *= (w2_n GT 2)

          ; grab cutouts, put them into cubes
          cut_w1 = fltarr(sidelen, sidelen)
          cut_w2 = fltarr(sidelen, sidelen)
          cut_w1[(ix_min_clip-ix_min):(sidelen-1-(ix_max-ix_max_clip)), $
                 (iy_min_clip-iy_min):(sidelen-1-(iy_max-iy_max_clip))] = $
                       w1[ix_min_clip:ix_max_clip, iy_min_clip:iy_max_clip]
          cut_w2[(ix_min_clip-ix_min):(sidelen-1-(ix_max-ix_max_clip)), $
                 (iy_min_clip-iy_min):(sidelen-1-(iy_max-iy_max_clip))] = $
                       w2[ix_min_clip:ix_max_clip, iy_min_clip:iy_max_clip]
          cube_w1[*,*,j] = cut_w1
          cube_w2[*,*,j] = cut_w2
      endfor
      writefits, 'qso_cubes/cube_' + string(i, format='(I03)') + '.fits', $
          cube_w1
      writefits, 'qso_cubes/cube_' + string(i, format='(I03)') + '.fits', $
          cube_w2, /append
  endfor

end

function calc_flux_single_exp, band

; calculate the flux in ab nanomaggies that corresponds to 
; the single-exposure detection threshold

; think it should be ~73 in w2 and ~63 in W1

  mag_lim_vega = ((band EQ 1) ? 15.3 : 14.5)
  offs_ab = ((band EQ 1) ? 2.699 : 3.339)
  flux = 10^((22.5-(mag_lim_vega + offs_ab))/2.5)

; these values should be in AB NANOMAGGIES !!

  return, flux
end

function color_variability_sample

; choose a relatively bright subsample to use in figuring out
; the histogram of (W1-W2) color changes

  w1_thresh = calc_flux_single_exp(1)
  w2_thresh = calc_flux_single_exp(2)

  fname = '/project/projectdirs/cosmo/data/legacysurvey/dr3/external/survey-dr3-DR12Q.fits'

; could cache this
  str = mrdfits(fname, 1)

;;; hack !!!
  fname_dr7 = '/project/projectdirs/cosmo/data/legacysurvey/dr3/external/survey-dr3-DR7Q.fits'

  str_dr7 = mrdfits(fname_dr7, 1)
  str = struct_append(str, str_dr7)

;;;

; this is the full-depth flux
  wise_flux = str.wise_flux

  w1_flux = reform(wise_flux[0,*], n_elements(str))
  w2_flux = reform(wise_flux[1,*], n_elements(str))

  w_bright = where((w1_flux GE w1_thresh) AND (w2_flux GE w2_thresh), n_bright)

; no zero or negative values in w_bright sample by construction, so i can
; look at the magnitudes of the bright objects without any problem

  mag_w1_bright = 22.5-2.5*alog10(w1_flux[w_bright]) ; AB
  mag_w2_bright = 22.5-2.5*alog10(w2_flux[w_bright]) ; AB

  plothist, mag_w1_bright-mag_w2_bright, bin=0.05

  outstr = str[w_bright]

  ; append a column with the full-depth color
  color = mag_w1_bright - mag_w2_bright

  addstr = replicate({color: 0.0, mag_w1: 0.0, mag_w2: 0.0}, n_bright)

  addstr.color = color
  addstr.mag_w1 = mag_w1_bright
  addstr.mag_w2 = mag_w2_bright

  outstr = struct_addtags(outstr, addstr)

  return, outstr

end

pro color_histogram, delta_color, delta_w1, delta_w2, sigmas_w1, sigmas_w2, $
                     frac_err_w1, frac_err_w2, fluxes_w1, fluxes_w2

  str = color_variability_sample()

  delta_color = [] ; in mags
  delta_w1 = [] ; in mags
  delta_w2 = [] ; in mags

  sigmas_w1 = []
  sigmas_w2 = []

  frac_err_w1 = []
  frac_err_w2 = []

  fluxes_w1 = []
  fluxes_w2 = []

  for i=0L, n_elements(str)-1 do begin
      print, i
      wise_lc_nobs = str[i].wise_lc_nobs
      wise_lc_flux = str[i].wise_lc_flux
      wise_lc_flux_ivar = str[i].wise_lc_flux_ivar

      color = str[i].color ; this is (w1-w2)
      w1_mag = str[i].mag_w1
      w2_mag = str[i].mag_w2
      for epoch=0, 4 do begin
          this_w1_flux = wise_lc_flux[epoch, 0]
          this_w2_flux = wise_lc_flux[epoch, 1]
          this_w1_nobs = wise_lc_nobs[epoch, 0]
          this_w2_nobs = wise_lc_nobs[epoch, 1]
          this_w1_flux_ivar = wise_lc_flux_ivar[epoch, 0]
          this_w2_flux_ivar = wise_lc_flux_ivar[epoch, 1]

          this_w1_sigma = sqrt(1.0/this_w1_flux_ivar)
          this_w2_sigma = sqrt(1.0/this_w2_flux_ivar)

          ; check that, for this epoch, both bands have (n_exp >= 8)
          if ((this_w1_nobs LT 8) OR (this_w2_nobs LT 8)) then continue
          ; check that, for this epoch, both fluxes are > 0 (dangerous)
          if ((this_w1_flux LE 0) OR (this_w2_flux LE 0)) then continue
          ; if not, then continue

          ; if so, append color change to delta_color list
          this_w1_mag = 22.5 - 2.5*alog10(this_w1_flux)
          this_w2_mag = 22.5 - 2.5*alog10(this_w2_flux)
          this_color = (this_w1_mag-this_w2_mag)
          this_delta_color = this_color-color
          delta_color = [delta_color, this_delta_color]
          delta_w1 = [delta_w1, this_w1_mag-w1_mag] 
          delta_w2 = [delta_w2, this_w2_mag-w2_mag]
          sigmas_w1 = [sigmas_w1, this_w1_sigma]
          sigmas_w2 = [sigmas_w2, this_w2_sigma]
          frac_err_w1 = [frac_err_w1, this_w1_sigma/this_w1_flux]
          frac_err_w2 = [frac_err_w2, this_w2_sigma/this_w2_flux]
          fluxes_w1 = [fluxes_w1, this_w1_flux]
          fluxes_w2 = [fluxes_w2, this_w2_flux]
      endfor
  endfor


  xvals = -1.0 + 0.01*lindgen(201)
  yvals = exp(-1.0*((xvals + 0.00397301)^2)/(0.058^2))

  xtitle = 'change in (W1-W2) color'
  ytitle = 'number'
  plothist,delta_color,bin=0.01,xrange=[-0.5,0.5], yrange=[0,13000], $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, /yst
  oplot,xvals,yvals*12150.0,color=djs_icolor('red')


  bitmap = tvrd(true=1)
  write_png, 'w1w2_color_variation.png', bitmap

  djs_contourpts, delta_w1, delta_w2, bin1=0.02,bin2=0.02, xrange=[-0.5,0.5], $
      yrange=[-0.5,0.5], xtitle=textoidl('\Delta') + 'W1 (mag)', $
      ytitle=textoidl('\Delta') + 'W2 (mag)', charsize=2
  oplot,[-0.5,0.5],[-0.5,0.5],color=djs_icolor('red')
  bitmap = tvrd(true=1)
  write_png, 'scatter_dW1_dW2.png', bitmap


  help, delta_color
end
