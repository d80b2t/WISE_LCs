'''
Why doesn't matplotlib look as good as SM??!!!
http://space.mit.edu/home/turnerm/python.html
'''

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

from astropy.io import ascii

path = '/cos_pc19a_npr/data/J110057/photometry/'


##
##  P h o t o m e t r y    d a t a
##
## CRTS
CRTS = ascii.read(path+'CRTS_LC.dat')
## LINEAR
LINEAR = ascii.read(path+'LINEAR_full.dat')

## Pan-STARRS
PS_g = ascii.read(path+'PanSTARRS_LC_g.dat')
PS_r = ascii.read(path+'PanSTARRS_LC_r.dat')
PS_i = ascii.read(path+'PanSTARRS_LC_i.dat')
PS_z = ascii.read(path+'PanSTARRS_LC_z.dat')
PS_y = ascii.read(path+'PanSTARRS_LC_y.dat')

## SDSS 
SDSS_u = ascii.read(path+'SDSS_orig_u.dat')
SDSS_g = ascii.read(path+'SDSS_orig_g.dat')
SDSS_r = ascii.read(path+'SDSS_orig_r.dat')
SDSS_i = ascii.read(path+'SDSS_orig_i.dat')
SDSS_z = ascii.read(path+'SDSS_orig_z.dat')

## DECaLS
DECaLS_g = ascii.read(path+'DECaLS_DR4_g.dat')
DECaLS_r = ascii.read(path+'DECaLS_DR4_r.dat')
DECaLS_z = ascii.read(path+'DECaLS_DR4_z.dat')

## UKIDSS LAS
ULAS_Y = ascii.read(path+'ULAS_Y.dat')
ULAS_J = ascii.read(path+'ULAS_J1.dat')
ULAS_H = ascii.read(path+'ULAS_H.dat')
ULAS_K = ascii.read(path+'ULAS_K.dat')

## From good ol' Peth et al. (2011)
#Y_AB  = Y  + 0.634
#J_AB  = J  + 0.938
#H_AB  = H  + 1.379
#Ks_AB = Ks + 1.84	
#K_AB  = K  + 1.900

ULAS_Y_AB  = ULAS_Y['YAperMag3']   + 0.634
ULAS_J_AB  = ULAS_J['J1_AperMag3'] + 0.938
ULAS_H_AB  = ULAS_H['HAperMag3']   + 1.379
ULAS_Ks_AB = ULAS_K['KAperMag3']   + 1.84	
   
## WISE    
WISE_W1 = ascii.read(path+'WISE_W1_LC.dat')
WISE_W2 = ascii.read(path+'WISE_W2_LC.dat')

## For WISE, we adopt 2.699 and 3.339 as the conversions to AB from W1 and W2 Vega magnitudes,
WISE_W1_AB = WISE_W1['W1_vega'] + 2.699
WISE_W2_AB = WISE_W2['W2_vega'] + 3.339


## SPECTRO-PHOTOMETRY
sdss_spec_u = ascii.read(path+'SDSS_spectrophot_u.dat')
sdss_specu = 22.5-2.5*(math.log10(sdss_spec_u['flux_nanomaggies']))
sdss_spec_g = ascii.read(path+'SDSS_spectrophot_g.dat')
sdss_specg = 22.5-2.5*(math.log10(sdss_spec_g['flux_nanomaggies']))
sdss_spec_r = ascii.read(path+'SDSS_spectrophot_r.dat')
sdss_specr = 22.5-2.5*(math.log10(sdss_spec_r['flux_nanomaggies']))
sdss_spec_i = ascii.read(path+'SDSS_spectrophot_i.dat')
sdss_speci = 22.5-2.5*(math.log10(sdss_spec_i['flux_nanomaggies']))
sdss_spec_z = ascii.read(path+'SDSS_spectrophot_z.dat')
sdss_specz = 22.5-2.5*(math.log10(sdss_spec_z['flux_nanomaggies']))

boss_spec_u = ascii.read(path+'BOSS_spectrophot_u.dat')
boss_specu = 22.5-2.5*(math.log10(boss_spec_u['flux_nanomaggies']))
boss_spec_g = ascii.read(path+'BOSS_spectrophot_g.dat')
boss_specg = 22.5-2.5*(math.log10(boss_spec_g['flux_nanomaggies']))
boss_spec_r = ascii.read(path+'BOSS_spectrophot_r.dat')
boss_specr = 22.5-2.5*(math.log10(boss_spec_r['flux_nanomaggies']))
boss_spec_i = ascii.read(path+'BOSS_spectrophot_i.dat')
boss_speci = 22.5-2.5*(math.log10(boss_spec_i['flux_nanomaggies']))
boss_spec_z = ascii.read(path+'BOSS_spectrophot_z.dat')
boss_specz = 22.5-2.5*(math.log10(boss_spec_z['flux_nanomaggies']))


##
## Making the plot
##
## http://matplotlib.org/examples/statistics/errorbar_limits.html
## http://matplotlib.org/examples/color/named_colors.html
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(figsize=(16, 7))

ls = 'dotted'
## CRTS
ax.errorbar(CRTS['MJD'], CRTS['mag'], yerr=CRTS['mag_err'], linestyle=ls, color='black')
## LINEAR
ax.errorbar(LINEAR['MJD'], LINEAR['mag'], yerr=LINEAR['mag_err'], linestyle='dashed', color='dimgray')

## Pan-STARRS
ax.errorbar(PS_g['MJD'], PS_g['mag'], yerr=PS_g['mag_err'], fmt='d', linestyle=ls, color='olivedrab')
ax.errorbar(PS_r['MJD'], PS_r['mag'], yerr=PS_r['mag_err'], fmt='d', linestyle=ls, color='red')
ax.errorbar(PS_i['MJD'], PS_i['mag'], yerr=PS_i['mag_err'], fmt='d', linestyle=ls, color='teal')
ax.errorbar(PS_z['MJD'], PS_z['mag'], yerr=PS_z['mag_err'], fmt='d', linestyle=ls, color='lightgrey')
ax.errorbar(PS_y['MJD'], PS_y['mag'], yerr=PS_y['mag_err'], fmt='d', linestyle=ls, color='darkmagenta')

## SDSS "original photometry"
ax.errorbar(SDSS_u['MJD'], SDSS_u['mag'], yerr=SDSS_u['mag_err'], fmt='o', linestyle=ls, color='b')
ax.errorbar(SDSS_g['MJD'], SDSS_g['mag'], yerr=SDSS_g['mag_err'], fmt='o', linestyle=ls, color='olivedrab')
ax.errorbar(SDSS_r['MJD'], SDSS_r['mag'], yerr=SDSS_r['mag_err'], fmt='o', linestyle=ls, color='red')
ax.errorbar(SDSS_i['MJD'], SDSS_i['mag'], yerr=SDSS_i['mag_err'], fmt='o', linestyle=ls, color='teal')
ax.errorbar(SDSS_z['MJD'], SDSS_z['mag'], yerr=SDSS_z['mag_err'], fmt='o', linestyle=ls, color='lightgrey')

## Photometry from SDSS spectrum...
plt.plot(sdss_spec_u['mjd'], sdss_specu, marker='s', markersize=7, color="b")
plt.plot(sdss_spec_g['mjd'], sdss_specg, marker='s', markersize=7, color="olivedrab")
plt.plot(sdss_spec_r['mjd'], sdss_specr, marker='s', markersize=7, color="red")
plt.plot(sdss_spec_i['mjd'], sdss_speci, marker='s', markersize=7, color="teal")
plt.plot(sdss_spec_z['mjd'], sdss_specz, marker='s', markersize=7, color="lightgrey")
## Photometry from BOSS spectrum...
plt.plot(boss_spec_u['mjd'], boss_specu, marker='s', markersize=7, color="b")
plt.plot(boss_spec_g['mjd'], boss_specg, marker='s', markersize=7, color="olivedrab")
plt.plot(boss_spec_r['mjd'], boss_specr, marker='s', markersize=7, color="red")
plt.plot(boss_spec_i['mjd'], boss_speci, marker='s', markersize=7, color="teal")
plt.plot(boss_spec_z['mjd'], boss_specz, marker='s', markersize=7, color="lightgrey")

##DECaLS
ax.errorbar(DECaLS_g['mjd'], DECaLS_g['aper_diam2pnt0'], yerr=DECaLS_g['err_diam2pnt0'], fmt='P', linestyle=ls, color='olivedrab')
ax.errorbar(DECaLS_r['mjd'], DECaLS_r['aper_diam2pnt0'], yerr=DECaLS_r['err_diam2pnt0'], fmt='P', linestyle=ls, color='red')
ax.errorbar(DECaLS_z['mjd'], DECaLS_z['aper_diam2pnt0'], yerr=DECaLS_z['err_diam2pnt0'], fmt='P', linestyle=ls, color='lightgrey')

## UKIDSS LAS
ax.errorbar(ULAS_Y['MJD'], ULAS_Y_AB,  yerr=ULAS_Y['YAperMag3Err'],   fmt='h', linestyle=ls, color='y', ms=10)
ax.errorbar(ULAS_J['MJD'], ULAS_J_AB,  yerr=ULAS_J['J1_AperMag3Err'], fmt='h', linestyle=ls, color='goldenrod', ms=10)
ax.errorbar(ULAS_H['MJD'], ULAS_H_AB,  yerr=ULAS_H['HAperMag3Err'],   fmt='h', linestyle=ls, color='darkgoldenrod', ms=10)
ax.errorbar(ULAS_K['MJD'], ULAS_Ks_AB, yerr=ULAS_K['KAperMag3Err'],   fmt='h', linestyle=ls, color='orangered', ms=10)

## WISE W1/W2
ax.errorbar(WISE_W1['MJD'], WISE_W1_AB, yerr=WISE_W1['W1_unc'], fmt='o',linestyle=ls, color='indigo')
ax.errorbar(WISE_W2['MJD'], WISE_W2_AB, yerr=WISE_W2['W2_unc'], fmt='o', linestyle=ls, color='brown')

## Plotting the SPECTRA as vertical lines
plt.axvline(x=51908, linewidth=2.5)
plt.axvline(x=55302, linewidth=2.5, linestyle='dotted', color='darkgoldenrod')
plt.axvline(x=57809, linewidth=2.5, linestyle='dashed', color='slategrey')




## Tidy up the figure
xmin= 49200
ax.set_xlim((xmin, 58250))
#ax.set_xlim((54000, 58500))
ax.set_ylim((20.05, 16.70))
ax.tick_params('x', direction='in')
ax.tick_params('y', direction='in')
#ax.minorticks_on('x', direction='in')
#ay.minorticks_on()
#ax.get_xaxis().set_tick_params(which='both', direction='out')
#ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.grid(True)

## https://matplotlib.org/api/legend_api.html
plt.legend([
    'specphot u', 'specphot g', 'specphot r', 'specphot i', 'specphot z',
    'specphot u', 'specphot g', 'specphot r', 'specphot i', 'specphot z',
    'SDSS spectrum','BOSS spectrum', 'Palomar spectrum',
            'CRTS', 'LINEAR',
            'PanSTARRS g', 'PanSTARRS r', 'PanSTARRS i', 'PanSTARRS z', 'PanSTARRS y',
            'SDSS u', 'SDSS g', 'SDSS r', 'SDSS i', 'SDSS z',
            'DECaLS g', 'DECaLS r', 'DECaLS z',
            'ULAS Y', 'ULAS J', 'ULAS H', 'ULAS K',
            'WISE W1', 'WISE W2'],
           loc="lower left", ncol=2, shadow=True, fancybox=True,
           fontsize=10, frameon=True)

plt.xlabel('MJD')
plt.ylabel('magnitude (AB)')
plt.show()


fig.savefig("plot.pdf")
plt.close(fig)

#savefig('foo.png', bbox_inches='tight')
#savefig('foo.pdf', bbox_inches='tight')
