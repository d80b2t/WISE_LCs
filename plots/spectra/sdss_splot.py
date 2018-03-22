'''
Example script for looking at BOSS spectra and redshift fits via Python.
Copyright Adam S. Bolton, Oct. 2009.
Licensed for free and unencumbered use in public domain.
No warranty express or implied.
'''

## Imports:
import numpy as np
import pyfits as pf
import matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt

mpl.use('TkAgg')
mpl.interactive(True)


## Reading in the DATA
topdir = '/cos_pc19a_npr/data/J110057/spectra/'

## SDSS spectrum
plate = '0277'
mjd = '51908'
spfile_sdss = topdir + '/spPlate-' + plate + '-' + mjd + '.fits'
zbfile_sdss = topdir + '/spZbest-' + plate + '-' + mjd + '.fits'

## Setting things up for the SDSS spectrum
hdulist = pf.open(spfile_sdss)
c0 = hdulist[0].header['coeff0']
c1 = hdulist[0].header['coeff1']
npix = hdulist[0].header['naxis1']
wave_sdss = 10.**(c0 + c1 * np.arange(npix))

bunit_sdss = hdulist[0].header['bunit']
flux_sdss = hdulist[0].data
ivar_sdss = hdulist[1].data
hdulist.close()
hdulist = 0

hdulist = pf.open(zbfile_sdss)
synflux_sdss = hdulist[2].data
zstruc_sdss = hdulist[1].data
hdulist.close()
hdulist = 0

fiberid_sdss = 212
fiberid_sdss = fiberid_sdss - 1 ## since 0-indexed

## Now the BOSS spectrum
plate = '3836'
mjd = '55302'
spfile_boss = topdir + '/spPlate-' + plate + '-' + mjd + '.fits'
zbfile_boss = topdir + '/spZbest-' + plate + '-' + mjd + '.fits'

hdulist = pf.open(spfile_boss)
c0 = hdulist[0].header['coeff0']
c1 = hdulist[0].header['coeff1']
npix = hdulist[0].header['naxis1']
wave_boss = 10.**(c0 + c1 * np.arange(npix))
bunit = hdulist[0].header['bunit']
flux_boss = hdulist[0].data
ivar_boss = hdulist[1].data
hdulist.close()
hdulist = 0

hdulist = pf.open(zbfile_boss)
synflux_boss = hdulist[2].data
zstruc_boss = hdulist[1].data
hdulist.close()
hdulist = 0


fiberid_boss = 258
fiberid_boss = fiberid_boss - 1 ## since 0-indexed

##
## Setting up the plot...
##
matplotlib.rc('text', usetex=True)
#matplotlib.rc('font', size=18, family='serif',                  style='normal', variant='normal',                  stretch='normal', weight='heavy')
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(16.0,8.0))

## BOSS   spec-3836-55302-0258.fits
ax.plot(wave_boss,        flux_boss[fiberid_boss,:] * (ivar_boss[fiberid_boss,:] > 0), 'r')
#ax.plot(wave_boss[:4634], synflux_boss[fiberid_boss,:], 'm')

## SDSS     spec-0277-51908-0212.fits
ax.plot(wave_sdss, flux_sdss[fiberid_sdss,:] * (ivar_sdss[fiberid_sdss,:] > 0), 'k')
#ax.plot(wave_sdss, synflux_sdss[fiberid_sdss,:], 'g')

plt.xlabel(r'Angstroms')
plt.ylabel(r'Flux $\times10^{-17}$ erg/cm$^2$/s/Ang')
plt.title('SDSS:  '+ zstruc_sdss[fiberid_sdss].field('class') + ', z = ' + str(zstruc_sdss[fiberid_sdss].field('z')) + 
          '   BOSS:  '+ zstruc_boss[fiberid_boss].field('class') + ', z = ' + str(zstruc_boss[fiberid_boss].field('z')))
plt.show()

print(wave_sdss.min(), wave_sdss.max())
print(wave_boss.min(), wave_boss.max())

w = wave_boss[np.where( (wave_boss > wave_sdss.min()) & (wave_boss < wave_sdss.max()))] 


