import matplotlib.pyplot as plt
import numpy as np

from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal

from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.constants import c
c_Angs = c*1e10

from astropy.io import fits

from collections import OrderedDict

#from .sqbase import datadir

## Setting a few things up...
softening_parameter = np.array([1.4,0.9,1.2,1.8,7.4])*1e-10

## Adjust accordingly... ;-) 
data_path = '/cos_pc19a_npr/data/WISE/LCs/'
data_filename = 'dr3_wise_lc_metrics_all_qso.fits'
data_in = data_path+data_filename

f = fits.open(data_in)   # open a FITS file

tbdata = f[3].data       # For the SDSS data of the unWISE QSOs, need the 3rd extension...

print(tbdata[:3])         # show the first three rows

def nmgy2abmag(f,df=None):
    ii = np.where(f>0.00001)
    mag = 99.99 + np.zeros_like(f)
    mag[ii] = 22.5 - 2.5*np.log10(f[ii])
    if df is None:
        return mag
    else:
        err = np.zeros_like(mag)
        err[ii] = 1.0857 * df[ii]/f[ii]
        return mag,err

 
def abmag2nmgy(_b,m):
    return 10**(-0.4*(m - 22.5))


def nmgy2asinhmag(_b,f,df=None):
    b = softening_parameter['ugriz'.find(_b)]
    mag = -1.0857*(np.arcsinh(1e-9*f/(2*b)) + np.log(b))
    if df is None:
        return mag
    else:
        err = 1.0857 * 1e-9*df/(2*b) / np.sqrt(1+((1e-9*f)/(2*b))**2)
    return mag,err


def asinhmag2nmgy(_b,m):
    b = softening_parameter['ugriz'.find(_b)]
    return 2*b*np.sinh(m/(-1.0857) - np.log(b)) / 1e-9



print(tbdata['u_flux'])    # The u-band flux

uflux = tbdata['u_flux']
gflux = tbdata['g_flux']
rflux = tbdata['r_flux']
iflux = tbdata['i_flux']
zflux = tbdata['z_flux']

#print(uflux)
#print(type(uflux), len(uflux))


#magAB = np.clip(nmgy2abmag(self.b,f_nmgy),0,magLim)
#magAB = np.clip(nmgy2abmag(self.b,f_nmgy),0,magLim)
#
#vegaMag = nmgy2abmag(self.band,f_nmgy) - self.vegaConv

u_band_mag = nmgy2abmag(uflux)
g_band_mag = nmgy2abmag(gflux)
r_band_mag = nmgy2abmag(rflux)
i_band_mag = nmgy2abmag(iflux)
z_band_mag = nmgy2abmag(zflux)

## Regular SDSS Colors!!
u_minus_g = u_band_mag  - g_band_mag
g_minus_r = g_band_mag  - r_band_mag
r_minus_i = r_band_mag  - i_band_mag
i_minus_z = i_band_mag  - z_band_mag

g_mins_i = g_band_mag  - i_band_mag


'''
Actual plotting...

'''

## from Fig 1 of Ross et al. (2012)
ug_min = -0.95   
ug_max = 4.05
gr_min = -0.75
gr_max = 2.1

fig, axs = plt.subplots(ncols=2, sharey=True, figsize=(7, 4))
fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)

ax = axs[0]
hb = ax.hexbin(u_minus_g, g_minus_r, gridsize=5000, cmap='inferno')
ax.axis([ug_min, ug_max, gr_min, gr_max])
ax.set_title("(u-g) vs. (g-r)")
cb = fig.colorbar(hb, ax=ax)
cb.set_label('counts')

ax = axs[1]
hb = ax.hexbin(u_minus_g, g_minus_r, gridsize=5000, bins='log', cmap='inferno')
ax.axis([ug_min, ug_max, gr_min, gr_max])
ax.set_title("(u-g) vs. (g-r) (log scale")
cb = fig.colorbar(hb, ax=ax)
cb.set_label('log10(N)')

plt.show()
