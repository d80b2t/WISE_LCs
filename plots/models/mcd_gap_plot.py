#######
#
#Program purpose:
# Plot spectrum of multi-color AGN disk as function of wavelength; allow for 'other than thin disk' regions
#  ---Plot up to 3 'altered' spectra and compare to unperturbed disk,
#  ---alter by artificially depressing flux interior to R_alt
####
#
# Author: K. E. Saavik Ford, modifications by N. P. Ross

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

from astropy.io import ascii


specfilenames=['sdss_spectrum.dat', 'boss_spectrum.dat', 'w1100m0052_b.flam.dat', 'w1100m0052_r.flam.dat']

# SDSS spectrum
data = ascii.read("values_sdss_spectrum.dat")
speclam_sdss   = data['speclam']
specflux_sdss  = data['specflux']
F_lam_all_sdss = data['F_lam_all']
F_lam_tot_sdss = data['F_lam_tot']

# BOSS spectrum
data = ascii.read("values_boss_spectrum.dat")
speclam_boss   = data['speclam']
specflux_boss  = data['specflux']
F_lam_all_boss = data['F_lam_all']
F_lam_tot_boss = data['F_lam_tot']
rayscat        = data['rayscat']

# Palomar Blue
data = ascii.read("values_w1100m0052_b.flam.dat")
speclam_Palb   = data['speclam']
specflux_Palb  = data['specflux']
F_lam_all_Palb = data['F_lam_all']
F_lam_tot_Palb = data['F_lam_tot']

# Palomar Red
data = ascii.read("values_w1100m0052_r.flam.dat")
speclam_Palr   = data['speclam']
specflux_Palr  = data['specflux']
F_lam_all_Palr = data['F_lam_all']
F_lam_tot_Palr = data['F_lam_tot']


speclamunits= [0.1, 0.1, 0.1, 0.1]
specfluxunits=[1.0, 1.0, 3.8e17, 3.8e17]

specflux_Palb = specflux_Palb*3.8e17
specflux_Palr = specflux_Palr*3.8e17
 

## Redshift of J110057
redshift=0.378

#cm = plt.get_cmap('gist_rainbow')

# Normalisation factor(s) for display purposes
normfactor    =4.5e42
normfactor_RS = rayscat[0]

# Setting up the plot
fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(18.5, 10.5))

#ax1.set_title('Sharing both axes')
ax1.set_ylim(-5,120.)

color_map = plt.cm.Spectral_r

# plot the DATA
#ax1.plot(speclam_sdss, specflux_sdss)
ax1.plot(speclam_sdss/(1+redshift),    specflux_sdss, linewidth=1)
ax1.plot(speclam_boss/(1+redshift),    specflux_boss, linewidth=1)
ax1.plot(speclam_Palomar/(1+redshift), specflux_Palo, linewidth=1)
ax1.plot(speclam_Palr/(1+redshift), specflux_Palr, color='r', linewidth=1)

#ax1.plot(speclam_Palb, specflux_Palomar_b)
#ax1.plot(speclam_Palr, specflux_Palomar_r)
# plot the MODELS
ax1.plot(speclam_sdss/(1+redshift), (F_lam_all_sdss/normfactor), ls='solid',  linewidth=2)
ax1.plot(speclam_boss/(1+redshift), (F_lam_tot_boss/normfactor), ls='solid',  linewidth=2)
ax1.plot(speclam_Palr/(1+redshift), (F_lam_tot_Palr/normfactor), ls='solid',  linewidth=2)
ax1.plot(speclam_Palb/(1+redshift), (F_lam_tot_Palb/normfactor), ls='solid',  linewidth=2)

ax1.plot(speclam_boss/(1+redshift), rayscat/normfactor_RS,                  ls='dashed', linewidth=2)
ax1.plot(speclam_boss/(1+redshift), ((F_lam_all_boss/normfactor)+(rayscat/normfactor_RS)),    ls='solid', linewidth=2)

# add axes & label them
ax1.set_ylabel(r"$Normalized \ F_{\lambda} (arb \ units)$", fontsize=22)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.tick_params(axis='both', which='minor', labelsize=8)

#legend, title
#ax1.legend(loc='upper right')
#plt.set_title(r"$Thermal \ Continuum; \ R_{alt}=%.1f \ r_{g}$" % (R_alt))

 
#USER: vertical lines for bandpass (optional)
#USER: choose 'y' to plot, set source redshift, min/max bandpass wavelength in nm
plot_bandpass='n'
redshift=0.378
band_min=380.0
band_max=750.0
 
#plot bandpass lines, if selected
if (plot_bandpass=='y'):
    #band_min_line=log10(band_min*1.0e-9/(1.0+redshift))
    #band_max_line=log10(band_max*1.0e-9/(1.0+redshift))
    band_min_line=(band_min/(1.0+redshift))
    band_max_line=(band_max/(1.0+redshift))
    ax1.axvline(x=band_min_line, color='k', linestyle='--')
    ax1.axvline(x=band_max_line, color='k', linestyle='--')
    
# Bottom panel plot
ax2.plot(speclam_sdss/(1+redshift), specflux_sdss/(F_lam_all_sdss/normfactor),            ls='dashed', linewidth=1.2)
ax2.plot(speclam_boss/(1+redshift), specflux_boss/((F_lam_all_boss+rayscat)/normfactor),  ls='dashed', linewidth=1.2)
ax2.plot(speclam_Palb/(1+redshift), specflux_Palb/(F_lam_all_Palb/normfactor), color='r', ls='dashed', linewidth=1.2)
ax2.plot(speclam_Palr/(1+redshift), specflux_Palr/(F_lam_all_Palr/normfactor), color='r', ls='dashed', linewidth=1.2)
#ax2.plot(speclam_Palb/(1+redshift), specflux_Palb/(F_lam_tot_Palb/normfactor), color='maroon', ls='dashed', linewidth=1.2)
#ax2.plot(speclam_Palr/(1+redshift), specflux_Palr/(F_lam_tot_Palr/normfactor), color='maroon', ls='dashed', linewidth=1.2)

#ax2.plot(speclam, specflux/rayscat, color='blue', ls='dashed', linewidth=2)

ax2.set_ylim(-1,7.)
ax2.set_ylabel(r"$F_{\lambda} / model$" , fontsize=22)
ax2.set_xlabel(r"$\lambda (nm)$",         fontsize=22)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=8)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

#output to files in appropriate formats
#USER: choose your filenames/formats
#plt.savefig('mcd_gap_v3_temp.eps')
plt.savefig('mcd_gap_v3_temp.png')
plt.savefig('mcd_gap_v3_temp.pdf')
plt.show()



