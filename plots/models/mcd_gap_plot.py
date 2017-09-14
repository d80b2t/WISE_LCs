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

# BOSS spectrum
data = ascii.read("values_boss_spectrum.dat")
speclam_boss   = data['speclam']
specflux_boss  = data['specflux']
F_lam_all_boss = data['F_lam_all']
F_lam_tot_boss = data['F_lam_tot']
rayscat        = data['rayscat']

speclamunits= [0.1, 0.1, 0.1, 0.1]
specfluxunits=[1.0, 1.0, 3.8e17, 3.8e17]

## Redshift of J110057
redshift=0.378


# Normalisation factor(s) for display purposes
normfactor    =4.5e42
normfactor_RS = rayscat[0]

# Setting up the plot
fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(18.5, 10.5))

ax1.set_title('Sharing both axes')
#ax1.set_ylim(-20,150.)
ax1.set_ylim(-10,40.)


# plot the DATA
#ax1.plot(speclam_sdss, specflux_sdss)
ax1.plot(speclam_boss/(1+redshift), specflux_boss)
#ax1.plot(speclam_Palb, specflux_Palomar_b)
#ax1.plot(speclam_Palr, specflux_Palomar_r)
# plot the MODELS
ax1.plot(speclam_boss/(1+redshift), F_lam_all_boss/normfactor, color='black', ls='solid',  linewidth=2)
ax1.plot(speclam_boss/(1+redshift), F_lam_tot_boss/normfactor, color='black', ls='solid',  linewidth=2)

ax1.plot(speclam_boss/(1+redshift), rayscat/normfactor_RS,                   color='y', ls='dashed', linewidth=2)
ax1.plot(speclam_boss/(1+redshift), ((F_lam_all_boss/normfactor)+(rayscat/normfactor_RS)),   color='r',     ls='solid', linewidth=2)

# add axes & label them
ax1.set_ylabel(r"$Normalized \ F_{\lambda} (arb \ units)$", fontsize=22)
#legend, title
ax1.legend(loc='upper right')
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
ax2.plot(speclam_boss/(1+redshift), specflux_boss/(F_lam_all_boss+rayscat), color='blue', ls='dashed', linewidth=2)
#ax2.plot(speclam, specflux/rayscat, color='blue', ls='dashed', linewidth=2)

ax2.set_ylabel(r"$F_{\lambda} / model$" , fontsize=22)
ax2.set_xlabel(r"$\lambda (nm)$", fontsize=22)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=8)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)


#output to files in appropriate formats
#USER: choose your filenames/formats
plt.savefig('mcd_gap_v3_temp.eps')
plt.savefig('mcd_gap_v3_temp.png')
plt.savefig('mcd_gap_v3_temp.pdf')
plt.show()



