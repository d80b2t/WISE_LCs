##################################################################################################
#
# Program purpose:
#   Plot spectrum of multi-color AGN disk as function of wavelength;
#   allow for 'other than thin disk' regions
#     ---Plot up to 3 'altered' spectra and compare to unperturbed disk,
#     ---alter by artificially depressing flux interior to R_alt
#
##################################################################################################
#
# Author: K. E. Saavik Ford, modifications by N. P. Ross


import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors


from astropy.io import ascii
from astropy.table import Table, Column


specfilenames=['sdss_spectrum.dat', 'boss_spectrum.dat', 'w1100m0052_b.flam.dat', 'w1100m0052_r.flam.dat']

# SDSS spectrum
sdss = ascii.read("data/values_sdss_spectrum.dat")
speclam_sdss   = sdss['speclam']
specflux_sdss  = sdss['specflux']
F_lam_all_sdss = sdss['F_lam_all']
F_lam_tot_sdss = sdss['F_lam_tot']

# BOSS spectrum
boss = ascii.read("data/values_boss_spectrum.dat")
speclam_boss   = boss['speclam']
specflux_boss  = boss['specflux']
F_lam_all_boss = boss['F_lam_all']
F_lam_tot_boss = boss['F_lam_tot']
rayscat        = boss['rayscat']

# Palomar Blue
palb = ascii.read("data/values_w1100m0052_b.flam.dat")
speclam_palb   = palb['speclam']
specflux_palb  = palb['specflux']
F_lam_all_palb = palb['F_lam_all']
F_lam_tot_palb = palb['F_lam_tot']

# Palomar Red
data = ascii.read("data/values_w1100m0052_r.flam.dat")
speclam_palr   = data['speclam']
specflux_palr  = data['specflux']
F_lam_all_palr = data['F_lam_all']
F_lam_tot_palr = data['F_lam_tot']


speclamunits  = [0.1, 0.1, 0.1, 0.1]
specfluxunits = [1.0, 1.0, 3.8e17, 3.8e17]

specflux_palb = specflux_palb*3.8e17
specflux_palr = specflux_palr*3.8e17
 

## Redshift of J110057
redshift=0.378

#cm = plt.get_cmap('gist_rainbow')

# Normalisation factor(s) for display purposes
normfactor    = 4.5e42
normfactor_RS = rayscat[0]


# Setting up the plot
xplotsize = 18.5
yplotsize = 10.5
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(xplotsize, yplotsize))
print('Size of the plot:: ', xplotsize, yplotsize)

#ax1.set_title('Sharing both axes')
ax1.set_xlim(2250,7350.)
ax1.set_ylim(-10,125.)

color_map = plt.cm.Spectral_r

##  T O P   P A N E L
##  plot the DATA
lw = 1.4
ax1.plot(speclam_sdss/(1+redshift),    specflux_sdss, color='b', linewidth=lw)
ax1.plot(speclam_boss/(1+redshift),    specflux_boss, color='r', linewidth=lw)
ax1.plot(speclam_palr/(1+redshift),    specflux_palr, color='k', linewidth=lw)
ax1.plot(speclam_palb/(1+redshift),    specflux_palb, color='k', linewidth=lw)

# plot the MODELS
lw = 3.4
ax1.plot(speclam_sdss/(1+redshift), (F_lam_all_sdss/normfactor*3.), ls='solid',   color='b', alpha=0.5, linewidth=lw)
#ax1.plot(speclam_boss/(1+redshift), (F_lam_tot_boss/normfactor*2), ls='solid',   color='r', alpha=0.75,  linewidth=lw)
ax1.plot(speclam_palr/(1+redshift), (F_lam_tot_palr/normfactor*1.2), ls='solid',  color='k', alpha=0.75, linewidth=lw)
ax1.plot(speclam_palb/(1+redshift), (F_lam_tot_palb/normfactor*1.2), ls='solid',  color='k', alpha=0.95, linewidth=lw)

ax1.plot(speclam_boss/(1+redshift),  rayscat/normfactor_RS,      color='r', alpha=0.75,       ls='dashed', linewidth=lw)
ax1.plot(speclam_boss/(1+redshift), ((F_lam_all_boss/normfactor)+(rayscat/normfactor_RS)),
         color='r', alpha=0.75,   ls='solid', linewidth=lw)

# add axes & label them
#ax1.set_ylabel(r"Normalized\nF$_{\lambda}$ (arb. units)", fontsize=20, multialignment='center')
ax1.set_ylabel("Normalized f$_{\lambda}$\n(arb. units)", fontsize=20, multialignment='center')
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.tick_params(axis='both', which='minor', labelsize=8)

legend_prop = {'weight':'bold'}
lw=2

ax1.legend(['SDSS - UT 2000 December 30',
            'BOSS - UT 2010 April 16',
            'Palomar/DBSP - UT 2017 February 25',
            'model', 'model', 'model+Rayleigh'], 
            loc='upper right', fontsize=20,
            ncol=2, 
#            shadow=True,
            fancybox=True
            #frameon=True
            )


#USER: vertical lines for bandpass (optional)
#USER: choose 'y' to plot, set source redshift, min/max bandpass wavelength in nm
plot_bandpass='n'
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

        
## MIDDLE  PANEL plot
## Dividing the spectra through by the model
lw=1.4
ax2.plot(speclam_sdss/(1+redshift), specflux_sdss/( F_lam_all_sdss/normfactor*3),          color='b', ls='solid', linewidth=lw)
ax2.plot(speclam_boss/(1+redshift), specflux_boss/((F_lam_all_boss+rayscat)/normfactor), color='r', ls='solid', linewidth=lw)
ax2.plot(speclam_palb/(1+redshift), specflux_palb/( F_lam_all_palb/normfactor*1.2),          color='k', ls='solid', linewidth=lw)
ax2.plot(speclam_palr/(1+redshift), specflux_palr/( F_lam_all_palr/normfactor*1.2),          color='k', ls='solid', linewidth=lw)

#ax2.set_ylim(-1,6.7)
ax2.set_ylim(-0.35,3.4)
ax2.set_ylabel(r"f$_{\lambda}$  /  model" , fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=18)
ax2.tick_params(axis='both', which='minor', labelsize=10)


## BOTTOM PANEL
## RATIO OF THE DATA

# Make the BOSS spectrum the same size as the SDSS one
boss_cut = boss[(boss['speclam']>3801.) & (boss['speclam']<9189.56)]
flux_ratio_a = (sdss['specflux']/boss_cut['specflux'])

lam_wave = sdss['speclam']
## f_lam \propto -(\alpha_nu + 2)
alpha_nu=(1/3.)
scaling = 8000
pwr_law_a = (lam_wave/scaling)**((alpha_nu+2.)*(-1.0))
alpha_nu=(1/2.)
pwr_law_b = (lam_wave/scaling)**((alpha_nu+2.)*(-1.0))
alpha_nu=(1/4.)
pwr_law_c = (lam_wave/scaling)**((alpha_nu+2.)*(-1.0))
alpha_nu=(1/3.)
pwr_law_d = (lam_wave)**((alpha_nu+2.)*(-1.0))

## https://matplotlib.org/examples/color/named_colors.html
lw=1.4
dsb = 'darkslateblue' 
ax3.plot(speclam_sdss/(1+redshift), flux_ratio_a, color=dsb, ls='solid', linewidth=lw)
#ax3.plot(lam_wave/(1+redshift),     pwr_law_a,      color=dsb,  ls='solid',  alpha=0.75, linewidth=lw)
#ax3.plot(lam_wave/(1+redshift),     pwr_law_b,      color=dsb,  ls='-.',  alpha=0.75, linewidth=lw)
#ax3.plot(lam_wave/(1+redshift),     pwr_law_c,      color=dsb,  ls='dotted',  alpha=0.75, linewidth=lw)
lw=2.4
ax3.plot(lam_wave/(1+redshift),    (pwr_law_d*6e8),   color=dsb,  ls='dotted',  alpha=0.75, linewidth=lw)

ax3.set_ylabel(r"|${\Delta}$ f$_{\lambda}$|" , fontsize=24)
ax3.set_xlabel("Rest Frame $\lambda$ (nm)",    fontsize=22)
ax3.tick_params(axis='both', which='major', labelsize=18)
ax3.tick_params(axis='both', which='minor', labelsize=10)

ax3.legend(['SDSS/BOSS',
           r'f$_{\nu}\propto\nu^{1/3}$ '],
            loc='upper right', fontsize=22,
            ncol=1, 
#            shadow=True,
            fancybox=True
            #frameon=True
            )




# Make subplots close to each other and hide x-ticks for all but bottom plot.
fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

# Output to files in appropriate formats
plt.savefig('mcd_gap_v3_temp.png')
plt.savefig('mcd_gap_v3_temp.pdf')
plt.show()



