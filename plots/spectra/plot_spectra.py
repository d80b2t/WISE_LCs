'''

My, NPR's attempt... :-)

'''

## J110057   165.24045 -0.88458
## J231742   349.42751  0.09309
##
##  http://astroweb.case.edu/jakub/TA/Interpolation.html
##

import numpy as np
import scipy.interpolate as interp
import matplotlib
import matplotlib.pyplot as plt
from astroML.plotting import setup_text_plots


from astroquery.sdss import SDSS
from astropy import coordinates as coords
import astropy.units as u

#Lets text in plots use latex
#setup_text_plots(usetex=True,fontsize=22)
plt.rcParams.update({'font.size': 14})
#plt.rcParams["font.weight"] = "bold"
#plt.rcParams["axes.labelweight"] = "bold"

#matplotlib.rc('font', size=18, family='serif',                  style='normal', variant='normal',                  stretch='normal', weight='heavy')

#Find a spectrum using astroquery
pos = coords.SkyCoord('165.24045d -0.88458d', frame='icrs')
xid = SDSS.query_region(pos, spectro=True,radius=2*u.arcsec)
print(xid)
j_eleven = SDSS.get_spectra(matches=xid)


pos = coords.SkyCoord('349.42751d  0.09309d', frame='icrs')
xxid = SDSS.query_region(pos, spectro=True,radius=2*u.arcsec)
print(xxid)
guo = SDSS.get_spectra(matches=xxid)

#Create a plot using the two spectra
#loop over spectra and plot
fig, ax = plt.subplots(figsize=(26.0,14.0))

ax.set_xlim(2400, 7700)

for i in np.arange(xid['ra'].size):
    print(i)
    ax.plot((10.**j_eleven[i][1].data['loglam']/(1+0.379)),  j_eleven[i][1].data['flux']) #,label=xid['instrument'][i])
#  ax.plot((10.**j_eleven[i][1].data['loglam']/(1+0.379)), (j_eleven[i][1].data['flux'] * j_eleven[i][1].data['and_mask'])) #,label=xid['instrument'][i])

   
#for i in np.arange(xxid['ra'].size):
#    print(i)
#    ax.plot((10.**     guo[i][1].data['loglam']/(1+0.321)),     guo[i][1].data['flux']) #,label=xxid['instrument'][i])
ax.plot((10.**     guo[2][1].data['loglam']/(1+0.321)),     guo[2][1].data['flux'], linewidth=1.0) 
ax.plot((10.**     guo[1][1].data['loglam']/(1+0.321)),     guo[1][1].data['flux'], linewidth=1.2)
ax.plot((10.**     guo[0][1].data['loglam']/(1+0.321)),     guo[0][1].data['flux'], linewidth=1.4)
    
ax.grid(True, linestyle='-.')
#ax.tick_params(labelcolor='k', labelsize='medium', width=3)
#ax.set_xlabel(r'Wavelength [\AA] ')
#ax.set_ylabel(r'Flux [10$^{-17}$ ergs/cm$^2$/s/\AA]', weight='bold')

ax.set_xlabel(r'Rest Wavelength / $\AA$')
ax.set_ylabel(r'Flux [10$^{-17}$ ergs/cm$^2$/s/$\AA$]')
   
ax.legend(['J110057  MJD 51908',
           'J110057  MJD 55302',
           'J231742  MJD 51816',
           'J231742  MJD 52177',
           'J231742  MJD 52200'], 
        loc='upper right',
#        linewidth = [4,4,4,4,4], 
       shadow=True, fancybox=True, fontsize=32, frameon=True)


plt.savefig("J110057_vs_Guo_temp.png", format='png')
plt.show()


