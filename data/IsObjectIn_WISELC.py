'''
http://pyastronomy.readthedocs.io/en/latest/index.html
'''
from __future__ import print_function, division

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from PyAstronomy import pyasl

import numpy as np

path = '/cos_pc19a_npr/data/WISE/LCs/'
filename = 'dr3_wise_lc_metrics_all_qso.fits'
#filename = 'dr3_wise_lc_sample.fits'
datafile = path+filename

hdulist = fits.open(datafile)
#print(hdulist)
#hdulist.info()

data_wise = hdulist[1].data 
data_wise.columns
data_wise['MAG_RANGE_W1']

data_decals = hdulist[2].data 
#data_decals.columns
#data_decals['WISE_LC_FLUX']

data_sdss = hdulist[3].data 
print(data_sdss.columns)


#name = input("What's the object's name? ")
#hd1 = input("Input object name, with gaps between HMS...")
#hd1 = "00:05:08.83239 -67:50:24.0135"
#hd1 = "11 01 03.66 +00 49 49.0"       ## in catalog; not notable
#hd1 = "11 00 57.70 -00 53 04.500"   ## in catalog; not notable
hd1 = "23 17 42.60 +00 05 35.1"


print("Coordinates of object: ", hd1, '\n')
##print("!!!  WARNING::   -00d dec wont be treated correctly  !!! \n")

# Obtain decimal representation
ra_in, dec_in = pyasl.coordsSexaToDeg(hd1)
print("Coordinates of HD 1 [deg]: %010.6f  %+09.6f" % (ra_in, dec_in))

# Convert back into sexagesimal representation
sexa = pyasl.coordsDegToSexa(ra_in, dec_in)
print("Coordinates of HD 1 [sexa]: ", sexa)

tol = 0.005
print("Matching tolerance is ", tol, " degress")

##dec_in = dec_in * -1.
print("ra_in, dec_in", ra_in, dec_in)

w = data_sdss[np.where( (data_sdss['RA']>=(ra_in-tol)) &
                        (data_sdss['RA']<=(ra_in+tol)) &
                        (data_sdss['DEC']>=(dec_in-tol)) &
                        (data_sdss['DEC']<=(dec_in+tol)) )]
print(w)

