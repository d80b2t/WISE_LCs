#######
#
#Program purpose:
# Plot spectrum of multi-color AGN disk as function of wavelength; allow for 'other than thin disk' regions
#  ---Plot up to 3 'altered' spectra and compare to unperturbed disk,
#  ---alter by artificially depressing flux interior to R_alt
#
#To use:
# User serviceable parts limited to:
#   --variables listed (and documented) between 'BEGIN PHYSICAL INPUTS' and 'END PHYSICAL INPUTS'
#   --variables preceded by 'USER' in the DISPLAY INPUTS section
#   --output file names (preceded by USER) near end of program
#
# To find planned future improvements, search '!!!'
#
# No warranty express or implied. Use at your own risk. Does not prevent drowning.
#
# Updates in v2_2:
#  --plot F_lambda, instead of lambda*F_lambda
#  --linearize x-axis (lambda)
#  --plot normalized, linearized flux on y-axis
#
# Updates in v2_3:
#  --plot two temp models on same plot, for comparison
#    -USER pick: tempmod=solid black, tempmod2=dashed black
#  --plot wavelength in nm
#  --added legend
#
# Updates in v2_4:
#  --read observed spectra from files and overplot!!!implementing
#
####
#
# Author: K. E. Saavik Ford

from pylab import *
from astropy.io import ascii

import sys, os, time, string, math, subprocess
import numpy as np
import matplotlib.pyplot as plt


## SI units
Msun=1.99e30        # kg per solar mass
Rsun=6.95e8         # meters per solar radius
G=6.67e-11          # Big G in m^3 kg^-1 s^-2
c=3e8               # 
sigma_SB=5.7e-8     # stefan-boltzmann const
yr=3.15e7           # seconds per year
pc=3.086e16         # meters per parsec
AU=1.496e11         # meters per AU
h=6.626e-34         # Planck const
kB=1.38e-23         # boltzmann const
m_p=1.67e-27        # mass of proton
sigma_T=6.65e-29    # Thomson xsec
PI=3.1415926        # Yum...
m_per_nm=1.0e-9     # meters per nanometer

#BEGIN PHYSICAL INPUTS:
M_SMBH   = 3.0e8    # mass of supermassive black hole in units of solar masses
dotm_edd = 0.01     # accretion rate in units of Eddington accretion 
R_alt    = 150.0    # outermost radius of altered disk; units of r_g of SMBH
epsilon  = 0.1      # radiative efficiency, depends on spin of SMBH, assume 0.1

#compute r_g for SMBH:
r_g_SMBH=G*M_SMBH*Msun/c**2

# Choose unperturbed temperature profile model, see also 'find_Temp' for more details
# options are:
#   SG  = Sirko & Goodman 2003
#   ZT  = zero-torque at ISCO model
#   NZT = non-zero-torque at ISCO model
# If none of the above, default is ZT with warning
tempmod='ZT'
#plot 2nd temp mod
tempmod2='NZT'

# What fraction of flux, compared to unperturbed disk, is left in altered disk region?
# f_depress=1.0 is unperturbed disk; f_depress=0.0 is no disk in altered region
f_depress=0.10

# USER can, but probably should not, alter the following physical variables
#  inner and outer disk radii in units of r_g of SMBH, unperturbed disk
#  (inner radius depends on spin, connect to epsilon later)!!!
radius_in  = 6.0
radius_out = 1.0e4
print('radius_in, radius_out', radius_in, radius_out)

#END PHYSICAL INPUTS

    
# divvy up disk radii
log_radius = np.arange(log10(radius_in),log10(radius_out),0.01)
radius     = pow(10,log_radius)

## make an iterable for multiple R_alt & loop
## flux_iter_R=[]
   
    
def find_B_lambda(Temp, lam):
    # Planck function
    # BB intensity=2hc^2/lam^5 * 1/(exp(hc/lamkT)-1)
    if Temp < 0.1:
        print('Temp<0.1...', )
        I = 0.01
    else:
        I=(2*h*c**2/pow(lam,5))/(exp(h*c/(lam*kB*Temp))-1)
    return I


def find_Temp(epsilon,M_SMBH,dotm_edd,radius,r_in,model):
    # Find temp as a fn of radius
    # Options:
    #  Sirko & Goodman 2003 (SG): T_SG
    #  Zero torque at ISCO (ZT): T_ZT  *****DEFAULT
    #  Non-zero torque at ISCO (NZT): T_NZT
    #
    # Sirko & Goodman 2003 is similar to non-zero torque at the ISCO
    # but differs by factor to approximate spectral hardening
    # assume
    #   sigma_SB T^4=(3/8pi) dotM Omega^2
    # use
    #   Omega=sqrt(GM/r^3);
    # put r in units r_g; dotM in units dotM_edd
    prefactor = pow(((3.0/2.0)*(c**5)* m_p / (epsilon*G*M_SMBH*Msun*sigma_T*sigma_SB)), 0.25)
    T_SG=prefactor*pow(dotm_edd, 0.25)*pow(radius,-0.75)
    
    # Non-zero-torque at ISCO temp profile:
    #  include factor f, O(1), which approximates spectral hardening from pure BB
    #  set f=2
    #  prefactor is otherwise identical to SG
    #  See e.g. Gammie 1999; Krolik & Agol 2000; Narayan et al. 1997; Ashfordi & Paczynski 2003
    f=2.0
    T_NZT=f*T_SG
    
    # Zero torque at the ISCO temp profile:
    #  see e.g. Zimmerman et al. 2005
    #  To implement:
    if (radius >= r_in):
        T_ZT = T_NZT*pow((1.0-pow(radius/r_in, -0.5)),0.25)
    else:
        T_ZT = 0.0

    if (model=='SG'):
        T = T_SG
    elif (model=='ZT'):
        T = T_ZT
    elif (model=='NZT'):
        T = T_NZT
    else:
        print('WARNING: temperature profile model not recognized, using zero-torque model')
        T=T_ZT
        
    return T


def find_area(r1,r2,r_g):
    # find area of annulus
    # 2.pi.R.deltaR
    area=2*PI*r1*r_g*(r2*r_g-r1*r_g)

    return area


def get_spectrum(path, filename):
    # Assumes file contains ONLY wavelength and flux, in appropriate units, no headers, etc.
    # open the file for reading
    file1 = open(path+filename, 'r')

    # Read in the data line by line, split into float lists of wavelength and flux
    lam1list=[]
    flux1list=[]
    for line in file1:
        line=line.strip()
        columns=line.split()
        lam1list.append(float(columns[0]))
        flux1list.append(float(columns[1]))

    #close file
    file1.close()

    #re-cast as arrays (from lists) for manipulation
    lam1 = np.array(lam1list)
    flux1 = np.array(flux1list)

    return lam1, flux1


# Get observed spectra and do stuff with it (optional)
#
# USER: path= specify relative path to ascii spectra files
# USER: specfilenames= make comma-separated list of filenames for arbitrary number of normalized spectra
#       keep inside square brackets
# USER: speclamunits= factor to multiply speclam by to obtain units of nm (0.1 for Angstroms)
#       specify for each file (can in principal be different), keep inside sqaure brackets
# USER: specfluxunits= factor to multiply specflux by to obtain consistent units
#       specify for each file, keep inside sqaure brackets
path='/cos_pc19a_npr/programs/WISE/WISE_LCs/data/'
#specfilenames=['sdss_spectrum.dat', 'boss_spectrum.dat', 'w1100m0052_b.flam.dat', 'w1100m0052_r.flam.dat']
specfilenames='w1100m0052_r.flam.dat'

speclamunits= [0.1, 0.1, 0.1, 0.1]
specfluxunits=[1.0, 1.0, 3.8e17, 3.8e17]
speccolors=['c','m','y','y']

speclam=[]
specflux=[]

fn=specfilenames
obslam, obsflux=get_spectrum(path, fn)
speclam.append(obslam)
specflux.append(obsflux)

    
# F_lam = sum over all blackbodies, assuming T(r), annuli of area=2.pi.R.deltaR
# integrating over solid angle=pi
#    BB intensity=2hc^2/lam^5 * 1/(exp(hc/lamkT)-1)

# lam = speclam[0]/10. ## puts lam into nm
lam = speclam[0]/1e10 ## puts lam into m

#initialize final spectrum arrays
F_lam_tot     = np.zeros(len(lam))
F_lam_resid   = np.zeros(len(lam))
#setup model for comparison
F_lam_compare = np.zeros(len(lam))

print('F_lam_tot',     F_lam_tot)
print('F_lam_resid',   F_lam_resid)
print('F_lam_compare', F_lam_compare)

myfile = open('xyz_temp_v2pnt5.txt', 'w')
#compute temp, area, emitted spectrum at each radius; then sum
for i in range((len(radius)-1)):

    myfile.write("i is " + str(i) + "  radius " + str(radius[i]) + "    lam " + str(lam) + " \n")

    #std
    Temp         = find_Temp(epsilon,M_SMBH,dotm_edd,radius[i],radius_in,tempmod)
    #comparison
    Temp_compare = find_Temp(epsilon,M_SMBH,dotm_edd,radius[i],radius_in,tempmod2)
      
    # std
    B_lambda         = find_B_lambda(Temp, lam)
    # comparison 
    B_lambda_compare = find_B_lambda(Temp_compare, lam)
    Area=find_area(radius[i],radius[i+1],r_g_SMBH)
    # flux is pi*A*B_lambda per annulus, pi from integrating over solid angle
    F_lam_ann         = PI*Area*B_lambda
    F_lam_ann_compare = PI*Area*B_lambda_compare
    print('i, F_lam_tot, F_lam_ann', i, F_lam_tot, F_lam_ann)
    myfile.write("Temp is " + str(Temp) + "   B_lambda is " + str(B_lambda) + " \n")

    if (radius[i] >= R_alt):
        #only add to total if outside R_alt
        F_lam_tot     = F_lam_tot     + F_lam_ann
        F_lam_compare = F_lam_compare + F_lam_ann_compare
    else:
        # if interior to gap, cut flux
        F_lam_tot   = F_lam_tot   + f_depress*F_lam_ann
        # add missing frac to resid to recover unperturbed spec
        F_lam_resid = F_lam_resid + (1-f_depress)*F_lam_ann
        # compare--just add it in for now
        F_lam_compare=F_lam_compare+F_lam_ann_compare
#    print('i, Temp, lam', i, Temp, lam, radius[i], R_alt, F_lam_tot )
#    print('i, F_lam_tot, F_lam_ann', i, F_lam_tot, F_lam_ann)
    myfile.write("F_lam_ann   " + str(F_lam_ann)   + " \n")
    myfile.write("F_lam_resid " + str(F_lam_resid) + " \n")
    myfile.write("F_lam_tot   " + str(F_lam_tot)   + " \n\n")
myfile.close()

#take log of total spectrum for plotting, dump to iterables for plotting
log_F_lam = log10(F_lam_tot)

F_lam_all    =F_lam_tot+F_lam_resid
log_F_lam_all=log10(F_lam_all)

## Setting things up to model the Rayleigh scattering
## USER: choose min/max wavelength in nm
lam_min= 200.      
lam_max= 800.  
log_lam_min=log10(lam_min*1.0e-9)
log_lam_max=log10(lam_max*1.0e-9)
log_lam=np.arange(log_lam_min, log_lam_max, 0.01)
lam_RaySc=pow(10,log_lam)

m_per_nm = 1e-09
# Rayleigh scattering
#rayscat=23.0*pow(((lam_RaySc/m_per_nm)/320.0),4.0)
rayscat=23.0*pow(((lam/m_per_nm)/320.0),4.0)


ascii.write([speclam[0], specflux[0], F_lam_all, F_lam_tot, rayscat], 'values_'+specfilenames,
            names=['speclam', 'specflux', 'F_lam_all', 'F_lam_tot', 'rayscat'], overwrite=True)


#print("type(speclam), type(specflux), type(F_lam_tot)")
#print(type(speclam), type(specflux), type(F_lam_tot), '\n')
#print("len(speclam), len(specflux), len(F_lam_tot)")
#print(len(speclam), len(specflux), len(F_lam_tot), '\n')
#print("(speclam), (specflux), (F_lam_tot)")
#print((speclam), (specflux), (F_lam_tot), '\n')
#print(len(speclam), len(specflux), len(F_lam_tot), '\n')
