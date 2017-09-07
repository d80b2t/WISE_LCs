#! /usr/bin/env python
#mcd_gap_v2_1_b1.py

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
####
#
# Author: K. E. Saavik Ford

from pylab import *

import sys, os, time, string, math, subprocess
import numpy as np
import matplotlib.pyplot as plt

#SI units
Msun=1.99e30 #kg per solar mass
Rsun=6.95e8 #meters per solar radius
G=6.67e-11
c=3e8
sigma_SB=5.7e-8 #stefan-boltzmann const
yr=3.15e7 #seconds per year
pc=3.086e16 #meters per parsec
AU=1.496e11 #meters per AU
h=6.626e-34 #planck const
kB=1.38e-23 #boltzmann const
m_p=1.67e-27 #mass of proton
sigma_T=6.65e-29 #Thomson xsec
PI=3.1415926

def find_B_lambda(Temp, lam):
    #Planck function
    #BB intensity=2hc^2/lam^5 * 1/(exp(hc/lamkT)-1)
    I=(2*h*c**2/pow(lam,5))/(exp(h*c/(lam*kB*Temp))-1)

    return I

def find_Temp(epsilon,M_SMBH,dotm_edd,radius,r_in,model):
    #find temp as a fn of radius
    #
    #Options:
    #Sirko & Goodman 2003 (SG): T_SG
    #Zero torque at ISCO (ZT): T_ZT  *****DEFAULT
    #Non-zero torque at ISCO (NZT): T_NZT
    #
    #Sirko & Goodman 2003 is similar to non-zero torque at the ISCO
    #---but differs by factor to approximate spectral hardening
    #assume sigma_SB T^4=(3/8pi) dotM Omega^2
    #use Omega=sqrt(GM/r^3); put r in units r_g; dotM in units dotM_edd
    prefactor=pow(((3.0/2.0)*(c**5)*m_p/(epsilon*G*M_SMBH*Msun*sigma_T*sigma_SB)), 0.25)
    T_SG=prefactor*pow(dotm_edd, 0.25)*pow(radius,-0.75)
    
    #Non-zero-torque at ISCO temp profile:
    #include factor f, O(1), which approximates spectral hardening from pure BB
    #set f=2
    #prefactor is otherwise identical to SG
    #See e.g. Gammie 1999; Krolik & Agol 2000; Narayan et al. 1997; Ashfordi & Paczynski 2003
    f=2.0
    T_NZT=f*T_SG
    
    #Zero torque at the ISCO temp profile:
    #see e.g. Zimmerman et al. 2005
    #To implement:
    if (radius >= r_in):
        T_ZT=T_NZT*pow((1.0-pow(radius/r_in, -0.5)),0.25)
    else:
        T_ZT=0.0

    if (model=='SG'):
        T=T_SG
    elif (model=='ZT'):
        T=T_ZT
    elif (model=='NZT'):
        T=T_NZT
    else:
        print 'WARNING: temperature profile model not recognized, using zero-torque model'
        T=T_ZT
        
    return T

def find_area(r1,r2,r_g):
    #find area of annulus
    #2pi R deltaR
    area=2*PI*r1*r_g*(r2*r_g-r1*r_g)

    return area

if __name__ == "__main__":
    #F_lam=sum over all blackbodies, assuming T(r), annuli of area=2piRdeltaR
    #integrating over solid angle=pi
    #BB intensity=2hc^2/lam^5 * 1/(exp(hc/lamkT)-1)
    #log_lam, lam=wavelength range of interest
    
    #BEGIN PHYSICAL INPUTS:
    #M_SMBH=mass of supermassive black hole in units of solar masses
    M_SMBH=3.0e8
    #dotm_edd=accretion rate in units of Eddington accretion
    dotm_edd=0.01
    #R_alt=outermost radius of altered disk
    #units of r_g of SMBH
    R_alt=[50.0, 100.0, 150.0]
    #Choose unperturbed temperature profile model
    #see also 'find_Temp' for more details
    #options are:
    #   SG=Sirko & Goodman 2003
    #   ZT=zero-torque at ISCO model
    #   NZT=non-zero-torque at ISCO model
    #   If none of the above, default is ZT with warning
    tempmod='NZT'
    #What fraction of flux, compared to unperturbed disk, is left in altered disk region?
    # !!!eventually would like to add more physical models!!!
    # f_depress=1.0 is unperturbed disk; f_depress=0.0 is no disk in altered region
    f_depress=0.01
    #
    #USER can, but probably should not, alter the following physical
    #variables without thinking really hard:
    #
    #inner and outer disk radii in units of r_g of SMBH, unperturbed disk
    #(inner radius depends on spin, connect to epsilon later)!!!
    radius_in=6.0
    radius_out=1.0e4
    #epsilon=radiative efficiency, depends on spin of SMBH, assume 0.1
    epsilon=0.1
    #END PHYSICAL INPUTS
    
    #BEGIN DISPLAY INPUTS:
    #divide display box for graphs
    #params for axes
    #format=left, bottom, width, height
    rect1=0.1,0.1,0.75,0.75
    
    #make figure
    fig1=plt.figure(1)
    #add axes & label them
    ax1=fig1.add_axes(rect1)
    ax1.set_ylabel(r"$log \ (\lambda F_{\lambda} (arb \ units))$")
    ax1.set_xlabel(r"$log \ (\lambda (m))$")
    #Title is set for variable R_alt
    ax1.set_title(r"$Thermal \ Continuum; \ R_{alt}=%.1f, \ %.1f, \ %.1f \ r_{g}$" %(R_alt[0],R_alt[1],R_alt[2]))
    #choose colors for different R_alt
    R_color_cycle=['r','g','b']

    #set up range for x-axis
    #USER: choose min/max wavelength in nm
    lam_min=10.0
    lam_max=1000.0
    #do arithmetic
    log_lam_min=log10(lam_min*1.0e-9)
    log_lam_max=log10(lam_max*1.0e-9)
    log_lam=np.arange(log_lam_min, log_lam_max, 0.01)
    lam=pow(10,log_lam)
    #USER: choose range for y-axis
    plt.ylim(37.0,40.0)
    #USER: vertical lines for bandpass (optional)
    #USER: choose 'y' to plot, set source redshift, min/max bandpass wavelength in nm
    plot_bandpass='y'
    redshift=0.378
    band_min=380.0
    band_max=750.0
    if (plot_bandpass=='y'):
        band_min_line=log10(band_min*1.0e-9/(1.0+redshift))
        band_max_line=log10(band_max*1.0e-9/(1.0+redshift))
        plt.axvline(x=band_min_line, color='k', linestyle='--')
        plt.axvline(x=band_max_line, color='k', linestyle='--')
    #END DISPLAY INPUTS
    
    #compute r_g for SMBH:
    r_g_SMBH=G*M_SMBH*Msun/c**2
    #set up y-axis variables
    #divvy up disk radii
    log_radius=np.arange(log10(radius_in),log10(radius_out),0.01)
    radius=pow(10,log_radius)
    #make an iterable for multiple R_alt & loop
    flux_iter_R=[]
    for j in range(len(R_alt)):
        #initialize final spectrum arrays
        F_lam_tot=np.zeros(len(lam))
        F_lam_resid=np.zeros(len(lam))
        #compute temp, area, emitted spectrum at each radius; then sum
        for i in range((len(radius)-1)):
            Temp=find_Temp(epsilon,M_SMBH,dotm_edd,radius[i],radius_in,tempmod)
            B_lambda=find_B_lambda(Temp, lam)
            Area=find_area(radius[i],radius[i+1],r_g_SMBH)
            #flux is pi*A*B_lambda per annulus, pi from integrating over solid angle
            F_lam_ann=PI*Area*B_lambda
            if (radius[i] >= R_alt[j]):
                #only add to total if outside R_alt
                F_lam_tot=F_lam_tot+F_lam_ann
            else:
                #if interior to gap, cut flux
                F_lam_tot=F_lam_tot+f_depress*F_lam_ann
                #add missing frac to resid to recover unperturbed spec
                F_lam_resid=F_lam_resid+(1-f_depress)*F_lam_ann
        #take log of total spectrum for plotting, dump to iterables for plotting
        log_F_lam=log10(F_lam_tot)
        log_lam_F_lam=log10(F_lam_tot*lam)
        F_lam_all=F_lam_tot+F_lam_resid
        log_lam_F_lam_all=log10(F_lam_all*lam)
        flux_iter_R.append(log_lam_F_lam)

    for j in range(len(R_alt)):
        ax1.plot(log_lam, flux_iter_R[j], color=R_color_cycle[j], ls='dashed', linewidth=2)
    ax1.plot(log_lam, log_lam_F_lam_all, color='black', ls='solid', linewidth=2)

    #output to files in appropriate formats
    #USER: choose your filenames/formats
    savefig('mcd_gap_v2_1_b1.eps')
    savefig('mcd_gap_v2_1_b1.png')
    show()
