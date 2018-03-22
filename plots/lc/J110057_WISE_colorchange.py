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

   
## WISE    
WISE_W1 = ascii.read(path+'WISE_W1_LC.dat')
WISE_W2 = ascii.read(path+'WISE_W2_LC.dat')

WISE_L1bs = ascii.read(path+'J110057_l1b.tbl')

## For WISE, we adopt 2.699 and 3.339 as the conversions to AB from W1 and W2 Vega magnitudes,
WISE_W1_ABave = WISE_W1['W1_vega'] + 2.699
WISE_W2_ABave = WISE_W2['W2_vega'] + 3.339
WISE_W1_ABs   = WISE_L1bs['w1mpro'] + 2.699
WISE_W2_ABs   = WISE_L1bs['w2mpro'] + 3.339

comb_err = np.sqrt((WISE_L1bs['w1sigmpro'] **2.) + (WISE_L1bs['w2sigmpro']**2))

##
## Making the plot
##
## http://matplotlib.org/examples/statistics/errorbar_limits.html
## http://matplotlib.org/examples/color/named_colors.html
plt.rcParams.update({'font.size': 14})

# 17.0 and 8.0 for paper
#fig, ax = plt.subplots(figsize=(17.0, 8.0))
# 12.0 and 8.0 for job app figure
fig, ax = plt.subplots(figsize=(12.0, 8.0))

## Tidy up the figure
xmin = 55000            ## 51259  ## 49200
xmax = 58250
ymin = -0.35   # 21.05
ymax = 1.0    #16.65

ls = 'solid'
lw = 1.0
## WISE W1/W2
#ax.errorbar(WISE_W2['MJD'], WISE_W2_ABave, yerr=WISE_W2['W2_unc'], fmt='o', ms=ms, linestyle=ls, linewidth=lw*2.5, color='brown')
ms=8
#ax.plot(WISE_L1bs['mjd'], (WISE_W1_ABs-WISE_W2_ABs), 'o', color='brown')
ax.errorbar(WISE_L1bs['mjd'], (WISE_W1_ABs-WISE_W2_ABs), yerr=comb_err, fmt='o', ms=ms, color='brown', zorder=1)
ms=18
plt.plot(WISE_W1['MJD'], (WISE_W1_ABave-WISE_W2_ABave),  'o', markersize=18, ls=ls, color='indigo', zorder=2)



## Just WISE LCs
#xmin = 55000            ## 51259  ## 49200
#xmax = 58250
#ymin = 17.05   
#ymax = 14.70    #16.65


ax.set_xlim((xmin, xmax))
#ax.set_xlim((54000, xmax))
ax.set_ylim((ymin, ymax))
ax.tick_params('x', direction='in')
ax.tick_params('y', direction='in')
#ax.minorticks_on('x', direction='in')
#ay.minorticks_on()
#ax.get_xaxis().set_tick_params(which='both', direction='out')
#ax.get_yaxis().set_tick_params(which='both', direction='out')
#ax.grid(True)

## https://matplotlib.org/api/legend_api.html
plt.legend([
            'J1100-0053'],
#           loc="lower left", ncol=3, shadow=True, fancybox=True,
           loc="upper left", ncol=3, shadow=True, fancybox=True,
           fontsize=22, frameon=True)

plt.xlabel('MJD')
plt.ylabel(' W1 - W2 ')


##plt.show()
#fig.savefig("plot.pdf",)
#fig.savefig("plot.eps",format='eps')
plt.savefig('J110057_WISEdeltacolor_temp.png', format='png')
plt.show()

plt.close(fig)

