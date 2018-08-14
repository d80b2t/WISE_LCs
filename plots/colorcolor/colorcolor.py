import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#import image

from matplotlib import colors as mcolors
from astropy.io import ascii


path = '/cos_pc19a_npr/data/J110057/photometry/'
data = ascii.read(path+'SDSS_DECaLSS_colorcolor.dat')


##
## Making the plot
##
## http://matplotlib.org/examples/statistics/errorbar_limits.html
## http://matplotlib.org/examples/color/named_colors.html
plt.rcParams.update({'font.size': 14})

# 17.0 and 8.0 for paper
#fig, ax = plt.subplots(figsize=(17.0, 8.0))
# 12.0 and 8.0 for job app figure
fig, ax = plt.subplots(figsize=(8.0, 8.0))

xmin = -0.5
xmax =  1.0
ymin =  0.0
ymax =  1.0

ax.set_xlim((xmin, xmax))
ax.set_ylim((ymin, ymax))
ax.tick_params('x', direction='in')
ax.tick_params('y', direction='in')

cmap=plt.get_cmap('plasma')

clrs =  ((data['MJD'] -50000))/32.

ms=128
ax.scatter((data['gMag']-data['rMag']), (data['rMag']-data['zMag']), c=clrs, marker="o", s=ms, cmap="viridis") 

ax.text(0.55, 0.25, 'MJD 51259', transform=ax.transAxes, color='midnightblue',  fontsize=24)
ax.text(0.55, 0.20, 'MJD 56555', transform=ax.transAxes, color='yellowgreen',             fontsize=24)
ax.text(0.55, 0.15, 'MJD 57580', transform=ax.transAxes, color='khaki',         fontsize=24)
ax.text(0.55, 0.10, 'MJD 57815', transform=ax.transAxes, color='gold',          fontsize=24)


#plt.legend(['51259', '56555', '57580', '57815'],           loc="bottom right", ncol=1, fontsize=24)

ax.set_xlabel(' ( g - r ) ')
ax.set_ylabel(' ( r - z ) ')

plt.savefig('J110057_colorcolor_201805temp.pdf', format='pdf')

plt.show()
