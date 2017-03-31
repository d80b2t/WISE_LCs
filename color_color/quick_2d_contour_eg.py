'''
MASSIVE H/T:
https://python4mpia.github.io/intro/quick-tour.html


Making a fancy plot from Monte-Carlo samples

Assume you have run an MCMC and you are left with two arrays X,Y of
MCMC samples of two fit parameters. You now want to use X,Y to
visualise the likelihood manifold. You can do that (a) as a simple
scatter plot or (b) in a more fancy way.

Instead of Monte-Carlo samples, you could also be faced with
distributions of any two parameters, such as effective temperature and
surface gravity of a set of stars, or redshift and magnitude of a set
of galaxies.

First, let us create some artificial toy data to mimick the output of
an MCMC algorithm in some science application:
'''

import numpy,math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Create artificial data mimicking some MCMC results.
N = 50000
X = numpy.random.normal(0.0, 1.5, N)  # Draw N samples from normal distribution
Y = numpy.random.gamma(2.0, 2.0, N)   # Draw N samples from Gamma distribution

'''
Second, let us create a simple plot by plainly plotting x vs. y. This is very easy and we can recap some of the basic Python plotting commands:
'''

# Define plot ranges at beginning, since used often later.
YRANGE = [-0.4,11.4]
XRANGE = [-6.4,6.4]

# Define figure size and formatting
fig = plt.figure(1, figsize=(7,7))
fig.subplots_adjust(left=0.10, bottom=0.09, top=0.98, right=0.98)

# Simply plot X vs. Y as data points.
plt.plot(X, Y, 'o', ms=4, alpha=0.1, color='blue')

# Set plot ranges, axes ticks and axes labels.
plt.xlim(XRANGE)                 # Set x plot range.
plt.ylim(YRANGE)                 # Set y plot range.
plt.xticks(fontsize=16)          # Set ticks x axis.
plt.yticks(fontsize=16)          # Set ticks y axis.
plt.xlabel(r'$x$', fontsize=24)  # Set label x axis.
plt.ylabel(r'$y$', fontsize=24)  # Set label y axis.

plt.savefig('plot_MCMC_samples_plain_temp.png') # Save png file.


'''

And now....

'''

fig = plt.figure(2, figsize=(7,7))
fig.subplots_adjust(hspace=0.001, wspace=0.001, left=0.10, bottom=0.095, top=0.975, right=0.98)
# gridspec enables you to assign different formats to panels in one plot.
gs = gridspec.GridSpec(2, 2, width_ratios=[1,4], height_ratios=[4,1])

plt.subplot(gs[1]) # Main panel top right contains full 2D histogram.
# Convert to 2d histogram.
Bins = 25
hist2D, xedges, yedges = numpy.histogram2d(X, Y, bins=[Bins,Bins], range=[XRANGE,YRANGE],
    normed=False)

# Plot Monte-Carlo samples as 2D histogram.
hist2D = numpy.transpose(hist2D)  # Beware: numpy switches axes, so switch back.
plt.pcolormesh(xedges, yedges, hist2D, cmap=plt.cm.gray)

# Overplot with error contours 1,2,3 sigma.
maximum    = numpy.max(hist2D)
[L1,L2,L3] = [0.5*maximum,0.25*maximum,0.125*maximum]  # Replace with a proper code!
# Use bin edges to restore extent.
extent = [xedges[0],xedges[-1], yedges[0],yedges[-1]]

## ARGH!!!! ARRGH!!!!!!
## ARGH!!!! ARRGH!!!!!!
## ARGH!!!! ARRGH!!!!!!
### *Something* isn't working quite right with this line:::

#cs = plt.contour(hist2D, extent=extent, levels=[L1,L2,L3], linestyles=['--','--','--'], colors=['orange','orange','orange'], linewidths=1)
cs = plt.contour(hist2D)

# use dictionary in order to assign your own labels to the contours.
fmtdict = {L1:r'$1\sigma$',L2:r'$2\sigma$',L3:r'$3\sigma$'}
plt.clabel(cs, fmt=fmtdict, inline=True, fontsize=20)

plt.xlim(XRANGE)
plt.ylim(YRANGE)

'''
Finally, add the two side panels showing the projected distributions of X and Y:
'''

# Bin X,Y separately. As 1D bin, can use more bins now.
S  = 101
LX = numpy.histogram(X, bins=S, range=XRANGE, normed=True)[0]
LY = numpy.histogram(Y, bins=S, range=YRANGE, normed=True)[0]
# Restore positions lost by binning.
X = XRANGE[0] + (XRANGE[1]-XRANGE[0])*numpy.array(range(0,len(LX)))/float(len(LX)-1)
Y = YRANGE[0] + (YRANGE[1]-YRANGE[0])*numpy.array(range(0,len(LY)))/float(len(LY)-1)

# bottom right panel: projected density of x.
plt.subplot(gs[3])
plt.plot(X, LX, '-', lw=3, color='black')

plt.xticks(fontsize=16)
plt.yticks([])
plt.xlabel(r'$x$', fontsize=24)
plt.ylabel(r'$\cal L$', fontsize=24)
plt.xlim(XRANGE)
plt.ylim(0.0, 1.1*numpy.max(LX))

# top left panel: projected density of y.
plt.subplot(gs[0])
plt.plot(LY, Y, '-', lw=3, color='black')

plt.yticks(fontsize=16)
plt.xticks([])
plt.xlabel(r'$\cal L$', fontsize=24)
plt.ylabel(r'$y$', fontsize=24)
plt.xlim(0.0, 1.1*numpy.max(LY))
plt.ylim(YRANGE)

plt.savefig('plot_MCMC_samples_fancy_temp.png')
plt.show()
