import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt("data.txt", skiprows = 1, unpack = True)
x = 10 ** x

fig, ax = plt.subplots()
ax.semilogx(x, y, 'o', label = "data")
ax.legend(loc = 2, title = "Legend")
ax.set_xlabel(r"Normal text vs. ${\rm math\, text}$")
ax.set_ylabel(r"A B C $\alpha$ $\beta$ $\gamma$")
#ax.minorticks_on()

fig.savefig("plot.pdf")
plt.close(fig)
