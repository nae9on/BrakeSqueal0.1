
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pylab

plt.ion()

fig = plt.figure()
fig2 = plt.figure()

ax = fig.add_subplot(1,1,1) # two rows, one column, first plot
rect = matplotlib.patches.Rectangle( (1,1), width=5, height=12)
ax.add_patch(rect)
ax.autoscale_view()
ax.figure.canvas.draw()


t = np.arange(0.0, 1.0, 0.01)
s = np.sin(2*np.pi*t)

ax2 = fig2.add_axes([0.15, 0.1, 0.7, 0.3])
ax2.plot(t,s)


plt.show()

