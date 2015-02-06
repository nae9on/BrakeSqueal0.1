import numpy as np
import matplotlib.pyplot as plt    

a = -2
b = 2

x = (b - a) * np.random.random(50) + a
y = (b - a) * np.random.random(50) + a
z = (b) * np.random.random(50)
print x.shape, y.shape, z.shape

fig = plt.figure()
hdl = plt.scatter(x,y,s=20,c=z,marker='o',vmin=0,vmax=2)
ax = plt.gca()
ax.set_xlim([-2,2])
ax.set_ylim([-2,2])
'''
#[(lower left x, lower left y), (upper right x, upper right y)] of the desired colorbar:
dat_coord = [(-1.5,1.5),(-0.5,1.75)]
#transform the two points from data coordinates to display coordinates:
tr1 = ax.transData.transform(dat_coord)
#create an inverse transversion from display to figure coordinates:
inv = fig.transFigure.inverted()
tr2 = inv.transform(tr1)
#left, bottom, width, height are obtained like this:
datco = [tr2[0,0], tr2[0,1], tr2[1,0]-tr2[0,0],tr2[1,1]-tr2[0,1]]
#and finally the new colorabar axes at the right position!
cbar_ax = fig.add_axes(datco)
#the rest stays the same:
clevs = [0, 1 , 2]
cb1 = plt.colorbar(hdl, cax=cbar_ax, orientation='horizontal', ticks=clevs)
'''
clevs = [0, 1 , 2]
cb1 = plt.colorbar(hdl, orientation='vertical', ticks=clevs)
plt.show()