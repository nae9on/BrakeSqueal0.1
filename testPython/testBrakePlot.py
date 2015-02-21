import brake
from matplotlib import cm
from numpy.random import randn
import createBrakeClassObject
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

obj = createBrakeClassObject.returnObject(0)

# Simple data to display in various forms
x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)
z = np.linspace(0, 2 * np.pi, 400)
fig = plt.figure(figsize=(24.0, 15.0))
'''
fig = plt.plot(x, y, 'ko',  label='Classical')
ax = plt.gca()
plt.plot(x, 5*y, 'r+', markersize=10, fillstyle='none', label='5times')
plt.plot(x, 10*y, 'go', markersize=10, fillstyle='none', label='10times')
plt.legend(loc='lower left', shadow=True)
plt.yscale('log')
plt.xlabel('dimension')
plt.title('relative error between computed and exact eigenvalues')
plt.grid(True)
'''
fig = plt.figure()
hdl = plt.scatter(x, y, s=100, c=z, marker='o', facecolors='none', linewidth='0', vmin=0, vmax=2, label='5times')
hdl = plt.scatter(x, 5*y, s=100, c=z, marker='o', facecolors='none', vmin=0,vmax=2, label='5times')
ax = plt.gca()
ax.set_xlim([-2,2])
ax.set_ylim([-2,2])
clevs = [0, 1 , 2]
cb1 = plt.colorbar(hdl, orientation='vertical', ticks=clevs)
plt.show()
