from ..tryPython import brake
import numpy as np
import matplotlib.pyplot as plt

xDim = np.array([10, 20, 30, 40, 50, 60, 70])
yDelta = np.array([1.66880306e-01, 1.94454583e-03, 5.43643504e-05, 2.05690946e-05, 2.24382233e-06, 1.07456482e-06, 3.39759376e-08])

plt.figure(1)
plt.xlabel('dimension')
plt.ylabel('relative error')
plt.title('logarithmic decay of relative error with dimension')
plt.plot(xDim, yDelta, '-go')
plt.yscale('log')
plt.grid(True)
plt.show()