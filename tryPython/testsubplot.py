import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
import brake

# Simple data to display in various forms
x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)

fig = plt.figure('logarithmic decay of different parameters with dimension')

plt.subplot(221)
plt.plot(x, y, '-ro')
plt.yscale('log')
plt.xlabel('dimension')
plt.title('relative error between computed and exact eigenvalues')
plt.grid(True)

plt.subplot(222)
plt.plot(x, y, '-kD')
plt.yscale('log')
plt.xlabel('dimension')
plt.title('maximum angle between exact eigenvector and the projection space')
plt.grid(True)

plt.subplot(223)
plt.plot(x, y, '-bs')
plt.yscale('log')
plt.xlabel('dimension')
plt.title('maximum angle between exact eigenvector and the computed eigenvector')
plt.grid(True)

plt.subplot(224)
plt.plot(x, y, '-go')
plt.yscale('log')
plt.xlabel('dimension')
plt.title('decay of singular values')
plt.grid(True)


plt.subplot_tool()
plt.show()
brake.save('erro90rDecay', ext="png", close=True, verbose=False)