import brake
import createBrakeClassObject
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

obj = createBrakeClassObject.returnObject(0)

# Simple data to display in various forms
x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)

fig = plt.figure(figsize=(24.0, 15.0))

plt.subplot(221)
plt.plot(x, y, '-ro', label='Classical')
plt.legend(loc='lower left', shadow=True)
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


#plt.subplot_tool()
#plt.show()

brake.save(obj.output_path+'testBrake', ext="png", close=False, verbose=True)
plt.show()
