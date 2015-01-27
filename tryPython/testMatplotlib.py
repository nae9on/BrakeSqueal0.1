
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pylab

fig = plt.figure()

def f(x, y):
    return np.sin(x) + np.cos(y)

x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)

plt.imshow(f(x, y), cmap=plt.get_cmap('jet'))

'''
def updatefig(*args):
    global x,y
    x += np.pi / 15.
    y += np.pi / 20.
    im.set_array(f(x,y))
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
'''
plt.ion()
plt.show()
plt.hold(True)
plt.pause(0.0001)