__author__ = 'mjohnpayne'


import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure

np.random.seed(5)
x = range(4,61,2)
y = range(3,60,2)
p =  plt.plot(x, y, "o")
plt.axvline(40, color='r')
#plt.vlines(40,100,250)
plt.show()

