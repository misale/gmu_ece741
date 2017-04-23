import numpy as np
import math
import matplotlib.pyplot as plt

def f(d):
    return  ((13*np.log10(d)) + (35*np.log10(1 + (d/90))))*(-1)

d = np.arange(0.0, 200.0, 1)

plt.plot(d, f(d),'b')
plt.xlabel('distance (meters)')
plt.ylabel('power (dB)')
plt.show()	
