#!/usr/bin/python
import numpy as np
import math
import matplotlib.pyplot as plt


def locationUp(radius,nla) :
	r = float(radius)
	n = float(nla)
	rate = ( 360 / (3 * (3**0.5) * math.pi * r) ) * ( 1 - ( (n - 1 )/(2*n) ))
	return rate	
	
if __name__ == "__main__" :

	# Plot rate of location update vs number of cells per location area for 
        # r = 0.1km, 1km , 10km  

	nla = np.arange(1.0, 10.0, 0.01)
	
	locationUpdate = np.vectorize(locationUp)
	plt.plot(nla, locationUpdate(0.1,nla),label="r=0.1km")
	plt.plot(nla, locationUpdate(1,nla),label="r=1km")
	plt.plot(nla, locationUpdate(10,nla),label="r=10km")
	plt.xlabel('Number of cells per location area per $km^{2}$, $N_{LA}$')
	plt.ylabel('Rate of location update per hour, $\lambda_{LU}$')
	plt.yscale('log',basex=10)
	plt.legend(loc=1, fontsize = 'small')
	plt.show()
