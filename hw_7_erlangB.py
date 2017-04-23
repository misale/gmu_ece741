#!/usr/bin/python
import numpy as np
import math
import matplotlib.pyplot as plt


def erlangB(A,N) :
	x = float(N)
	a = float(A)
	#print "number of chan : ",N
	if x == 1.0 :
		return a/(1.0 + a)
	else :
		invE = 1.0 + ((x/a)*(1/erlangB(a,x-1)))
		#print "N ", N, " P_B " ,1.0/invE 
		return 1.0/invE
		
	
if __name__ == "__main__" :
	plt.clf()

	# Plot P_b for N=1, 2, 5, 10, 20, 50,
	# 100 and infinity. (Note for N= infinity P_b was derived in class as A^N/N!(e^-A). 


	a = np.arange(1.0, 100.0, 0.1)
	#a = np.array([ x for x in range(100)])
	erlang = np.vectorize(erlangB)
	plt.plot(a, erlang(a,1),label="N=1")
	plt.plot(a, erlang(a,2),label="N=2")
	plt.plot(a, erlang(a,5),label="N=5")
	plt.plot(a, erlang(a,10),label="N=10")
	plt.plot(a, erlang(a,20),label="N=20")
	plt.plot(a, erlang(a,50),label="N=50")
	plt.plot(a, erlang(a,100),label="N=100")
	plt.plot(a, erlang(a,150),label="N=150")
	plt.xlabel('A')
	plt.ylabel('$P_{B}$')
	plt.text(10, erlangB(10,1), 'N=1')
	plt.text(10, erlangB(10,2), 'N=2')
	plt.text(10, erlangB(10,5), 'N=5')
	plt.text(10, erlangB(10,10), 'N=10')
	plt.text(20, erlangB(20,20), 'N=20')
	plt.text(50, erlangB(50,50), 'N=50')
	plt.text(100, erlangB(100,100), 'N=100')
	plt.text(150, erlangB(150,150), 'N=150')
	#plt.legend(loc=1, fontsize = 'x-small')
	plt.show()
