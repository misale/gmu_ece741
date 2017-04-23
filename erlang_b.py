#!/usr/bin/python
import math

def erlangB(A,N) :
	x = float(N)
	a = float(A)
	#print "number of chan : ",N
	if x == 1.0 :
		return a/(1.0 + a)
	else :
		invE = 1.0 + ((x/a)*(1/erlangB(a,x-1)))
		print "N ", N, " P_B " ,1.0/invE 
		return 1.0/invE
		
	
if __name__ == "__main__" :

	print erlangB(37.9,50)
