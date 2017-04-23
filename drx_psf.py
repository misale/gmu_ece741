#!/usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt

#fixed Variables
RHO = 0.01 
LAMBDA = 0.1
N = 2
C_S = 8
C_L = 2 * C_S
C_T = 16
TAU = 0.1



def drxPsfCs(C_S,C_T) :
	psf = ( (1 - RHO)* ( ( (1 - math.e**(-1*LAMBDA*C_S*N)) / (1 - math.e**(-1*LAMBDA*C_S))  )*C_S + ( (math.e**(-1*LAMBDA*C_S*N)) / (1 - math.e**(-1*LAMBDA*C_L)) )*C_L ) ) / ( (( (1 - math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_S)) )*C_S)     + (( (math.e**(-1*LAMBDA*C_S*N)) / (1 - math.e**(-1*LAMBDA*C_L))  )*C_L) + ((1/LAMBDA)*( (math.e**(LAMBDA*C_T)) - 1)) )
	return psf

def PacketDelayMs(C_S,C_T) :
	delay = ((1/(LAMBDA*(1-RHO)))  + ( (0.5*( ( (1 - math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_S)) )*(C_S**2) + ( (math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_L)) )*(C_L**2)  )) / ( (( (1 - math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_S)) )*C_S) + (( (math.e**(-1*LAMBDA*C_S*N)) / (1 - math.e**(-1*LAMBDA*C_L))  )*C_L) + ((1/LAMBDA)*( (math.e**(LAMBDA*C_T)) - 1)) ) )
	#delay = ((LAMBDA*(TAU + TAU**2))/(2*(1-RHO)))  + ( (0.5*( ( (1 - math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_S)) )*(C_S**2) + ( (math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_L)) )*(C_L**2)  )) / ( (( (1 - math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_S)) )*C_S) + (( (math.e**(-1*LAMBDA*C_S*N)) / (1 - math.e**(-1*LAMBDA*C_L))  )*C_L) + ((1/LAMBDA)*( (math.e**(LAMBDA*C_T)) - 1)) ) )
	#delay = (((LAMBDA*(TAU + TAU**2))/(2*(1-RHO)))  +  0.5*( ( (1 - math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_S)) )*(C_S**2) + ( (math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_L)) )*(C_L**2)  )) / ( (( (1 - math.e**(-1*LAMBDA*C_S*N))/(1 - math.e**(-1*LAMBDA*C_S)) )*C_S) + (( (math.e**(-1*LAMBDA*C_S*N)) / (1 - math.e**(-1*LAMBDA*C_L))  )*C_L) + ((1/LAMBDA)*( (math.e**(LAMBDA*C_T)) - 1)) ) 
	return delay
def f(d):
    return  ((13*np.log10(d)) + (35*np.log10(1 + (d/90))))*(-1)

d = np.arange(1.0, 500.0, 1)
#d = np.array([ 2**x for x in range(2,9)])

#plt.plot(d, PacketDelayMs(d,0),'b.',label="0")
#plt.plot(d, PacketDelayMs(d,4),'r',label="4")
#plt.plot(d, PacketDelayMs(d,8),'y',label="8")
#plt.plot(d, PacketDelayMs(d,16),'c',label="16")
#plt.plot(d, PacketDelayMs(d,32),'g',label="32")
#plt.plot(d, PacketDelayMs(d,64),'-',label="64")

plt.plot(d, drxPsfCs(d,0),'b',linestyle="dashed",label="0")
plt.text(10,PacketDelayMs(500,0),"C_T=0",fontsize=12,color="b")
plt.plot(d, drxPsfCs(d,4),'r',linestyle="dashed",label="4")
plt.plot(d, drxPsfCs(d,8),'y',label="8")
plt.plot(d, drxPsfCs(d,16),'c',label="16")
plt.plot(d, drxPsfCs(d,32),'g',label="32")
plt.plot(d, drxPsfCs(d,64),'-',label="64")

plt.xscale('log',basex=2)
plt.xlabel('Short Cycle ms')
plt.ylabel('psf')
#plt.legend()
plt.show()
