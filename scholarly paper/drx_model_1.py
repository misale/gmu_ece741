#!/usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt


def drxPSF(RHO,LAMBDA,N,t_DS,t_DL,t_I) :
	BETA = ( (1.0 - (math.e**(-1.0*LAMBDA*t_DS))**N) / (1.0 - (math.e**(-1.0*LAMBDA*t_DS))) )*t_DS
	GAMMA = ( ( (math.e**(-1.0*LAMBDA*t_DS))**N ) / (1.0 - math.e**(-1.0*LAMBDA*t_DL) ) )*t_DL
	PSF = ( 1.0 - RHO ) * ( (BETA) + (GAMMA) ) / (  BETA + GAMMA + ((1.0/LAMBDA)*((math.e**(LAMBDA*t_I)) - 1.0)) )
	return PSF

def PacketDelay(RHO,LAMBDA,N,t_DS,t_DL,t_I) :
	BETA = ( (1.0 - (math.e**(-1.0*LAMBDA*t_DS))**N) / (1.0 - (math.e**(-1.0*LAMBDA*t_DS))) )*t_DS
        GAMMA = ( ( (math.e**(-1.0*LAMBDA*t_DS))**N ) / (1.0 - math.e**(-1.0*LAMBDA*t_DL) ) )*t_DL
	D = (1.0 / (LAMBDA*(1.0 - RHO))) + (0.5*( (BETA*t_DS) + (GAMMA*t_DL) ) / ( BETA + GAMMA + ((1.0/LAMBDA)*((math.e**(LAMBDA*t_I)) - 1)) ))
	return D

def PacketDelaySimple(RHO,LAMBDA,N,t_DS,t_DL,t_I) :
        BETA = ( (1.0 - (math.e**(-1.0*LAMBDA*t_DS))**N) / (1.0 - (math.e**(-1.0*LAMBDA*t_DS))) )*t_DS
        GAMMA = ( ( (math.e**(-1.0*LAMBDA*t_DS))**N ) / (1.0 - math.e**(-1.0*LAMBDA*t_DL) ) )*t_DL
        D = (1.0 / (LAMBDA*(1.0 - RHO))) + (0.25*( (t_DS**2) + (t_DL**2) ) / ( BETA + GAMMA + ((1.0/LAMBDA)*((math.e**(LAMBDA*t_I)) - 1)) ))
        return D

def f(d):
    return  ((13*np.log10(d)) + (35*np.log10(1 + (d/90))))*(-1)

#simulation Variables - effect of t_I

LAMBDA = 0.1
N = 2
t_DS = 10
t_DL = 2 * t_DS
tau = 0.1
RHO = LAMBDA*tau

t_DS = np.arange(1.0, 256.0, 0.1)
#t_DS = np.array([ 2**x for x in range(0,9)])

#plt.plot(t_DS, drxPSF(RHO,LAMBDA,N,t_DS,2*t_DS,0),'-.',label="$t_{I}$=0msec")
#plt.plot(t_DS, drxPSF(RHO,LAMBDA,N,t_DS,2*t_DS,4),'g',label="$t_{I}$=4msec")
#plt.plot(t_DS, drxPSF(RHO,LAMBDA,N,t_DS,2*t_DS,8),'r',label="$t_{I}$=8msec")
#plt.plot(t_DS, drxPSF(RHO,LAMBDA,N,t_DS,2*t_DS,16),'y',label="$t_{I}$=16msec")
#plt.plot(t_DS, drxPSF(RHO,LAMBDA,N,t_DS,2*t_DS,32),'b',label="$t_{I}$=32msec")
#plt.plot(t_DS, drxPSF(RHO,LAMBDA,N,t_DS,2*t_DS,64),'black',label="$t_{I}$=64msec")
#plt.xscale('log',basex=2)
#plt.xlabel('$t_{DS}$(msec)')
#plt.ylabel('PS')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()

#plt.plot(t_DS, PacketDelaySimple(RHO,LAMBDA,N,t_DS,2*t_DS,0),'-.',label="$t_{I}$=0msec")
#plt.plot(t_DS, PacketDelaySimple(RHO,LAMBDA,N,t_DS,2*t_DS,4),'g',label="$t_{I}$=4msec")
#plt.plot(t_DS, PacketDelaySimple(RHO,LAMBDA,N,t_DS,2*t_DS,8),'r',label="$t_{I}$=8msec")
#plt.plot(t_DS, PacketDelaySimple(RHO,LAMBDA,N,t_DS,2*t_DS,16),'y',label="$t_{I}$=16msec")
#plt.plot(t_DS, PacketDelaySimple(RHO,LAMBDA,N,t_DS,2*t_DS,32),'b',label="$t_{I}$=32msec")
#plt.plot(t_DS, PacketDelaySimple(RHO,LAMBDA,N,t_DS,2*t_DS,64),'black',label="$t_{I}$=64msec")
#plt.xscale('log',basex=2)
#plt.xlabel('$t_{DS}$(msec)')
#plt.ylabel('Packet transmission delay(msec)')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()


#simulation Variables - effect of N

LAMBDA = 0.1
#N = 2
t_DS = 10
t_DL = 2 * t_DS
tau = 0.1
RHO = LAMBDA*tau
t_I = 10
N = np.array([ 2**x for x in range(0,5)])

#plt.plot(N, drxPSF(RHO,LAMBDA,N,t_DS,2*t_DS,t_I),'-.',label="$t_{DL}$=20msec")
#plt.plot(N, drxPSF(RHO,LAMBDA,N,t_DS,4*t_DS,t_I),'g',label="$t_{DL}$=40msec")
#plt.plot(N, drxPSF(RHO,LAMBDA,N,t_DS,8*t_DS,t_I),'r',label="$t_{DL}$=80msec")
#plt.plot(N, drxPSF(RHO,LAMBDA,N,t_DS,16*t_DS,t_I),'y',label="$t_{DL}$=160msec")
#plt.plot(N, drxPSF(RHO,LAMBDA,N,t_DS,32*t_DS,t_I),'b',label="$t_{DL}$=320msec")
#plt.plot(N, drxPSF(RHO,LAMBDA,N,t_DS,64*t_DS,t_I),'black',label="$t_{DL}$=640msec")
#plt.xscale('log',basex=2)
#plt.xlabel('N(number of DRX short cycles)')
#plt.ylabel('PS')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()

#plt.plot(N, PacketDelay(RHO,LAMBDA,N,t_DS,2*t_DS,t_I),'-.',label="$t_{DL}$=20msec")
#plt.plot(N, PacketDelay(RHO,LAMBDA,N,t_DS,4*t_DS,t_I),'g',label="$t_{DL}$=40msec")
#plt.plot(N, PacketDelay(RHO,LAMBDA,N,t_DS,8*t_DS,t_I),'r',label="$t_{DL}$=80msec")
#plt.plot(N, PacketDelay(RHO,LAMBDA,N,t_DS,16*t_DS,t_I),'y',label="$t_{DL}$=160msec")
#plt.plot(N, PacketDelay(RHO,LAMBDA,N,t_DS,32*t_DS,t_I),'b',label="$t_{DL}$=320msec")
#plt.plot(N, PacketDelay(RHO,LAMBDA,N,t_DS,64*t_DS,t_I),'black',label="$t_{DL}$=640msec")
#plt.xscale('log',basex=2)
#plt.xlabel('N(number of DRX short cycles)')
#plt.ylabel('Packet transmission delay(msec)')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()



#simulation Variables - effect of N on t_DS

LAMBDA = 0.1
N = 2
t_I = 40
t_DL = 128	
tau = 0.1
RHO = LAMBDA*tau

t_DS = np.arange(1.0, 128.0, 0.1)
#t_DS = np.array([ 2**x for x in range(2,20)])

plt.plot(t_DS, drxPSF(RHO,LAMBDA,2,t_DS,t_DL,t_I),'-.',label="N=2")
plt.plot(t_DS, drxPSF(RHO,LAMBDA,4,t_DS,t_DL,t_I),'g',label="N=4")
plt.plot(t_DS, drxPSF(RHO,LAMBDA,8,t_DS,t_DL,t_I),'r',label="N=8")
plt.plot(t_DS, drxPSF(RHO,LAMBDA,16,t_DS,t_DL,t_I),'y',label="N=16")
plt.plot(t_DS, drxPSF(RHO,LAMBDA,32,t_DS,t_DL,t_I),'b',label="N=32")
plt.plot(t_DS, drxPSF(RHO,LAMBDA,64,t_DS,t_DL,t_I),'black',label="N=64")
plt.xscale('log',basex=2)
plt.xlabel('$t_{DS}$(msec)')
plt.ylabel('PS')
plt.legend(loc=1, fontsize = 'x-small')
plt.show()

#plt.plot(t_DS, PacketDelay(RHO,LAMBDA,2,t_DS,t_DL,t_I),'-.',label="N=2")
#plt.plot(t_DS, PacketDelay(RHO,LAMBDA,4,t_DS,t_DL,t_I),'g',label="N=4")
#plt.plot(t_DS, PacketDelay(RHO,LAMBDA,8,t_DS,t_DL,t_I),'r',label="N=8")
#plt.plot(t_DS, PacketDelay(RHO,LAMBDA,16,t_DS,t_DL,t_I),'y',label="N=16")
#plt.plot(t_DS, PacketDelay(RHO,LAMBDA,32,t_DS,t_DL,t_I),'b',label="N=32")
#plt.plot(t_DS, PacketDelay(RHO,LAMBDA,64,t_DS,t_DL,t_I),'black',label="N=64")
#plt.xscale('log',basex=2)
#plt.xlabel('$t_{DS}$(msec)')
#plt.ylabel('Packet transmission delay(msec)')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()


#simulation Variables - effect of N on t_DL

LAMBDA = 0.1
N = 2
t_I = 10
t_DS = 10
tau = 0.1
RHO = LAMBDA*tau

t_DL = np.arange(1.0, 500.0, 0.1)
#t_DS = np.array([ 2**x for x in range(2,20)])

#plt.plot(t_DL, drxPSF(RHO,LAMBDA,2,t_DS,t_DL,t_I),'-.',label="N=2")
#plt.plot(t_DL, drxPSF(RHO,LAMBDA,4,t_DS,t_DL,t_I),'g',label="N=4")
#plt.plot(t_DL, drxPSF(RHO,LAMBDA,8,t_DS,t_DL,t_I),'r',label="N=8")
#plt.plot(t_DL, drxPSF(RHO,LAMBDA,16,t_DS,t_DL,t_I),'y',label="N=16")
#plt.plot(t_DL, drxPSF(RHO,LAMBDA,32,t_DS,t_DL,t_I),'b',label="N=32")
#plt.plot(t_DL, drxPSF(RHO,LAMBDA,64,t_DS,t_DL,t_I),'black',label="N=64")
#plt.yscale('log',basex=10)
#plt.xlabel('$t_{DL}$(msec)')
#plt.ylabel('PS')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()

#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,2,t_DS,t_DL,t_I),'-.',label="N=2")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,4,t_DS,t_DL,t_I),'g',label="N=4")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,8,t_DS,t_DL,t_I),'r',label="N=8")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,16,t_DS,t_DL,t_I),'y',label="N=16")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,32,t_DS,t_DL,t_I),'b',label="N=32")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,64,t_DS,t_DL,t_I),'black',label="N=64")
##plt.xscale('log',basex=2)
#plt.xlabel('$t_{DL}$(msec)')
#plt.ylabel('Packet transmission delay(msec)')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()

# PS vs t_I , for different N

LAMBDA = 0.1
N = 2
t_I = 10
t_DS = 10
tau = 0.1
RHO = LAMBDA*tau

t_I = np.arange(1.0, 20.0, 0.1)
#t_DS = np.array([ 2**x for x in range(2,20)])
t_DL = 64*t_DS

#plt.plot(t_I, drxPSF(RHO,LAMBDA,2,t_DS,t_DL,t_I),'-.',label="N=2")
#plt.plot(t_I, drxPSF(RHO,LAMBDA,4,t_DS,t_DL,t_I),'g',label="N=4")
#plt.plot(t_I, drxPSF(RHO,LAMBDA,8,t_DS,t_DL,t_I),'r',label="N=8")
#plt.plot(t_I, drxPSF(RHO,LAMBDA,16,t_DS,t_DL,t_I),'y',label="N=16")
#plt.plot(t_I, drxPSF(RHO,LAMBDA,32,t_DS,t_DL,t_I),'b',label="N=32")
#plt.plot(t_I, drxPSF(RHO,LAMBDA,64,t_DS,t_DL,t_I),'black',label="N=64")
#plt.xscale('log',basex=10)
#plt.xlabel('$t_{I}$(msec)')
#plt.ylabel('PS')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()

#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,2,t_DS,t_DL,t_I),'-.',label="N=2")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,4,t_DS,t_DL,t_I),'g',label="N=4")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,8,t_DS,t_DL,t_I),'r',label="N=8")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,16,t_DS,t_DL,t_I),'y',label="N=16")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,32,t_DS,t_DL,t_I),'b',label="N=32")
#plt.plot(t_DL, PacketDelay(RHO,LAMBDA,64,t_DS,t_DL,t_I),'black',label="N=64")
##plt.xscale('log',basex=2)
#plt.xlabel('$t_{DL}$(msec)')
#plt.ylabel('Packet transmission delay(msec)')
#plt.legend(loc=1, fontsize = 'x-small')
#plt.show()
