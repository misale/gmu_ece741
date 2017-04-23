#!/usr/bin/python
import numpy as np
import math
import matplotlib.pyplot as plt

# for infinite iteration defining infinity
inf = 1000
lambda_ipc = 1.0/30.0
lambda_is = 1.0/2000.0
# inter packet arrival time
lambda_x = 10.0 
#number of packets calls within a packet service session
mu_p = 25.0
mu_pc = 5.0
k = 1.0 
tDS =  2.0
tDL = 10.0 
tau = 0.1

#The packet calls may be the interpacket call idle time (tipc) with probability Ppc
p_pc = 1.0 - (1.0/mu_pc)
#or the intersession idle time (tis) with probability Ps
p_s = 1.0 /mu_pc

#Adjustable DRX Timer State for Light Sleep
def n_ds(k) :
	n = []
	N_DS = 0
	if k == 1 :
        	n = range(1,10)
        	N_DS = 9
	if k == 5 :
        	n = range(1,7)
        	N_DS = 6

	return n,N_DS

# The duration of the kth sleep cycle, which consists of a sleep interval and a listening
def c_ds(x,k) :
	n,N_DS = n_ds(k)
	max = n[-1]
	if x < max :
		return k * (2**x)
	else :
		return k * (2**max)
#the transition from Light Sleep to Deep Sleep consist of t_CS
def tn(n,k) :
	t_CS = 0.0
	for i in n :
       		t_CS += c_ds(i,k)

	t_CS = k*t_CS

	tN = t_CS
	return tN

#State 1 to State 1 and State 1 to State 2
#The probability that a new packet call begins before the expiration of tI is q_1 - condition1

def q_1(tI) :
	return 1.0 - (math.e**((-1)*lambda_ipc*tI))

# condition2 - q_2

def q_2(tI) : 
	return 1.0 - (math.e**((-1)*lambda_is*tI))

# ipc prod value

def ipc_prod(k) :
	n,N_DS = n_ds(k)
	result = 1.0
	for i in n :
		result *= math.e**(-1*k*lambda_ipc*c_ds(i,k))
	return result

def is_prod(k) :
	n,N_DS = n_ds(k)
	result = 1.0
	for i in n :
		result*= math.e**(-1*k*lambda_is*c_ds(i,k))
	return result


# holding times

def mean_holding_state_1(tI):
    result = (mu_p/lambda_x) + (lambda_ipc**(-1))*(1.0 - math.e**((-1)*lambda_ipc*tI)) + (lambda_is**-1)*(1 - math.e**((-1)*lambda_is*tI))
    return result

def prob_s(j,tI,tDS) :
	return p_pc*(math.e**(-1*lambda_ipc*tI))*(math.e**(-1*lambda_ipc*(j - 1.0)*tDS))*(1.0 - math.e**(-1*lambda_ipc*tDS)) + p_s*(math.e**(-1*lambda_is*tI))*(math.e**(-1*lambda_is*(j - 1.0)*tDS))*(1.0 - math.e**(-1*lambda_is*tDS))	

def prob_l(j,k,tDS,tDL) :
	n,N_DS = n_ds(k)
	return p_pc*(math.e**(-1*lambda_ipc*(tI + N_DS*tDS + (j - N_DS - 1.0)*tDL)))*(1.0 - math.e**(-1*lambda_ipc*tDL)) + p_s*(math.e**(-1*lambda_is*(tI + N_DS*tDS + (j - N_DS - 1.0)*tDL)))*(1.0 - math.e**(-1*lambda_is*tDL))

def prob_adj_s(j,tI,k) :
	return p_pc*(math.e**(-1*lambda_ipc*tI))*ipc_prod(j)*(1.0 - math.e**(-1*k*lambda_ipc*c_ds(j,k))) + p_s*(math.e**(-1*lambda_is*tI))*is_prod(j)*(1.0 - math.e**(-1*lambda_is*c_ds(j,k)))

def prob_adj_l(j,k,tI,tDL,tN) :
	n,N_DS = n_ds(k)
	return p_pc*(math.e**(-1*lambda_ipc*(tI + tN + (j - N_DS - 1)*tDL)))*(1.0 - math.e**(-1*lambda_ipc*tDL)) + p_s*(math.e**(-1*lambda_is*(tI + tN + (j - N_DS - 1.0)*tDL)))*(1.0 - math.e**(-1*lambda_is*tDL))

def sum_cds_over_n(n,N_DS,k) :
	result = 0.0
	for i in n:
		result += (c_ds(i,k) / i )
	return result

def mean_holding_state_2(k,n,N_DS): 
	return ( ( p_pc*ipc_prod(k) + p_s*is_prod(k) )*N_DS + ( p_pc*(1 - ipc_prod(k)) + p_s*(1.0 - is_prod(k)) )*( (p_pc/(1.0 - ipc_prod(k))) +  (p_s/(1.0 - is_prod(k))) ) )*k*sum_cds_over_n(n,N_DS,k)

def mean_holding_state_3(k,tDL): 
	return ( (p_pc/(1.0 - math.e**(-1*lambda_ipc*tDL)))   + ( p_s/(1.0 - math.e**(-1*lambda_is*tDL)) )  )*tDL

def mean_holding_state_3_effective(k,tDL): 
	return ( (p_pc/(1.0 - math.e**(-1*lambda_ipc*tDL)))   + ( p_s/(1.0 - math.e**(-1*lambda_is*tDL)) )  )*(tDL - tau)

def mean_holding_state_2_effective(k,n,N_DS):
	return ( ( p_pc*ipc_prod(k) + p_s*is_prod(k) )*N_DS + ( p_pc*(1 - ipc_prod(k)) + p_s*(1.0 - is_prod(k)) )*( (p_pc/(1.0 - ipc_prod(k))) +  (p_s/(1.0 - is_prod(k))) ) )*(k*sum_cds_over_n(n,N_DS,k) - tau)


# psf 

def psf(tI,k,n,N_DS,tDL):
	# markov chain p_11

	p_11 = p_pc*q_1(tI) + p_s*q_2(tI)

	# markov chain p_12

	p_12 = p_pc*(1 - q_1(tI)) + p_s*(1-q_2(tI))

	# markov chain p_21
	#print ipc_prod(k)
	p_21 = p_pc*(1.0 - ipc_prod(k)) + p_s*(1.0 - is_prod(k))

	# markov chain p_23

	p_23 = p_pc*ipc_prod(k) + p_s*is_prod(k)

	# markov chain p_31

	p_31 = 1


	# balance equations

	pi_1 = 1.0 / ( 1.0 + p_12 + (p_12*p_23) )

	pi_2 = p_12 / ( 1.0 + p_12 + (p_12*p_23) )

	pi_3 = (p_12*p_23) / (1 + p_12 + (p_12*p_23))
 
	return ((pi_2*mean_holding_state_2_effective(k,n,N_DS)) + (pi_3*mean_holding_state_3_effective(k,tDL)) ) / (pi_1*mean_holding_state_1(tI) + pi_2*mean_holding_state_2(k,n,N_DS) + pi_3*mean_holding_state_3(k,tDL))


def mean_delay_adj(k,n,N_DS,tI,tDL,tN) :
	s_sum = 0.0
	l_sum = 0.0
	for i in range(1,N_DS+1) :
		s_sum += (prob_adj_s(i,tI,k)*c_ds(i,k) / (2.0*i))
	for i in range(1,inf) :
                l_sum += (prob_adj_l(i,k,tI,tDL,tN)*tDL / (2.0*i))

	return s_sum + l_sum

def mean_delay(k,n,N_DS,tI,tDS,tDL) :
        s_sum = 0.0
        l_sum = 0.0
        for i in range(1,N_DS+1) :
                s_sum += (prob_s(i,tI,tDS)*c_ds(i,k) / (2.0*i))
        for i in range(1,inf) :
                l_sum += (prob_l(i,k,tDS,tDL)*tDL / (2.0*i))

        return s_sum + l_sum

#tN = tn(n,k)

### effect of tI#
#tI = np.arange(1.0, 400.0, 0.5)
#plt.plot(tI, psf(tI,1,n_ds(1)[0],n_ds(1)[1],tDL),'b',linestyle="dashed",label="$c_{DS}$=$2^N$")
#plt.plot(tI, psf(tI,5,n_ds(5)[0],n_ds(5)[1],tDL),'r',linestyle="dashed",label="$c_{DS}$=k*$2^N$")
#plt.plot(tI, psf(tI,1,n_ds(1)[0],10,tDL),'g',linestyle="solid",label="$T_N$=10")
#plt.plot(tI, psf(tI,1,n_ds(1)[0],20,tDL),'r',linestyle="solid",label="$T_N$=20")
#plt.plot(tI, psf(tI,1,n_ds(1)[0],40,tDL),'b',linestyle="solid",label="$T_N$=40")

#plt.plot(tI, mean_delay_adj(1,n_ds(1)[0],n_ds(1)[1],tI,tDL,tn(n_ds(1)[0],1)),'b',linestyle="dashed",label="$c_DS$=$2^N$")
#plt.plot(tI, mean_delay_adj(5,n_ds(5)[0],n_ds(5)[1],tI,tDL,tn(n_ds(5)[0],5)),'r',linestyle="dashed",label="$c_(DS)$=k*$2^N$")
#plt.plot(tI, mean_delay(1,n_ds(1)[0],20,tI,tDS,tDL),'g',linestyle="solid",label="$T_N$=20")
#plt.plot(tI, mean_delay(1,n_ds(1)[0],30,tI,tDS,tDL),'r',linestyle="solid",label="$T_N$=30")
#plt.plot(tI, mean_delay(1,n_ds(1)[0],40,tI,tDS,tDL),'b',linestyle="solid",label="$T_N$=40")
#plt.xlabel('tI(sec)')
#plt.ylabel('psf')
#plt.legend()
#plt.show()

### effect of N_DS
#tI = 2.0
#N_DS = tI = np.arange(1.0, 100.0, 0.5)
#plt.plot(tI, psf(tI,1,n_ds(1)[0],N_DS,tDL),'b',linestyle="dashed",label="$c_{DS}$=$2^N$")
#plt.plot(tI, psf(tI,5,n_ds(5)[0],N_DS,tDL),'r',linestyle="dashed",label="$c_{DS}$=k*$2^N$")
#plt.plot(tI, psf(tI,1,n_ds(1)[0],10,tDL),'g',linestyle="solid",label="$T_N$=10")
#plt.plot(tI, psf(tI,1,n_ds(1)[0],20,tDL),'r',linestyle="solid",label="$T_N$=20")
#plt.plot(tI, psf(tI,1,n_ds(1)[0],40,tDL),'b',linestyle="solid",label="$T_N$=40")
#plt.xlabel('T_N')
#plt.ylabel('psf')
#plt.legend()
#plt.show()

### effect of tDL
tI = 2.0
tDS = 2.0
tau = 0.1
tDL = np.arange(1, 10.0, 0.1)
#
plt.plot(tDL, psf(tI,1,n_ds(1)[0],n_ds(1)[1],tDL),'b',linestyle="dashed",label="$c_{DS}$=$2^N$")
plt.plot(tDL, psf(tI,5,n_ds(5)[0],n_ds(5)[1],tDL),'r',linestyle="dashed",label="$c_{DS}$=k*$2^N$")
plt.plot(tDL, psf(tI,1,n_ds(1)[0],10,tDL),'g',linestyle="solid",label="$T_N$=10")
plt.plot(tDL, psf(tI,1,n_ds(1)[0],20,tDL),'r',linestyle="solid",label="$T_N$=20")
plt.plot(tDL, psf(tI,1,n_ds(1)[0],40,tDL),'b',linestyle="solid",label="$T_N$=40")
plt.xlabel('$t_{DL}$')
plt.ylabel('psf')
plt.legend()
plt.show()


#plt.plot(tDL, mean_delay_adj(1,n_ds(1)[0],n_ds(1)[1],tI,tDL,tn(n_ds(1)[0],1)),'b',linestyle="dashed",label="$c_DS$=$2^N$")
#plt.plot(tDL, mean_delay_adj(5,n_ds(5)[0],n_ds(5)[1],tI,tDL,tn(n_ds(5)[0],5)),'r',linestyle="dashed",label="$c_(DS)$=k*$2^N$")
#plt.plot(tDL, mean_delay(1,n_ds(1)[0],20,tI,tDS,tDL),'g',linestyle="solid",label="$T_N$=20")
#plt.plot(tDL, mean_delay(1,n_ds(1)[0],30,tI,tDS,tDL),'r',linestyle="solid",label="$T_N$=30")
#plt.plot(tDL, mean_delay(1,n_ds(1)[0],40,tI,tDS,tDL),'b',linestyle="solid",label="$T_N$=40")
#plt.xlabel('$t_{DL}$')
#plt.ylabel('Delay')
#plt.legend()
#plt.show()
