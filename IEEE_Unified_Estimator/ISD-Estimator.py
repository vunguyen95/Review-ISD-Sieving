import math
import numpy as np
from math import inf, floor, ceil, log2, exp
import matplotlib.pyplot as plt
from functools import lru_cache, cache
import sys

"""This estimation script computes the bit-complexity of various ISD algorithms, including: 
- Stern/Dumer
- MMT : 
- BJMM (depth 2): For most (sparse error) parameters considered, depth 2 gives the best results.

MO/B-M ISDs are not included due to unclear polynomial overhead. 
-------------------------------------------------------------
The complexity is estimated as:
	T = Number of iteration x bit_ops per iteration.
	M = (total list) x log(n)
Note that bit_ops is taking into account that the vector operations can be armotized with the low-weight.
This script does not include:
- the Syndrome Decoding Estimator due to different framework.
- Peters Gaussian elimination tricks.
"""


def truncate(x, precision):
	"""
	Truncates a float

	INPUT:

	- ``x`` -- value to be truncated
	- ``precision`` -- number of decimal places to after which the ``x`` is truncated

	"""

	return float(int(x * 10 ** precision) / 10 ** precision)


def binomial(n, k):
	return math.comb(int(n), int(k))
	
def binom(n, k):
	return math.comb(int(n), int(k))
def _mem_matrix(n, k, r):
	"""
	Memory usage of parity check matrix in vector space elements

	INPUT:

	- ``n`` -- length of the code
	- ``k`` -- dimension of the code
	- ``r`` -- block size of M4RI procedure

	EXAMPLES::

		>>> from .estimator import _mem_matrix
		>>> _mem_matrix(n=100,k=20,r=0) # doctest: +SKIP

	"""
	return n - k + 2 ** r

def _gaussian_elimination_complexity(n, k, r): 
	"""
	Complexity estimate of Gaussian elimination routine

	INPUT:

	- ``n`` -- Row additons are perfomed on ``n`` coordinates
	- ``k`` -- Matrix consists of ``n-k`` rows
	- ``r`` -- Blocksize of method of the four russian for inversion, default is zero

	[Bar07]_ Bard, G.V.: Algorithms for solving linear and polynomial systems of equations over finite fields
	with applications to cryptanalysis. Ph.D. thesis (2007)

	[BLP08] Bernstein, D.J., Lange, T., Peters, C.: Attacking and defending the mceliece cryptosystem.
	In: International Workshop on Post-Quantum Cryptography. pp. 31–46. Springer (2008)

	EXAMPLES::

		>>> from .estimator import _gaussian_elimination_complexity
		>>> _gaussian_elimination_complexity(n=100,k=20,r=1) # doctest: +SKIP

	"""

	if r != 0:
		return (r ** 2 + 2 ** r + (n - k - r)) * int(((n + r - 1) / r))

	return (n - k) ** 2

#copied from Sieving ISD repo
def _optimize_m4ri(n, k, mem=1e10):
	"""
	Find optimal blocksize for Gaussian elimination via M4RI

	INPUT:

	- ``n`` -- Row additons are perfomed on ``n`` coordinates
	- ``k`` -- Matrix consists of ``n-k`` rows

	"""

	(r, v) = (0, 1e10)
	for i in range(n - k):
		tmp = log2(_gaussian_elimination_complexity(n, k, i))
		if v > tmp and r < mem:
			r = i
			v = tmp
	return r

#copied from Sieving ISD repo
def gauss(n,k):
	r = _optimize_m4ri(n,k, 30)
	return  _gaussian_elimination_complexity(n, k, r)*n
	
#From the reviewer  
def opt_dumer(n,k,w):
	"""
	this function calculates the cost of Stern-Dumer in a similar way as the cost of Sieving ISD is computed in the paper
	
	not all operations are performed on vectors of size n, indeed most of the operations happen only on the small instance use sparse test patterns.

	The claim is NOT that the cost should be computed in this way.
	This is only for comparison.
	"""

	time = 1e6
	r = _optimize_m4ri(n, k, mem = 30)
	# naive cost of M4RI, might be slightly faster to reuse pivots as suggested by Bernstein-Lange-Peters
	#Tg = gauss(n,k)
	Tg = _gaussian_elimination_complexity(n, k, r)*n

	for l in range(90):
		if (k+l)%2: continue # let's use only even-sized small instances

		k1 = (k+l)/2
		for p in range(1,6): # weight on half of small instance

			# time required for permutations
			Tp = max(log2(binomial(n, w)) - log2(binomial(n - k - l, w - 2 * p)) - 2*log2(binomial(k1, p)), 0)

			# list size
			L = binomial(k1,p)

			# List1: 
			# for each pattern calculate syndrome of small instance -> cost p*l per pattern (actually p*l is the naive cost of syndrome computation, can be done slightly faster with revolving door algorithm)
			# then put index of pattern into hashmap -> log2(L) bits stored per pattern
	
			C1 = L*(p*l + log2(L))

			# List2: calculate syndrome of small instance and check hashmap
			# for each collision in the hashamp extend to complete instance -> calculate remaining syndrome elements only for collisions
			C2 = C1 +  L**2 / 2**l * (n-k-l)*2*p

			# overall cost in RAM 
			tmp = Tp + log2(Tg + C1 + C2)

			if tmp < time: 
				time = tmp
				Memory = log2(2*L + _mem_matrix(n,k,r))
				params = {"p": p, "l" : l, "Memory:" : Memory}
				# required memory: store list 1 in hashmap -> store index of test pattern, requires log2(L) bits
				# for ease of implmentation overdimension hashmap by factor 2
				mem = log2(L*log2(L) * 2)

			time = min(time, tmp)

	return time, params, mem 

	
#Stern complexity modified according to the reviewer
def opt_stern(n,k,w):


	"""this function calculates the cost of Stern-Dumer in a similar way as the cost of Sieving ISD is computed in the paper
	
	not all operations are performed on vectors of size n, indeed most of the operations happen only on the small instance use sparse test patterns.

	The claim is NOT that the cost should be computed in this way.
	This is only for comparison.
	"""

	time = 1e6
	r = _optimize_m4ri(n, k, mem = 30)
	# naive cost of M4RI, might be slightly faster to reuse pivots as suggested by Bernstein-Lange-Peters
	#Tg = gauss(n,k)
	Tg = _gaussian_elimination_complexity(n, k, r)*n

	for l in range(90):

		k1 = k/2
		for p in range(1,10): # weight on half of small instance

			# time required for permutations
			Tp = max(log2(binomial(n, w)) - log2(binomial(n - k - l, w - 2 * p)) - 2*log2(binomial(k1, p)), 0)

			# list size
			L = binomial(k1,p)

			# List1: 
			# for each pattern calculate syndrome of small instance -> cost p*l per pattern (actually p*l is the naive cost of syndrome computation, can be done slightly faster with revolving door algorithm)
			# then put index of pattern into hashmap -> log2(L) bits stored per pattern
	
			C1 = L*(p*l + log2(L))

			# List2: calculate syndrome of small instance and check hashmap
			# for each collision in the hashamp extend to complete instance -> calculate remaining syndrome elements only for collisions
			C2 = C1 +  L**2 / 2**l * (n-k-l)*2*p

			# overall cost in RAM 
			tmp = Tp + log2(Tg + C1 + C2)

			if tmp < time: 
				time = tmp
				Memory = log2(2*L + _mem_matrix(n,k,r))
				params = {"p": p, "l" : l, "Memory" : Memory, "Memory (bits)" : Memory + log2(n)}
				# required memory: store list 1 in hashmap -> store index of test pattern, requires log2(L) bits
				# for ease of implmentation overdimension hashmap by factor 2
				mem = log2(L*log2(L) * 2)

			time = min(time, tmp)

	return time, params, mem
	
#mmt algorithm

	
def mmt_baldi(n,k,w, mem = inf):
	"""
	MMT concrete complexity analysis according to Baldi et al + Esser et al. fixes.
	
	This is only for comparison.
	"""

	time = 1e6
	r = _optimize_m4ri(n, k, 30)
	# naive cost of M4RI, might be slightly faster to reuse pivots as suggested by Bernstein-Lange-Peters
	#Tg = gauss(n,k)
	Tg = _gaussian_elimination_complexity(n, k, r)*n

	for p in range(4,w,4):
		for l in range(0,n-k-(w-p)):
			k1 = int((k+l)/2)
			#print(p)
			l1 = int(log2(binom(k1,p/4)))
			l2 = l - l1

			# time required for permutations
			Tp = max(log2(binomial(n, w)) - log2(binomial(n - k - l, w -  p)) - log2(binomial((k+l), p)), 0)

			# list size
			L = binomial(k1,p/4)
			L12 = max(1, L ** 2 // 2 ** l1)
			tmp_mem = log2((2 * L + L12) + _mem_matrix(n, k, r))
			if tmp_mem > mem:
				continue 
			
			T_L = L* p/4 * l1
			#T_tree = T_L + min(L, binom(k1, p/2)/ binom(p, p/2)) * (p/4 * l1 + L/2**l1 * p/2 * l2) + 2* T_L + L**2/2**l1*(p/2* l2 + binom(k+l, p/2)*binom(p, p/2)/ 2**l2 * (n-k-l) * p) 
			T_tree = min(L, (binom(k1,p/2)/binom(p,p/2)))*(p/4 * l1 + L/(2**l1) * p/2 * l2) + L*p/2 * l1 + L*(p/4 * l1 + L/(2**l1) * p/2 * l2 + L* binom(k+l, p/2)/ 2**(l1+l2) / binom(p,p/2) * p * (n-k-l))
			
			#B-E fix:            
			random = L / binom(p, p/2)
			# overall cost in RAM 
			tmp = Tp + log2(Tg + random*T_tree)

			if tmp < time: 
				time = tmp
				params = {"p": p, "l" : l, "memory (list size)" : tmp_mem}
				# required memory: Don't care right now, remember to adjust later

			time = min(time, tmp)

	return time, params

def mmt_new(n,k,w,mem = inf):
	"""
	MMT concrete complexity analysis, using for comparision with GJN algorithm, very close to mmt_baldi
	REF: MMT, Decoding Random Linear code in 2^{0.054n}
	"""
	time = inf
	memory = 0
	params = {}
	r = _optimize_m4ri(n,k,mem)
	Tg = _gaussian_elimination_complexity(n,k,r)*n
	
	for p in range(4, 40, 4):
		rep = binom(p, p/2) #The number of possible way to split the error vectors.
		l2 = int(floor(log2(rep))) # the number of bits to match in the middle layer
		
		for l1 in range(2, min(n-k - w + p - l2, 500)):
			l = l1 + l2
			
			#Top list, index set of size (k+l)/2, weight precisely p/4
			L2 = binom((k+l)/2, p/4)
			#Middle list
			L1 = L2**2 /2**l2
			#Bottom
			L = L1**2 / 2**l1
			
			#Memory, in terms of list size
			tmp_mem = (log2(2*L2 + L1 + _mem_matrix(n, k, r)))
			if tmp_mem > mem:
				continue
			#complexity of the tree
			T2 = 4 * (L2) * (p/4 * l2 + log2(L2)) 	#compute l2 syndrome in 4 lists, each vector weight p/4
			T1 = 2 * L1 * (p/2 * l1 + log2(L1)) 	# for each matching (weight p/2) in the middle layer, compute l1 syndrome.
			T0 = L * p * (n-k-l) 		            # For each solution candidate, check the weight condition.
			T_tree = T2 + T1 + T0

			#Split Probability : The way of constructing the top list lower the number of representation for each solution
			#in MMT case, each vector in the middle layer list has weight precisely p/4 in each half of the interval 0, k+l 
			#(as opposed to uniformly sampling p/2 over the whole interval.)
			P_split = binom(p/2, p/4)**2 / binom(p,p/2)	
			#Number of iteration:
			# 1) Probability of good weight distribution
			# 2) Splitting probability due to constructing base list
			# 3) probability of at least one representation survive: For MMT algorithm, this probability is about 1/2 (heuristic, cit. Theorem 2, page 14,
			# Decoding Random Linear code in 2^0.054n)

			Tp =  binom((k+l)/2, p/2)**2 * binom(n-k-l, w - p) / binom(n,w) #Good Permutation
			Tp = Tp * P_split * 1/2

			
			tmp_complexity = log2((Tg + T_tree)/Tp)
			if tmp_complexity < time:
					time = tmp_complexity
					params = {"p": p, "l" : l, "l2": l2, "l1": l1, "memory (list size)" : tmp_mem, "Memory (bits)" : tmp_mem + log2(n)}
	return time, params
	
def bjmm_depth_2_new(n,k,w,mem = inf):
	"""
	BJMM in depth 2 concrete complexity analysis, using for comparision with GJN algorithm.
	Optimal BJMM depth for interested parameters is often 2.
	[BJMM12] Becker, A., Joux, A., May, A., Meurer, A.: Decoding random binary linear codes in 2^(n/20): How 1+ 1= 0
		improves information set decoding. In: Annual international conference on the theory and applications of
		cryptographic techniques. pp. 520–536. Springer (2012)
	Alexander Meurer Dissertation: A coding-theoretic Approach to Cryptanalysis.

	"""
	#solutions = max(0, log2(binom(n, w)) - (n - k))
	time = inf
	memory = 0
	params = {}
	r = _optimize_m4ri(n, k, mem)
	Tg = _gaussian_elimination_complexity(n, k, r)*n
	
	for eps in range(0, 12, 2): #epsilon in layer 1, weight p/2 + eps
		for p in range (2, 30, 2):
			p1 = int(p//2) + eps
			for l in range (50, min(n-k- w + p, 500)):
				rep = binom(p, p/2) * binom(k+l - p, eps)
				l2 = int(floor(log2(rep)))
				if l2 > l:
					continue
				l1 = l - l2
				
				#Base list, pick random index set of size (k+l)/2, weight p1/2
				L2  = binom((k+l)/2, p1/2)
				#Middle list
				L1 = binom(k+l, p1) / rep
				
				#expected matchings
				C2 = L2**2 // 2**l2 #C2 == L1, if l2 is set correctly.
				C1 = L1**2 // 2**l1
				#print(log2(C2), p, p1, l, l1, l2)
				
				#Memory, in terms of list size
				tmp_mem = log2(2*L2 + L1 + _mem_matrix(n, k, r))
				if tmp_mem > mem:
					continue
				
				#complexity of the tree
				T2 = 4 * L2 * ((p1/2) * l2 + log2(L2)) 	#compute l2 syndrome bits
				T1 = 2*C2 * ((p1) * l1 + log2(L1)) 		#compute l1 syndrome bits for each matching		
				T0 = C1 * p * (n-k-l)		            #checking solution
				T_tree = T2 + T1 + T0

				#Splitting probability of middle layer
				# In BJMM case, each of the 4 top list is sampled by partitions Pi, Pj of size k+l/2, each contains p1/2 ones.
				# For each random partition, a vector in the middle layer can be represented with prob P_split
				# ``we can guarantee independent splitting conditions for all the e(1)_i yielding a total splitting probability of P_slit^2 ''.
				P_split = binom((k+l)/2, p1/2)**2 / binom(k+l, p1) 
				
				#Number of iterations:
				#1) Good permutation: log2(binom(k+l,p)) + log2(binom(n-k-l, w-p)) - log2(binom(n,w))
				#2) Splitting Probability, as above
				#3) PRobability that at least one representation survive the algorithm
				# Cit. Alexander Meuer Dissertation, p.134, formula 8.15: (1-exp(-2^{-l2}*rep)), where l2 = floor(log2(rep)).
				#Tp = P_split**2 * binom(k+l, p) * binom(n-k-l, w-p) / binom(n,w)
				Tp = 2*log2(P_split) + log2(binom(k+l,p)) + log2(binom(n-k-l, w-p)) - log2(binom(n,w)) + log2(1-exp(-2**(-l2)*rep))
				tmp_complexity = log2((Tg + T_tree)) - Tp
				
				if tmp_complexity < time:
					time = tmp_complexity
					params = {"p": p, "p1": p1, "eps" : eps, "l" : l, "l2": l2, "l1": l1, "P_split": 2*log2(P_split), "reps sucess": log2(1- exp(-2**(-l2)*rep )), "memory (list size)": tmp_mem, "Memory (bits)" : tmp_mem + log2(n)}
	return time, params
					
print('#----# TEST #----#')

# ignoring DOOM, i.e., quasi-cyclic structure
McEliece1 = {"n": 3488, "k" : 2720, "w" : 64}
McEliece3 = {"n": 4608, "k" : 3360, "w" : 96}
McEliece5_0 = {"n": 6688, "k" : 5024, "w" : 128}
McEliece5_1 = {"n": 6960, "k" : 5413, "w" : 119}
McEliece5_2 = {"n": 8192, "k" : 6527, "w" : 128}
MDPC1_0 = {"n": 24646, "k" : 12323, "w" : 134}
MDPC2_0 = {"n": 49318, "k" : 24659, "w" : 199}
MDPC3_0 = {"n": 81946, "k" : 40973, "w" : 264}

MDPC1_1 = {"n": 24646, "k" : 12323, "w" : 142}
MDPC2_1 = {"n": 49318, "k" : 24659, "w" : 206}
MDPC3_1 = {"n": 81946, "k" : 40973, "w" : 274}

MDPC1_2 = {"n": 35338, "k" : 17669, "w" : 132}
MDPC2_2 = {"n": 71702, "k" : 35851, "w" : 200}
MDPC3_2 = {"n": 115274, "k" : 57637, "w" : 262}
params = MDPC1_0
n = params["n"]
k = params["k"]
w = params["w"]

mem = 55 - log2(n)
#mem  = inf
#print("MEMORY LIMIT: ", mem +log2(n) )

print(n)
print(k)
print(w)

dumer_GJN = opt_dumer(n,k,w)
T_dumer_GJN, params_GJN = dumer_GJN[0], dumer_GJN[1]

print("Complexities Analysis of ISD algorithms")
print(f'T_dumer_reviewer = {truncate(T_dumer_GJN,2)}')
print("params: ", params_GJN)

stern_GJN = opt_stern(n,k,w)
T_Stern_GJN, params_stern_GJN = stern_GJN[0], stern_GJN[1]

print(f'T_Stern_reviewer = {truncate(T_Stern_GJN,2)}')
print("params: ", params_stern_GJN)

print("----------------------")
#mmt_test = mmt_baldi(n,k,w, mem)
#T_mmt_test = mmt_test[0]
#params_mmt_test = mmt_test[1]
#print("mmt_Baldi =: ", truncate(T_mmt_test,2))
#print(params_mmt_test)


print("----------------------------")
bjmm = bjmm_depth_2_new(n,k,w, mem)
T_bjmm = bjmm[0]
params_bjmm = bjmm[1]
print(f'BJMM = {truncate(T_bjmm,2)}')
print(params_bjmm)

print("----------------------------")
mmt = mmt_new(n,k,w, mem)
T_mmt = mmt[0]
params_mmt = mmt[1]
print(f'MMT = {truncate(T_mmt,2)}')
print(params_mmt)
