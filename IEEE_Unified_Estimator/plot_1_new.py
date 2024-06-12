import math
import numpy as np
from math import inf, floor, ceil, log2, exp
import matplotlib.pyplot as plt
from functools import lru_cache
import sys

# parameters
c_label = 2
McEliece1 = {"n": 3488, "k" : 2720, "w" : 64}
McEliece3 = {"n": 4608, "k" : 3360, "w" : 96}
McEliece5_0 = {"n": 6688, "k" : 5024, "w" : 128}
params = McEliece1
n = params["n"]
k = params["k"]
w = params["w"]
M_overhead = 1.5



@lru_cache(maxsize = 12800000)
def comb(N,k):
    val= math.factorial(N)/(math.factorial(k)*math.factorial(N-k))
    return val
def binomial(n, k):
    return math.comb(int(n), int(k))
def binom(n, k):
    return math.comb(int(n), int(k))
def log2(N):
    return math.log(N,2)

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
@lru_cache(maxsize = 12800000)
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

@lru_cache(maxsize = 12800000)
def _optimize_m4ri(n, k, mem=inf):
    """
    Find optimal blocksize for Gaussian elimination via M4RI

    INPUT:

    - ``n`` -- Row additons are perfomed on ``n`` coordinates
    - ``k`` -- Matrix consists of ``n-k`` rows

    """

    (r, v) = (0, inf)
    for i in range(n - k):
        tmp = log2(_gaussian_elimination_complexity(n, k, i))
        if v > tmp and r < mem:
            r = i
            v = tmp
    return r


    
def gauss(n,k,w,p,l):
	r = _optimize_m4ri(n,k, 30)
	return  _gaussian_elimination_complexity(n, k, r)*n

delta = 0
@lru_cache(maxsize = 12800000)
def pr_dual_left(k, w, p, l):
	return (comb(floor((float)(k+l)/2)+delta ,p)*comb(k-floor((float)(k+l)/2)-delta, w//2 - p)) / comb(k,w//2)


@lru_cache(maxsize = 12800000)
def pr_dual_right(k, w, p, l):
	return (comb(ceil((float)(k+l)/2)-delta,p)*comb(k-ceil((float)(k+l)/2)+delta, w//2 - p)) / comb(k,w//2)

def pr_success(n,k,w,p,l,dual=0):
	pr = comb(n,w)/(comb(k+l,p)*comb(n-k-l, w - p)) 
	return log2(pr)

def opt_stern(M_UPPER_BOUND):


    """this function calculates the cost of Stern-Dumer in a similar way as the cost of Sieving ISD is computed in the paper
    
    not all operations are performed on vectors of size n, indeed most of the operations happen only on the small instance use sparse test patterns.

    The claim is NOT that the cost should be computed in this way.
    This is only for comparison.
    """

    time = 1e6
    r = _optimize_m4ri(n, k, M_UPPER_BOUND)
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
            tmp_mem = log2(2*L + _mem_matrix(n,k,r))
            if tmp_mem + log2(n) > M_UPPER_BOUND:
                continue
            if tmp < time: 
                time = tmp
                Memory = log2(2*L + _mem_matrix(n,k,r))
                params = {"p": p, "l" : l, "Memory" : Memory}
                # required memory: store list 1 in hashmap -> store index of test pattern, requires log2(L) bits
                # for ease of implmentation overdimension hashmap by factor 2
                mem = log2(L*log2(L) * 2)

            time = min(time, tmp)

    return time

def mmt_new(M_UPPER_BOUND):
	
	time = inf
	memory = 0
	params = {}
	r = _optimize_m4ri(n,k,M_UPPER_BOUND)
	Tg = _gaussian_elimination_complexity(n,k,r)*n
	
	for p in range(4, 60, 4):
		rep = binom(p, p/2) #The number of possible way to split the error vectors.
		l2 = int(floor(log2(rep))) # the number of bits to match in the middle layer
		
		for l1 in range(2, min(n-k - w + p - l2, 100)):
			l = l1 + l2
			
			#Top list, index set of size (k+l)/2, weight precisely p/4
			L2 = binom((k+l)/2, p/4)
			#Middle list
			L1 = L2**2 /2**l2
			#Bottom
			L = L1**2 / 2**l1
			
			#Memory, in terms of list size
			tmp_mem = (log2(2*L2 + L1 + _mem_matrix(n, k, r)))
			if tmp_mem + log2(n) > M_UPPER_BOUND:
				continue
			#complexity of the tree
			T2 = 4 * (L2) * p/4 * l2 	#compute l2 syndrome in 4 lists, each vector weight p/4
			T1 = 2 * L1 * p/2 * l1 		# for each matching (weight p/2) in the middle layer, compute l1 syndrome.
			T0 = L * p * (n-k-l) 		# For each solution candidate, check the weight condition.
			T_tree = T2 + T1 + T0
			
			##Split Probability : The way of constructing the top list lower the number of representation for each solution.
			P_split = binom(p/2, p/4)**2 / binom(p,p/2)
			#Probability
			Tp = binom((k+l)/2, p/2)**2 * binom(n-k-l, w - p) / binom(n,w)
			Tp = Tp * P_split * 1/2 # Last term: reps surviving probability for MMT.
			tmp_complexity = log2((Tg + T_tree)/Tp)
			if tmp_complexity < time:
					time = tmp_complexity
					#params = {"p": p, "l" : l, "l2": l2, "l1": l1, "memory (list size)" : tmp_mem}
	return time

def bjmm_depth_2_new(M_UPPER_BOUND):
	"""
	BJMM in depth 2 concrete complexity analysis, using for comparision with GJN algorithm.
	[BJMM12] Becker, A., Joux, A., May, A., Meurer, A.: Decoding random binary linear codes in 2^(n/20): How 1+ 1= 0
    	improves information set decoding. In: Annual international conference on the theory and applications of
    	cryptographic techniques. pp. 520–536. Springer (2012)


	"""
	#solutions = max(0, log2(binom(n, w)) - (n - k))
	time = inf
	memory = 0
	params = {}
	r = _optimize_m4ri(n, k, M_UPPER_BOUND)
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
				C2 = L2**2 // 2**l2 # == L1
				C1 = L1**2 // 2**l1
				#print(log2(C2), p, p1, l, l1, l2)
				
				#Memory, in terms of list size
				tmp_mem = log2(2*L2 + L1 + _mem_matrix(n, k, r))
				if tmp_mem + log2(n) > M_UPPER_BOUND:
					continue
				

				#complexity of the tree
				T2 = 4 * L2 * ((p1/2) * l2) 	#compute l2 syndrome bits
				T1 = 2*C2 * ((p1) * l1) 		#compute l1 syndrome bits for each matching		
				T0 = C1 * p * (n-k-l)		#checking solution
				T_tree = T2 + T1 + T0
				
				P_split = binom((k+l)/2, p1/2)**2 / binom(k+l, p1) #Splitting probability of middle layer (i.e, probability that the solution can be represented with the choice of the base lists)
				
				#Tp = P_split**2 * binom(k+l, p) * binom(n-k-l, w-p) / binom(n,w)
				Tp = 2*log2(P_split) + log2(binom(k+l,p)) + log2(binom(n-k-l, w-p)) - log2(binom(n,w)) + log2(1-exp(-2**(-l2)*rep)) # P_split^2 in depth-2, P_split^4 in depth 3, last term reps surviving probability.
				tmp_complexity = log2((Tg + T_tree)) - Tp
				
				if tmp_complexity < time:
					time = tmp_complexity
					params = {"p": p, "p1": p1, "eps" : eps, "l" : l, "l2": l2, "l1": l1, "P_split": 2*log2(P_split), "memory (list size)": tmp_mem}
	return time

def C_final_check(n,k,p,l,M, c_label = c_label):
	sol = comb(k+l,p)/(2**l) # number of solutions
	return p*(n-k-l)*sol

def C_check(n,k,p,l,M, c_label = c_label):
	return p*l*M  

def C_label(n,k,p,l,M, c_label = c_label):
	return c_label*comb(p, int(p/2)) *M *l

def C_move(n,k,p,l,M, p_prime, c_label = c_label):
	p_pprime = int(p/2) - p_prime
	return comb(p-p_pprime, p_prime) *M *l

def C_combine(n,k,p,l,M, c_label = c_label):
	return (2*p  + 4* log2(M))*M *l

def C_inner_prime(n,k,p,l,M, p_prime, c_label = c_label):
	p_pprime = int(p/2) - p_prime
	inner = l*M*(p + c_label*comb(p, int(p/2))  + comb(p - p_pprime, p_prime) + 2*p  + 4* log2(M))
	# p is from the check; c_label is label; C_mov is comb(p-p_pprime, p_prime); 2*p + 4* log2(M) is from combine. 
	final_check = C_final_check(n,k,p,l,M, c_label)
	return final_check + inner

def C_sd(n,k,p,l,M, p_prime, c_label = c_label):
	p_pprime = int(p/2) - p_prime
	inner = l*M*(p + c_label*comb(p, int(p/2))  + comb(p - p_pprime, p_prime) + 2*p  + 4* log2(M))
	# p is from the check; c_label is label; C_mov is comb(p-p_pprime, p_prime); 2*p + 4* log2(M) is from combine. 
	return inner


def ISD_prime(n,k,w,p,p_prime, l,M,dual):
	# print(k)
	return pr_success(n,k,w,p,l, dual) + log2(float(gauss(n,k,w,p,l) + C_inner_prime(n,k,p,l,M,p_prime)))
	

def  CM_attack_key(M_UPPER_BOUND):
		'''
		Attack complexity of CM key.
		'''
		init = 1000.0
		mem = 0
		for p in range (2, 28+1,2):
			for p_prime in range (0, int(p/2)):
				p_pprime = int(p/2) - p_prime
				for l in range (0,151):
					q = comb(p,int(p/2))*comb(k+l-p,int(p/2))/comb(k+l,p)
					M = 2/q
					M *= M_overhead
					sol = (comb(k+l,p))/(np.longdouble)(2**l)
					temp = ISD_prime(n,k,w,p,p_prime, l,M, dual=0)
					if sol >= 1 and temp <= init and M/2**(5)>= sol and comb(k+l, p_pprime) <= M and log2(M) + log2(n) <M_UPPER_BOUND: 
						# The constraint of M/2**(5)>= sol is added to ensure the probability is higher than (50%)
						init = temp
		return init


if __name__=='__main__':
	M_UPPER_BOUND = np.arange(40, 100, 1)
	
	
	y1 = np.vectorize(CM_attack_key)
	y2 = np.vectorize(mmt_new)
	y3 = np.vectorize(opt_stern)
	y4 = np.vectorize(bjmm_depth_2_new)
	#y5 = np.vectorize(both_may_depth_2_complexity)
	#y6 = np.vectorize(bjmm_depth_2_complexity)
	plt.plot(M_UPPER_BOUND, y1(M_UPPER_BOUND), linestyle = 'solid', label ='Our attack')
	plt.plot(M_UPPER_BOUND,y2(M_UPPER_BOUND), linestyle = 'dotted', label = 'MMT')
	plt.plot(M_UPPER_BOUND,y3(M_UPPER_BOUND), linestyle = 'dashed', label = 'Stern')
	plt.plot(M_UPPER_BOUND,y4(M_UPPER_BOUND), linestyle = 'dashdot',label = 'BJMM')
	#plt.plot(M_UPPER_BOUND,y5(M_UPPER_BOUND), linestyle = 'dotted',label = 'Both-May')
	#plt.plot(M_UPPER_BOUND,y6(M_UPPER_BOUND), linestyle = 'dashed',label = 'BJMM')
	
	
	
	plt.legend(loc="upper right")
	plt.xlabel(r'Memory restriction')
	plt.ylabel(r'Complexity')
	plt.savefig("time-memTO1.png")
	
	#M_UPPER_BOUND = inf
	#print(CM_attack_key(M_UPPER_BOUND))
	#print(n,k,w)		
