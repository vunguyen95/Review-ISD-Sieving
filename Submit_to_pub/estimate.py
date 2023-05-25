# import libraries
import mpmath as math
import numpy as np
from math import inf, floor, ceil
from functools import lru_cache, cache
import sys

DOOM = 1 # 0 for the normal attack; 1 for the doom attack.




# parameters
c_label = int(sys.argv[1])
M_UPPER_BOUND = float(sys.argv[2])
M_overhead = 1.5
print("C_label:",c_label)
print("M_upper_bound:",M_UPPER_BOUND)
print("DOOM:", DOOM)
print("  ")
print("================================================")

@lru_cache(maxsize = 12800000)
def comb(N,k):
    val= math.factorial(N)/(math.factorial(k)*math.factorial(N-k))
    return val

def log2(N):
    return math.log(N,2)


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



def C_final_check_doom(n,k,p,l,M, n_0, c_label = c_label):
	sol = n_0*comb(k+l,p)/(2**l) # number of solutions
	return p*(n-k-l)*sol


@cache
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

def ISD_doom(n,k,w,p,p_prime, l,M,dual, n_0):
	'''
	@summary: ISD_doom is the complexity for the ISD for the doom attack. (decoding out of the many syndromes)
	@param n_0: n_0 is the number of solutions in the doom attack.
	'''
	def C_final_check_doom(n,k,p,l,M, n_0, c_label = c_label):
		sol = n_0*comb(k+l,p)/(2**l) # number of solutions
		return p*(n-k-l)*sol
	def C_inner_doom(n,k,p,l,M, p_prime, n_0, c_label = c_label):
		p_pprime = int(p/2) - p_prime
		inner = l*M*(p + c_label*comb(p, int(p/2))  + comb(p - p_pprime, p_prime) + 2*p  + 4* log2(M))
		# p is from the check; c_label is label; C_mov is comb(p-p_pprime, p_prime); 2*p + 4* log2(M) is from combine. 
		final_check = C_final_check_doom(n,k,p,l,M,n_0, c_label)
		return final_check + inner
	# print(k)
	return pr_success(n,k,w,p,l, dual)- log2(n_0) + log2(float(gauss(n+n_0,k,w,p,l) + C_inner_doom(n,k,p,l,M,p_prime,n_0)))

class BIKE_security:
	"""
	_summary_
	Estimate the security of the BIKE proposal against the key-recovery attacks.
	"""
	def __init__(self, level):
		if level == 1: 
			self.k = 12323
			self.n = self.k * 2
			self.w = 142
			self.t = 134
			self.lev = 1
		if level == 3:
			self.k = 24659
			self.n = self.k * 2
			self.w = 206
			self.t = 199
			self.lev = 3
		if level == 5:
			self.k = 40973
			self.n = self.k * 2
			self.w = 274
			self.t = 264
			self.lev = 5

	def BIKE_attack_key(self):
		'''
		Attack complexity of BIKE
		'''
		params = [0,0]
		init = 1000.0
		mem = 0
		params.append(init)
		params.append(mem)
		params.append(0)
		for p in range (2, 10+1,2):
			for p_prime in range (0, int(p/2)):
				p_pprime = int(p/2) - p_prime
				for l in range (0,101):
					q = comb(p,int(p/2))*comb(self.k+l-p,int(p/2))/comb(self.k+l,p)
					M = 2/q
					M *= M_overhead
					sol = (comb(self.k+l,p))/(np.longdouble)(2**l)
					temp = ISD_prime(self.n,self.k,self.w,p,p_prime, l,M, dual=0) - log2(self.k) # cyclic shift of bike key.
					if sol >= 1 and temp <= init and M/2>= sol and comb(self.k+l, p_pprime) <= M and log2(M) <120:
						init = temp
						params[0] = p
						params[1]= l
						params[2] = temp
						params[3] = log2(M)
						params[4] = p_prime
					
		print(params)
		print("Gauss cost:", log2(gauss(self.n,self.k,self.w,params[0],params[1])))
		print("Merge set cost:", log2(C_sd(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("check cost:", log2(C_check(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_label cost:", log2(C_label(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_move cost:", log2(C_move(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("C_combine cost:", log2(C_combine(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("final check cost:", log2(C_final_check(self.n,self.k,params[0],params[1],2**params[3], params[4])))

	def BIKE_attack_message(self):
		'''
		Attack complexity of BIKE
		'''
		params = [0,0]
		init = 1000.0
		mem = 0
		params.append(init)
		params.append(mem)
		params.append(0)
		for p in range (2, 10+1,2):
			for p_prime in range (0, int(p/2)):
				p_pprime = int(p/2) - p_prime
				for l in range (0,101):
					q = comb(p,int(p/2))*comb(self.k+l-p,int(p/2))/comb(self.k+l,p)
					M = 2/q
					M *= M_overhead					
					if DOOM == 1: 
						M *= 2 # from the doom attack
						sol = self.k * (comb(self.k+l,p))/(np.longdouble)(2**l)
						temp = ISD_doom(self.n,self.k,self.t,p,p_prime, l,M,0, self.k) # out of many.
					if DOOM == 0:
						sol =  (comb(self.k+l,p))/(np.longdouble)(2**l)
						temp = ISD_prime(self.n,self.k,self.t,p,p_prime, l,M,0)
					if sol >= 1 and temp <= init and M/2>= sol and comb(self.k+l, p_pprime) <= M and log2(M) <M_UPPER_BOUND:
						init = temp
						params[0] = p
						params[1]= l
						params[2] = temp
						params[3] = log2(M)
						params[4] = p_prime
					
		print(params)
		print("Gauss cost:", log2(gauss(self.n,self.k,self.t,params[0],params[1])))
		print("Merge set cost:", log2(C_sd(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("check cost:", log2(C_check(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_label cost:", log2(C_label(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_move cost:", log2(C_move(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("C_combine cost:", log2(C_combine(self.n,self.k,params[0],params[1],2**params[3])))
		print("final check cost:", log2(C_final_check_doom(self.n,self.k,params[0],params[1],2**params[3], self.k)))
		print("probability:", pr_success(self.n,self.k,self.t,params[0],params[1], dual=0)-log2(self.k))


class HQC_security:
	"""
	_summary_
	Estimate the security of the HQC proposal against the key-recovery attacks.
	"""
	def __init__(self, level):
		if level == 1: 
			self.k = 17669
			self.n = self.k * 2
			self.w = 132
			self.lev = 1
		if level == 3:
			self.k = 35851
			self.n = self.k * 2
			self.w = 200
			self.lev = 3
		if level == 5:
			self.k = 57637
			self.n = self.k * 2
			self.w = 262
			self.lev = 5

	def HQC_attack_key(self):
		'''
		Attack complexity of HQC key.
		'''
		params = [0,0]
		init = 1000.0
		mem = 0
		params.append(init)
		params.append(mem)
		params.append(0)
		for p in range (2, 12+1,2):
			for p_prime in range (0, int(p/2)):
				p_pprime = int(p/2) - p_prime
				for l in range (0,101):
					q = comb(p,int(p/2))*comb(self.k+l-p,int(p/2))/comb(self.k+l,p)
					M = 2/q
					M *= M_overhead
					if DOOM == 1:
						M *= 2 # from the doom attack 
						sol = self.k * (comb(self.k+l,p))/(np.longdouble)(2**l)
						temp = ISD_doom(self.n,self.k,self.w,p,p_prime, l,M,0, self.k) # out of many.
					if DOOM == 0:
						sol =  (comb(self.k+l,p))/(np.longdouble)(2**l)
						temp = ISD_prime(self.n,self.k,self.w,p,p_prime, l,M,0)
					if sol >= 1 and temp <= init and M/2>= sol and comb(self.k+l, p_pprime) <= M and log2(M) <M_UPPER_BOUND:
						init = temp
						params[0] = p
						params[1]= l
						params[2] = temp
						params[3] = log2(M)
						params[4] = p_prime
					
		print(params)
		print("Gauss cost:", log2(gauss(self.n,self.k,self.w,params[0],params[1])))
		print("Merge set cost:", log2(C_sd(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("check cost:", log2(C_check(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_label cost:", log2(C_label(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_move cost:", log2(C_move(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("C_combine cost:", log2(C_combine(self.n,self.k,params[0],params[1],2**params[3])))
		print("final check cost:", log2(C_final_check_doom(self.n,self.k,params[0],params[1],2**params[3], self.k)))
		print("probability:", pr_success(self.n,self.k,self.w,params[0],params[1], dual=0)-log2(self.k))


class CM_security:
	"""
	_summary_
	Estimate the security of the CM proposal against the key-recovery attacks.
	"""
	def __init__(self, level):
		if level == 1: 
			self.k = 2720
			self.n = 3488
			self.w = 64
			self.lev = 1
		if level == 3:
			self.k = 3360
			self.n = 4608
			self.w = 96
			self.lev = 3
		if level == 5:
			self.k = 5024
			self.n = 6688
			self.w = 128
			self.lev = 5
		if level == 7:
			self.k = 5413
			self.n = 6960
			self.w = 119
			self.lev = 5
		if level == 9:
			self.k = 6528
			self.n = 8192
			self.w = 128
			self.lev = 5

	def CM_attack_key(self):
		'''
		Attack complexity of CM key.
		'''
		params = [0,0]
		init = 1000.0
		mem = 0
		params.append(init)
		params.append(mem)
		params.append(0)
		for p in range (2, 28+1,2):
			for p_prime in range (0, int(p/2)):
				p_pprime = int(p/2) - p_prime
				for l in range (0,151):
					q = comb(p,int(p/2))*comb(self.k+l-p,int(p/2))/comb(self.k+l,p)
					M = 2/q
					M *= M_overhead
					sol = (comb(self.k+l,p))/(np.longdouble)(2**l)
					temp = ISD_prime(self.n,self.k,self.w,p,p_prime, l,M, dual=0)
					if sol >= 1 and temp <= init and M/2**(5)>= sol and comb(self.k+l, p_pprime) <= M and log2(M) <M_UPPER_BOUND: 
						# The constraint of M/2**(5)>= sol is added to ensure the probability is higher than (50%)
						init = temp
						params[0] = p
						params[1]= l
						params[2] = temp
						params[3] = log2(M)
						params[4] = p_prime
					
		print(params)
		print("Gauss cost:", log2(gauss(self.n,self.k,self.w,params[0],params[1])))
		print("Merge set cost:", log2(C_sd(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("check cost:", log2(C_check(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_label cost:", log2(C_label(self.n,self.k,params[0],params[1],2**params[3])))
		print("C_move cost:", log2(C_move(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("C_combine cost:", log2(C_combine(self.n,self.k,params[0],params[1],2**params[3], params[4])))
		print("final check cost:", log2(C_final_check(self.n,self.k,params[0],params[1],2**params[3], params[4])))

def outputAttackBikeHqc():	
	for level in range(1,6,2):
		print("level:", level)
		print(" 	")
		print("BIKE key recovery:")
		bikePara = BIKE_security(level)
		print("n,k,w is:", bikePara.n, bikePara.k, bikePara.w)
		bikePara.BIKE_attack_key()
		print(" ")
		print("BIKE message recovery attacks:")
		print("n,k,w is:", bikePara.n, bikePara.k, bikePara.t)
		bikePara.BIKE_attack_message()
		print(" ")
		print("----------------------------------------------------------------")
		print(" ")
		print("HQC key recovery attacks:")
		hqcPara = HQC_security(level)
		print("n,k,w is:", hqcPara.n, hqcPara.k, hqcPara.w)
		hqcPara.HQC_attack_key()
		print(" ")
		print("=====================================================")	

def outputAttackCM():	
	for level in range(1,10,2):
		print("level:", level)
		print("------------------")
		print(" ")
		print("CM key recovery attacks:")
		CM_Para = CM_security(level)
		print("n,k,w is:",CM_Para.n, CM_Para.k, CM_Para.w)
		CM_Para.CM_attack_key()
		print(" ")
		print("=====================================================")	




if __name__=='__main__':
	math.mp.dps = 20
	outputAttackBikeHqc()
	outputAttackCM()

	




	
