import numpy as np
from math import comb, sqrt,log2,ceil,floor, isnan, exp
from functools import cache


def h(x):
    if x==0.0 or x==1: 
        return 0
    else: 
        return -x*np.log2(x) - (1-x)*np.log2(1-x)

@cache
def com(x,y): #x>y
	c = 5.0 # coefficient for approximations of the binomial function.
	if x < y: return 0
	elif x==0 and y==0:
		return 0.0
	elif y == 1:
		return x
	elif y == 0:
		return 1
	elif y == x:
		return 1
	elif x < 10:
		return comb(x,y)
	elif x <= 10000:
		return com(x-1, y-1) + com(x-1,y)
	else:
		if c * y *(x-y) <= 0:
			print(y)
			print(x-y)
			print(c * y *(x-y))
		return 2**((h(float(y/x)))*x) * sqrt(x/(c * y *(x-y)))



def successProb_new(k,l,p):
	q = comb(p,int(p/2))*comb(k+l-p,int(p/2))/comb(k+l,p)
	M = 2/q
	delta = 2.0/3
	M /= delta
	print("Memory requirements:", log2(M))
	initial = comb(k+l, p)
	# success = (1-1/initial)**M
	success = exp(-M/initial)
	M_i = M
	for i in range (1, l+1):
		print("run: ", i)
		#2 times M_prime i because we consider the syndrome and 0.
		M_primei = 2*(comb(k+l,p)/(2**i))
		print("M'_i: ", M_primei)
		M_1i = M_i/2
		print("M_1i ", log2(M_1i))
		new = 9/16 * M_i**2/ M
		M_2i = min(new*(1- (M_1i + new/2)/M_primei),M/2)
		print("M_2i: ", log2(M_2i))
		M_i = min(M, M_1i + M_2i)
		if (M_i < 2*M_primei and M_2i > 1):
			# success *= (1-1/(M_primei))**(M_2i)
			success *= exp(-M_2i/M_primei)
		print("M_i: ", log2(M_i))
		print("Success probability:" ,1- success) 
	print("Success probability:" ,1- success) 

def testC():
	k= 1000
	l = 27
	p = 4
	successProb_new(k,l,p)

def testD(para):
	if para == 0:
		k= 300
		l = 28
		p = 6
	elif para == 1:
		k= 12323
		l = 48
		p = 6
	elif para == '1H':
		k= 17669
		l = 51
		p = 6
	elif para == '3B':
		k= 24659
		l = 52
		p = 6
	elif para == 2:
		k= 57637
		l = 57
		p = 6
	elif para == 3:
		k= 3488
		l = 83
		p = 14
	elif para == 4:
		k= 3360
		l = 96
		p = 16
	elif para == 5:
		k= 6528
		l = 92
		p = 14
	elif para == 6:
		k= 5024
		l = 101
		p = 16
	elif para == 7:
		k= 5413
		l = 102
		p = 16
	successProb_new(k,l,p)


def main():
	testC()
	testD(4)
	


if __name__ == '__main__':
	main()