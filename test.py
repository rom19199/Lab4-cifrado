# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 22:38:33 2021

@author: hugo_
"""

import math
import scipy.special as ss
from fractions import Fraction
import numpy as np

bit = "100010101010101001010101101011001010100111000111010001111000101010010010011101010101000010110001010100010101011010100001110000010101010100111100010"


def test(string, th=0.01):
    n=len(string)
    ones = string.count('1') #number of ones

    zeroes = string.count('0')    #number of zeros

    s = abs(ones - zeroes)  

    p = math.erfc(float(s)/(math.sqrt(float(n)) * math.sqrt(2.0))) #p-value

    success = ( p >= th)  # success = true if p-value >= 0.01

    return [ p, success]

# print("test1",test(bit))


def test2(string, th=0.01, M=32):
   
    n =len(string)
    N = int(math.floor(n/M))
    
    if N > 99:
        N=99
        M = int(math.floor(n/N))

    if n < 100:
        # Too little data for test. Input of length at least 100 bits required
        return [0.0, 0.0, False]

    num_of_blocks = N

    block_size = M 

    proportions = list()

    for i in range(num_of_blocks):
        
        block = string[i*(block_size):((i+1)*(block_size))]
        
        ones = block.count('1')

        zeroes = block.count('0') 
        
        proportions.append(Fraction(ones,block_size))

    chisq = 0.0

    for prop in proportions:
        chisq += 4.0*block_size*((prop - Fraction(1,2))**2)
    
    p = ss.gammaincc((num_of_blocks/2.0),float(chisq)/2.0) # p-value
    
    success = (p>= th)

    return [ p, success]

# print("test2", test2(bit))


def test3(cadena, th = 0.01 ):
    n =len(cadena)
    ones = cadena.count('1') #number of ones

    zeroes = cadena.count('0')    #number of zeros

    prop = float(ones)/float(n)

    tau = 2.0/math.sqrt(n)

    vobs = 0.0

    if abs(prop-0.5) > tau:
        p = 0
    else:

        vobs = 1.0
        for i in range(n-1):
            if cadena[i] != cadena[i+1]:
                vobs += 1.0

        p = math.erfc(abs(vobs - (2.0*n*prop*(1.0-prop)))/(2.0*math.sqrt(2.0*n)*prop*(1-prop) ))
    
    success = (p>=th)


    return [p, success]

# print("test3",test3(bit))



import copy
import gf2matrix

def test4(cadena, th = 0.01, M=32, Q=32):
    n = len(cadena)
    N = int(math.floor(n/(M*Q))) #Number of blocks
    
    if N < 38:
        print("  Number of blocks must be greater than 37")
        p = 0.0
        return [0]*9
        
    # Compute the reference probabilities for FM, FMM and remainder 
    r = M
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FR_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    r = M-1
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FRM1_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    LR_prob = 1.0 - (FR_prob + FRM1_prob)
    
    FM = 0      # Number of full rank matrices
    FMM = 0     # Number of rank -1 matrices
    remainder = 0
    for blknum in range(N):
        block = [None] * (M*Q)
        
        for i in range(M*Q):
            block[i] = int(cadena[blknum*M*Q + i],2)
            
        # Put in a matrix
        matrix = gf2matrix.matrix_from_bits(M,Q,block,blknum) 
        rank = gf2matrix.rank(M,Q,matrix,blknum)


        if rank == M: # count the result
            FM += 1
        elif rank == M-1:
            FMM += 1  
        else:
            remainder += 1

    chisq =  (((FM-(FR_prob*N))**2)/(FR_prob*N))
    chisq += (((FMM-(FRM1_prob*N))**2)/(FRM1_prob*N))
    chisq += (((remainder-(LR_prob*N))**2)/(LR_prob*N))
    
    p = math.e **(-chisq/2.0)

    success = (p >= th)

    return [ p, success]

# print("test4",test4(bit))


import scipy.special as ss
import random
def test5(cadena, th=0.01):
    n = len(cadena)

    # The templates provdided in SP800-22rev1a
    templates = [None for x in range(7)]
    templates[0] = [[0,1],[1,0]]
    templates[1] = [[0,0,1],[0,1,1],[1,0,0],[1,1,0]]
    templates[2] = [[0,0,0,1],[0,0,1,1],[0,1,1,1],[1,0,0,0],[1,1,0,0],[1,1,1,0]]
    templates[3] = [[0,0,0,0,1],[0,0,0,1,1],[0,0,1,0,1],[0,1,0,1,1],[0,0,1,1,1],[0,1,1,1,1],
                    [1,1,1,0,0],[1,1,0,1,0],[1,0,1,0,0],[1,1,0,0,0],[1,0,0,0,0],[1,1,1,1,0]]
    templates[4] = [[0,0,0,0,0,1],[0,0,0,0,1,1],[0,0,0,1,0,1],[0,0,0,1,1,1],[0,0,1,0,1,1],
                    [0,0,1,1,0,1],[0,0,1,1,1,1],[0,1,0,0,1,1],
                    [0,1,0,1,1,1],[0,1,1,1,1,1],[1,0,0,0,0,0],
                    [1,0,1,0,0,0],[1,0,1,1,0,0],[1,1,0,0,0,0],
                    [1,1,0,0,1,0],[1,1,0,1,0,0],[1,1,1,0,0,0],
                    [1,1,1,0,1,0],[1,1,1,1,0,0],[1,1,1,1,1,0]]
    templates[5] = [[0,0,0,0,0,0,1],[0,0,0,0,0,1,1],[0,0,0,0,1,0,1],[0,0,0,0,1,1,1],
                    [0,0,0,1,0,0,1],[0,0,0,1,0,1,1],[0,0,0,1,1,0,1],[0,0,0,1,1,1,1],
                    [0,0,1,0,0,1,1],[0,0,1,0,1,0,1],[0,0,1,0,1,1,1],[0,0,1,1,0,1,1],
                    [0,0,1,1,1,0,1],[0,0,1,1,1,1,1],[0,1,0,0,0,1,1],[0,1,0,0,1,1,1],
                    [0,1,0,1,0,1,1],[0,1,0,1,1,1,1],[0,1,1,0,1,1,1],[0,1,1,1,1,1,1],
                    [1,0,0,0,0,0,0],[1,0,0,1,0,0,0],[1,0,1,0,0,0,0],[1,0,1,0,1,0,0],
                    [1,0,1,1,0,0,0],[1,0,1,1,1,0,0],[1,1,0,0,0,0,0],[1,1,0,0,0,1,0],
                    [1,1,0,0,1,0,0],[1,1,0,1,0,0,0],[1,1,0,1,0,1,0],[1,1,0,1,1,0,0],
                    [1,1,1,0,0,0,0],[1,1,1,0,0,1,0],[1,1,1,0,1,0,0],[1,1,1,0,1,1,0],
                    [1,1,1,1,0,0,0],[1,1,1,1,0,1,0],[1,1,1,1,1,0,0],[1,1,1,1,1,1,0]]
    templates[6] = [[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,1],[0,0,0,0,0,1,0,1],[0,0,0,0,0,1,1,1],
                    [0,0,0,0,1,0,0,1],[0,0,0,0,1,0,1,1],[0,0,0,0,1,1,0,1],[0,0,0,0,1,1,1,1],
                    [0,0,0,1,0,0,1,1],[0,0,0,1,0,1,0,1],[0,0,0,1,0,1,1,1],[0,0,0,1,1,0,0,1],
                    [0,0,0,1,1,0,1,1],[0,0,0,1,1,1,0,1],[0,0,0,1,1,1,1,1],[0,0,1,0,0,0,1,1],
                    [0,0,1,0,0,1,0,1],[0,0,1,0,0,1,1,1],[0,0,1,0,1,0,1,1],[0,0,1,0,1,1,0,1],
                    [0,0,1,0,1,1,1,1],[0,0,1,1,0,1,0,1],[0,0,1,1,0,1,1,1],[0,0,1,1,1,0,1,1],
                    [0,0,1,1,1,1,0,1],[0,0,1,1,1,1,1,1],[0,1,0,0,0,0,1,1],[0,1,0,0,0,1,1,1],
                    [0,1,0,0,1,0,1,1],[0,1,0,0,1,1,1,1],[0,1,0,1,0,0,1,1],[0,1,0,1,0,1,1,1],
                    [0,1,0,1,1,0,1,1],[0,1,0,1,1,1,1,1],[0,1,1,0,0,1,1,1],[0,1,1,0,1,1,1,1],
                    [0,1,1,1,1,1,1,1],[1,0,0,0,0,0,0,0],[1,0,0,1,0,0,0,0],[1,0,0,1,1,0,0,0],
                    [1,0,1,0,0,0,0,0],[1,0,1,0,0,1,0,0],[1,0,1,0,1,0,0,0],[1,0,1,0,1,1,0,0],
                    [1,0,1,1,0,0,0,0],[1,0,1,1,0,1,0,0],[1,0,1,1,1,0,0,0],[1,0,1,1,1,1,0,0],
                    [1,1,0,0,0,0,0,0],[1,1,0,0,0,0,1,0],[1,1,0,0,0,1,0,0],[1,1,0,0,1,0,0,0],
                    [1,1,0,0,1,0,1,0],[1,1,0,1,0,0,0,0],[1,1,0,1,0,0,1,0],[1,1,0,1,0,1,0,0],
                    [1,1,0,1,1,0,0,0],[1,1,0,1,1,0,1,0],[1,1,0,1,1,1,0,0],[1,1,1,0,0,0,0,0],
                    [1,1,1,0,0,0,1,0],[1,1,1,0,0,1,0,0],[1,1,1,0,0,1,1,0],[1,1,1,0,1,0,0,0],
                    [1,1,1,0,1,0,1,0],[1,1,1,0,1,1,0,0],[1,1,1,1,0,0,0,0],[1,1,1,1,0,0,1,0],
                    [1,1,1,1,0,1,0,0],[1,1,1,1,0,1,1,0],[1,1,1,1,1,0,0,0],[1,1,1,1,1,0,1,0],
                    [1,1,1,1,1,1,0,0],[1,1,1,1,1,1,1,0]]
    
    # Randomly choose the template B 
    r = random.SystemRandom()
    template_list = r.choice(templates)
    B = r.choice(template_list)
    
    m = len(B)
    
    N = 8  #number of block
    M = int(n/N) #length of each block
    
    blocks = list() # Split into N blocks of M bits
    for i in range(N):
        block = list()
        for j in range(M):
            block.append(int(cadena[i*M+j],2))
        blocks.append(block)

    W=list() # Count the number of matches of the template in each block Wj
    for block in blocks:
        position = 0
        count = 0
        while position < (M-m):

            if block[position:position+m] == B:
                position += m
                count += 1
            else:
                position += 1
        W.append(count)

    mu = float(M-m+1)/float(2**m) # Compute mu and sigma
    sigma = M * ((1.0/float(2**m))-(float((2*m)-1)/float(2**(2*m))))

    chi_sq = 0.0  # Compute Chi-Square
    for j in range(N):
        chi_sq += ((W[j] - mu)**2)/sigma

    p = ss.gammaincc(N/2.0, chi_sq/2.0) # Compute P value

    success = ( p >= th)

    return [ p, success]



import scipy.special as ss

def padding(string, n):
	while len(string) <n:
		string = '0' + string
	return input

def berlekamp_massey(string):
    n = len(string)

    b = '1' + '0'*(n-1)
    c = '1' + '0'*(n-1)

    L = 0
    m = -1
    N = 0
    while (N < n):
        #compute discrepancy
        d = int(string[N],2)
        if L>0:
            k = c[1:L+1]
            h = string[N-L:N][::-1]

            k = int(k,2)  #binary to integer

            h = int(h,2)    #binary to integer

            r = k&h #bitwise-and

            r = bin(r)[2:] 

            r = r.count('1')

            d = d ^ (r%2)

        if d != 0:  # If d is not zero, adjust poly
            t = c
            k = c[N-m:n]
            k = int(k, 2)
            h = b[0:n-N+m]
            h = int(h,2)
            k = k^h
            c = c[0:N-m] + padding(bin(k)[2:], n-N+m)
            # print(c)
            if (L <= (N/2)):
                L = N + 1 - L
                m = N
                b = t 
        N = N +1
    # Return length of generator and the polynomial
    return L , c[0:L]
    

def int2patt(n,m):
    pattern = list()
    for i in range(m):
        pattern.append((n >> i) & 1)
    return pattern
    
def countpattern(patt,padded_input,n):
    thecount = 0
    for i in range(n):
        match = True
        for j in range(len(patt)):
            if str(patt[j]) != padded_input[i+j]:
                match = False
                break
        if match:
            thecount += 1
    return thecount

def psi_sq_mv1(m, n, padded_input):
    counts = [0 for i in range(2**m)]
    for i in range(2**m):

        pattern = int2patt(i,m)
        count = countpattern(pattern,padded_input,n)
        counts.append(count)

        # pattern = padding(bin(i)[2:], m)


        # count = padded_input.count(pattern)
        
        # counts.append(count)
        
    psi_sq_m = 0.0
    for count in counts: 
        psi_sq_m += (count**2)
    psi_sq_m = psi_sq_m * (2**m)/n 
    psi_sq_m -= n
    return psi_sq_m            
         
def test6(string, th = 0.01, patternlen=None):
    n = len(string)
    # Pattern length
    if patternlen != None:
        m = patternlen  
    else:  
        m = int(math.floor(math.log(n,2)))-2
    
        if m < 4:
            print("Error. Not enough data for m to be 4")
            return [0]*8
        m = 4

    # Step 1
    padded_input=string[0:n]+string[0:m-1]
    
    # Step 2
    psi_sq_m   = psi_sq_mv1(m, n, padded_input)
    psi_sq_mm1 = psi_sq_mv1(m-1, n, padded_input)
    psi_sq_mm2 = psi_sq_mv1(m-2, n, padded_input)    
    
    delta1 = psi_sq_m - psi_sq_mm1
    delta2 = psi_sq_m - (2*psi_sq_mm1) + psi_sq_mm2

    p1 = ss.gammaincc(2**(m-2),delta1/2.0)
    p2 = ss.gammaincc(2**(m-3),delta2/2.0)
     
    success = (p1 >= 0.01) and (p2 >= 0.01)

    return [p1, p2, (p1+p2)/2, success]



def LGC(a,b,mod,k=1,shape=1):
    
    s = ''
    x = np.random.choice(mod)
    for i in range(0,shape//k):
        x = (a*x + b) % mod
        bts = '{0:08b}'.format(x)
        s = s + bts[-k:]
    return s

mod = 2**15 - 1
a = 23550
b = 12964
n = 40000
s = LGC(a,b,mod,8,n)
print(len(s))

# tests = [test,test2,test3,test4,test5,test6,test7,test8,test9,test10]
tests = [test,test2,test3,test4,test5,test6]


results = lambda x: 'success' if x else 'fail'
for i in range(0,len(tests)):
    p,result = tests[i](s)
    print('test{}: p = {},{}'.format(str(i+1).zfill(2),p,results(result)))
    
    
n = 200
frequ = np.zeros((n,len(tests))).astype(int)

for i in range(0,n):
    s = LGC(a,b,mod,8,n)
    for j in range(0,len(tests)):
        _,res = tests[j](s,th=0.01)
        frequ[i,j]  = int(res)
        

    