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




import scipy.special as ss

def test4(string, th = 0.01):
    n = len(string)
    M8 = [0.2148, 0.3672, 0.2305, 0.1875]

    # Length of blocks
    M = 8 
                
    K = 3

    N = 16
            
    # Table of frequencies
    v = [0,0,0,0,0,0,0]

    for i in range(N): # over each block
        #find the longest run
        block = string[i*M:((i+1)*M)] # Block i
        
        run = 0
        longest = 0
        for j in range(M): # Count the bits.
            if block[j] == '1':
                run += 1
                if run > longest:
                    longest = run
            else:
                run = 0

        if longest <= 1:    v[0] += 1
        elif longest == 2:  v[1] += 1
        elif longest == 3:  v[2] += 1
        else:               v[3] += 1
    
    # Compute Chi-Sq
    chi_sq = 0.0
    for i in range(K+1):
        p_i = M8[i]
        upper = (v[i] - N*p_i)**2
        lower = N*p_i
        chi_sq += upper/lower
    # p-value
    p = ss.gammaincc(K/2.0, chi_sq/2.0)
    
    success = (p>=th)

    return [p, success]
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

    return [p1, success]


def test7(string, th = 0.01):
    n = len(string)
    m = int(math.floor(math.log(n,2)))-6
    if m < 2:
        m = 2
    if m >3 :
        m = 3
    
    Cmi = list()
    phi_m = list()
    for iterm in range(m,m+2):
        # Step 1 
        padded_input=string+string[0:iterm-1]
    
        # Step 2
        counts = list()
        for i in range(2**iterm):
            count = 0
            for j in range(n):
                if int(padded_input[j:j+iterm],2) == i:
                    count += 1
            counts.append(count)
    
        # step 3
        Ci = list()
        for i in range(2**iterm):
            Ci.append(float(counts[i])/float(n))
        
        Cmi.append(Ci)
    
        # Step 4
        sum = 0.0
        for i in range(2**iterm):
            if (Ci[i] > 0.0):
                sum += Ci[i]*math.log((Ci[i]/10.0))
        phi_m.append(sum)
        
    # Step 5 - let the loop steps 1-4 complete
    
    # Step 6
    appen_m = phi_m[0] - phi_m[1]

    chisq = 2*n*(math.log(2) - appen_m)

    # Step 7
    p = ss.gammaincc(2**(m-1),(chisq/2.0))
    
    success = (p >= th)

    return [ p, success]

def normcdf(n):
    return 0.5 * math.erfc(-n * math.sqrt(0.5))

def p_value(n,z):
    sum_a = 0.0
    startk = int(math.floor((((float(-n)/z)+1.0)/4.0)))
    endk   = int(math.floor((((float(n)/z)-1.0)/4.0)))
    for k in range(startk,endk+1):
        c = (((4.0*k)+1.0)*z)/math.sqrt(n)
        #d = scipy.stats.norm.cdf(c)
        d = normcdf(c)
        c = (((4.0*k)-1.0)*z)/math.sqrt(n)
        #e = scipy.stats.norm.cdf(c)
        e = normcdf(c)
        sum_a = sum_a + d - e

    sum_b = 0.0
    startk = int(math.floor((((float(-n)/z)-3.0)/4.0)))
    endk   = int(math.floor((((float(n)/z)-1.0)/4.0)))
    for k in range(startk,endk+1):
        c = (((4.0*k)+3.0)*z)/math.sqrt(n)
        #d = scipy.stats.norm.cdf(c)
        d = normcdf(c)
        c = (((4.0*k)+1.0)*z)/math.sqrt(n)
        #e = scipy.stats.norm.cdf(c)
        e = normcdf(c)
        sum_b = sum_b + d - e 

    p = 1.0 - sum_a + sum_b
    return p
    
def test8(string, th = 0.01):
    n = len(string)
    # Step 1
    x = list()             # Convert to +1,-1
    for i in range(n):
        #if bit == 0:
        x.append(int(string[i],2)*2-1)
        
    # Steps 2 and 3 Combined
    # Compute the partial sum and records the largest excursion.
    pos = 0
    forward_max = 0
    for e in x:
        pos = pos+e
        if abs(pos) > forward_max:
            forward_max = abs(pos)
    pos = 0
    backward_max = 0
    for e in reversed(x):
        pos = pos+e
        if abs(pos) > backward_max:
            backward_max = abs(pos)
     
    # Step 4
    p_forward  = p_value(n, forward_max)
    p_backward = p_value(n,backward_max)
    
    success = ((p_forward >= 0.01) and (p_backward >= 0.01))

    return [p_backward, success]

def test9(string, th=0.01):
    n = len(string)
    # Convert to +1,-1
    x = list()
    for i in range(n):
        x.append(int(string[i],2)*2 -1 )

    # Build the partial sums
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)    
    sprime = [0]+s+[0] # Add 0 on each end
    

    # Build the list of cycles
    pos = 1
    cycles = list()
    while (pos < len(sprime)):
        cycle = list()
        cycle.append(0)
        while sprime[pos]!=0:
            cycle.append(sprime[pos])
            pos += 1
        cycle.append(0)
        cycles.append(cycle)
        pos = pos + 1
    
    J = len(cycles)  
    
    vxk = [['a','b','c','d','e','f'] for y in [-4,-3,-2,-1,1,2,3,4] ]

    # Count Occurances  
    for k in range(6):
        for index in range(8):
            mapping = [-4,-3,-2,-1,1,2,3,4]
            x = mapping[index]
            cyclecount = 0
            #count how many cycles in which x occurs k times
            for cycle in cycles:
                oc = 0
                #Count how many times x occurs in the current cycle
                for pos in cycle:
                    if (pos == x):
                        oc += 1
                # If x occurs k times, increment the cycle count
                if (k < 5):
                    if oc == k:
                        cyclecount += 1
                else:
                    if k == 5:
                        if oc >=5:
                            cyclecount += 1
            vxk[index][k] = cyclecount
    
    # Table for reference random probabilities 
    pikx=[[0.5     ,0.25   ,0.125  ,0.0625  ,0.0312 ,0.0312],
          [0.75    ,0.0625 ,0.0469 ,0.0352  ,0.0264 ,0.0791],
          [0.8333  ,0.0278 ,0.0231 ,0.0193  ,0.0161 ,0.0804],
          [0.875   ,0.0156 ,0.0137 ,0.012   ,0.0105 ,0.0733],
          [0.9     ,0.01   ,0.009  ,0.0081  ,0.0073 ,0.0656],
          [0.9167  ,0.0069 ,0.0064 ,0.0058  ,0.0053 ,0.0588],
          [0.9286  ,0.0051 ,0.0047 ,0.0044  ,0.0041 ,0.0531]]
    
    success = True
    plist = list()
    chi_sq = list()
    p_total = 0.0
    for index in range(8):
        #list of states
        mapping = [-4,-3,-2,-1,1,2,3,4] 
        x = mapping[index]
        chisq = 0.0
        for k in range(6):
            top = float(vxk[index][k]) - (float(J) * (pikx[abs(x)-1][k]))
            top = top*top
            bottom = J * pikx[abs(x)-1][k]
            chisq += top/bottom

        p = ss.gammaincc(5.0/2.0,chisq/2.0)
        p_total += p
        plist.append(p)
        chi_sq.append(chisq)
        if p < 0.01:
            success = False

    return [ p_total/8, success]

def test10(string, th = 0.01):
    n = len(string)

    x = list()             # Convert to +1,-2
    for i in range(n):
        x.append(int(string[i],2)*2-1)

    # Build the partial sums
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)  
    # print(s)  
    sprime = [0]+s+[0] # Add 0 on each end

    # Count the number of cycles J
    J = 0
    for value in sprime[1:]:
        if value == 0:
            J += 1
            
    # Build the counts of offsets
    count = [0 for x in range(-9,10)]
    for value in sprime:
        if (abs(value) < 10):
            count[value] += 1

    # Compute P values
    success = True
    plist = list() # list of p-values for each state
    p_average = 0.0
    for x in range(-9,10): 
        if x != 0:
            top = abs(count[x]-J)
            bottom = math.sqrt(2.0 * J *((4.0*abs(x))-2.0))
            p = ss.erfc(top/bottom)

            # print("p[" + str(x) +"] = " + str(p))

            p_average +=p
            plist.append(p)
            if p < 0.01:
                success = False

    return [p_average/19, success]




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
tests = [test,test2,test3,test4,test5,test6,test7,test8,test9,test10]


results = lambda x: 'success' if x else 'fail'
for i in range(0,len(tests)):
    p,result = tests[i](s, th = 0.01)
    print('test {}: p = {}, {}'.format(str(i+1).zfill(2),p,results(result)))
    
    
n = 200
frequ = np.zeros((n,len(tests))).astype(int)

for i in range(0,n):
    s = LGC(a,b,mod,8,n)
    for j in range(0,len(tests)):
        _,res = tests[j](s,th=0.01)
        frequ[i,j]  = int(res)
        

    