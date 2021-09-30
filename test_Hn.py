#!/usr/bin/env python3
import numpy as np
from fractions import Fraction as frac
import math

def Hn(x,x1,x2,n):
    Pi=math.pi
    y=np.minimum(Pi,np.maximum(0,(x-x1)*Pi/(x2-x1)))
    f=y 
    print('n:{} '.format(n), end = '')
    for j in range(1,n+1):
        if j==1 : c = frac(-n,n+1)
        else:     c *= frac((j-1)*(j-1-n),j*(n+j))
        f += c*np.sin(2*j*y)
        print('{:10s}'.format('{}'.format(c)), end = '')
    print('')
    return(f/Pi)


x=0.4
x1=0.1
x2=0.8              
for n in range(1,10): Hn(x,x1,x2,n)
              
