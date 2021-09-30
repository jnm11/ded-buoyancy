#!/usr/bin/env python3
#import matplotlib
#matplotlib.use("TkAgg")
"""
Test spectral aliasing of front 

Usage: 
   test_alias.py [options]

Options:
  -L, --Lambda=<>  : Inverse time period    [default:  1]
  -S, --Sigma=<>   : Noise variance         [default:  1]
  -d, --dt=<>      : Time step              [default:  0.01]
  -n, --Num=<>     : Number of trajecties   [default:  1000]
  -T, --Time=<>    : Time to integrate      [default:  10]

"""

from docopt import docopt
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys
import jutil as ju
import logging
import math
import copy
plt.ion()

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

lmbd  = np.float64(args['--Lambda'])
Sigma = np.float64(args['--Sigma'])
dt    = np.float64(args['--dt'])
num   = np.int64(  args['--Num'])
Time  = np.float64(args['--Time'])

rand = np.random.RandomState()#(seed=42)

m=np.int64(Time/dt)
x=np.zeros([m,num])
x[0,:]=np.sqrt(Sigma)*rand.standard_normal([num])
(a,b) = ju.noisedt(dt,lmbd,np.sqrt(Sigma))
for j in range(m-1): x[j+1,:] = a*x[j,]+b*rand.standard_normal([num])

s=np.sum(x**2,axis=1)/num

c=np.zeros(m)
for j in range(m-1): c[j]=np.sum(x[0:m-j,:]*x[j:m,:])/(num*(m-j+1))

t=range(m)*dt   
fg, (ax1,ax2,ax3)=plt.subplots(nrows=3,ncols=1)
ax1.plot(t,s,color='r')
ax1.plot([0,Time],[Sigma,Sigma],color='b')
ax3.set_xlim((0,Time))
ax2.hist(x.flatten(),bins=100)
ax3.set_xlim((-4,4))
#ax3.hist(np.log10(x.flatten()**2/Sigma),bins=100)
#ax3.set_xlim((-5,5))
ax3.plot(t,c,color='r')
ax3.plot(t,Sigma*np.exp(-lmbd*t),color='b')
ax3.set_xlim((0,Time))

plt.draw()
plt.waitforbuttonpress()

