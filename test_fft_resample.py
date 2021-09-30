#!/usr/bin/env python3
#import matplotlib
#matplotlib.use("TkAgg")
# test_alias
"""
Test spectral aliasing of front 

Usage: 
   test_alias.py [options]

Options:
  -L, --Length=<>  : Length of domain       [default:  1]
  -N, --Number=<>  : Number of grid points  [default:  32]
  -W, --Width=<>   : Width of front region  [default:  0.01]
  -O, --Order=<>   : Polynomial order       [default:  4]
  -F, --Front=<>   : Front position         [default:  2]
  -A, --Alias=<>   : Anti-aliasing factor   [default:  50]
      --y1=<>      : y1                     [default:  0.2]
      --y2=<>      : y2                     [default:  0.4]

"""
#import matplotlib
#matplotlib.use("TkAgg")

from docopt import docopt
from mpi4py import MPI
import dedalus.public as de
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys
import jutil as ju
import logging
import math
import copy
import scipy.special as sp
import logging

plt.ion()

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

ju.logger=logging.getLogger(__name__)

L=np.float64(args['--Length'])
N=np.int64(  args['--Number'])
W=np.float64(args['--Width'])
M=np.int64(  args['--Order'])
x=np.float64(args['--Front'])
A=np.float64(args['--Alias'])
y1=np.float64(args['--y1'])
y2=np.float64(args['--y2'])

NN=np.int64(A*N)

d=list()
ff=list()
fs=list()
g=list()
w=list()



M=100
N1=np.int64(M-2)
N2=np.int64(M*3)
x  = (np.arange(M)+0.5)/M*2*math.pi
x1 = (np.arange(N1)+0.5)/N1*2*math.pi
x2 = (np.arange(N2)+0.5)/N2*2*math.pi
y  = (np.arange(M))/M*2*math.pi
y1 = (np.arange(N1))/N1*2*math.pi
y2 = (np.arange(N2))/N2*2*math.pi

for n in range(3):
    s0=np.sin((2*n+1)*x/2)
    c0=np.cos((2*n+1)*x/2)
    f0=np.sin((n+1)*y+0.2)
    
    s0a=ju.fft_resample(s0,N1,-1)
    s0b=ju.fft_resample(s0,N2,-1)
    s0c=ju.fft_resample(s0, M,-1)
    f0a=ju.fft_resample(f0,N1,0)
    f0b=ju.fft_resample(f0,N2,0)
    f0c=ju.fft_resample(f0, M,0)
    c0a=ju.fft_resample(c0,N1,1)
    c0b=ju.fft_resample(c0,N2,1)
    c0c=ju.fft_resample(c0, M,1)

    fg, a=plt.subplots(nrows=3,ncols=1)                    
    
    a[0].plot(x,s0,x1,s0a,x2,s0b,x,s0c)
    a[1].plot(y,f0,y1,f0a,y2,f0b,y,f0c)
    a[2].plot(x,c0,x1,c0a,x2,c0b,x,c0c)
    
    print('s0.max  {:7.3f}'.format(np.abs(s0).max()))
    print('s0a.max {:7.3f}'.format(np.abs(s0a).max()))
    print('s0b.max {:7.3f}'.format(np.abs(s0b).max()))
    print('s0c.max {:7.3f}'.format(np.abs(s0c).max()))
    print('c0.max  {:7.3f}'.format(np.abs(c0).max()))
    print('c0a.max {:7.3f}'.format(np.abs(c0a).max()))
    print('c0b.max {:7.3f}'.format(np.abs(c0b).max()))
    print('c0c.max {:7.3f}'.format(np.abs(c0c).max()))
    
plt.draw()
plt.waitforbuttonpress()
quit(1)
