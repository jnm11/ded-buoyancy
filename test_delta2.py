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

g.append(de.Fourier('x', N, interval = (0,L), dealias=A))
g.append(de.SinCos( 'x', N, interval = (0,L), dealias=A))

for k in range(2):
    d.append(de.Domain([g[k]], np.float64))
    ff.append(d[k].new_field())
    ff[k].meta['x']['parity']=k
    fs.append(d[k].new_field())
    fs[k].meta['x']['parity']=k
    
w=np.abs([d[0].bases[0].wavenumbers,d[1].bases[0].wavenumbers])*L/math.pi
w[0]=w[0]/(N-2)
w[1]=w[1]/(N-1)

MM=np.arange(2,7)
e1=np.ndarray((2,))
e2=np.ndarray((2,))
e3=np.ndarray((2,))
YY=[0,0.1,0.3,0.5]
fg, ax=plt.subplots(nrows=len(YY),ncols=4)
fg, bx=plt.subplots(nrows=len(YY),ncols=4)


float_formatter = lambda x: '{:8.1e}'.format(x)
np.set_printoptions(formatter={'float_kind':float_formatter})
for k in range(len(YY)):
    y=YY[k]
    ax[k,0].set_ylabel('y {:4.2f}'.format(y))
    for j in range(2):
        x1 = g[j].grid(scale=1)
        x2 = g[j].grid(scale=A)
        P=j
        x1 = g[j].grid(scale=1)
        x2 = g[j].grid(scale=A)
        f=ff[j]
        h=fs[j]
        f.set_scales(1)
        h.set_scales(1)
        z,I=ju.ffta_delta2(x1,y,P,N,M)
        f['g']=z
        h['g']=z**2
        fc=copy.copy(np.abs(f['c']))
        hc=copy.copy(np.abs(h['c']))
        f.set_scales(A)
        h.set_scales(A)
        f1=copy.copy(f['g'])
        f2=copy.copy(h['g'])
        if   P==0:
            f1=np.concatenate((f1,f1,f1,f1))
            f2=np.concatenate((f2,f2,f2,f2))
        else:
            f1=np.concatenate((np.flip(f1),f1,np.flip(f1),f1))
            f2=np.concatenate((np.flip(f2),f2,np.flip(f2),f2))
        x2=np.concatenate((x2-L,x2,x2+L,x2+2*L))

        #print('j={}, w.shape {}, fc.shape {}, hc.shape {}'.format(j,w[j].shape,fc.shape,hc.shape))
        ax[k,j].plot(x2,f1)
        ax[k,j].set_xlim((-L,3*L))
        ax[k,j+2].plot(x2,f2)
        ax[k,j+2].set_xlim((-L,3*L))
        bx[k,j].plot(w[j],fc)
        bx[k,j].set_xlim((0,1))
        bx[k,j+2].plot(w[j],hc)
        bx[k,j+2].set_xlim((0,1))
        
        e3[j]=f1.sum()/NN/I/4-1
        e1[j]=max(f1.max()-1,-f1.min())
        e2[j]=max(f2.max()-1,-f2.min())
    print('M: {:2d}, I=[{:8.1e} {:8.1e}] e1=[{:8.1e} {:8.1e}] e2=[{:8.1e} {:8.1e}]'.format(M,e3[0],e3[1],e1[0],e1[1],e2[0],e2[1]))
    ax[0,0].set_title('P=0 f')
    ax[0,1].set_title('P=1 f')
    bx[0,0].set_title('P=0 fc')
    bx[0,1].set_title('P=1 fc')
    ax[0,2].set_title('P=0 f**2')
    ax[0,3].set_title('P=1 f**2')
    bx[0,2].set_title('P=0 f**2 c')
    bx[0,3].set_title('P=1 f**2 c')

plt.draw()
plt.waitforbuttonpress()
