#!/usr/bin/env python3
#import matplotlib
#matplotlib.use("TkAgg")
# test_alias
"""
Test spectral aliasing of front 

Usage: 
   test_alias.py [options]

Options:
  -L, --Length=<>  : Length of domain       [default:  4]
  -N, --Number=<>  : Number of grid points  [default:  400]
  -W, --Width=<>   : Width of front region  [default:  0.01]
  -O, --Order=<>   : Polynomial order       [default:  10]
  -F, --Front=<>   : Front position         [default:  2]
  -A, --Alias=<>   : Anti-aliasing factor   [default:  10]

"""

from docopt import docopt
from mpi4py import MPI
import dedalus.public as de
#import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys
import jutil as ju
import logging
import math
import copy
import scipy.special as sp

#plt.ion()

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

L=np.float64(args['--Length'])
N=np.int64(  args['--Number'])
W=np.float64(args['--Width'])
O=np.int64(  args['--Order'])
x=np.float64(args['--Front'])
A=np.float64(args['--Alias'])

#xb = de.Fourier('x', N, interval = (0,L), dealias=A)
xb = de.SinCos( 'x', N, interval = (0,L), dealias=A)
domain  = de.Domain([xb], np.float64)
x1      = xb.grid(scale=1)
x2      = xb.grid(scale=A)
ff      = domain.new_field()
ff.meta['x']['parity']=1
dx=x1[1]-x1[0]

nW=50
WW=np.power(10,np.linspace(0, 0.5, num=nW))
WW=np.linspace(2.8, 3.2, num=nW)

def ftanh(x,W):
    return((1-np.tanh(4*x/W))/2)
           
def ferf(x,W):
    return((1-sp.erf(2*np.sqrt(math.pi)*x/W))/2)

def fsmooth(x,W):
    return(1-ju.smoothwf(x,O,-W/2,W/2))


f=ju.ffta_single(ff,domain,'x',0)    


emax=np.zeros([nW])
eint=np.zeros([nW])
eabs=np.zeros([nW])
tmax=np.zeros([nW])
tint=np.zeros([nW])
tabs=np.zeros([nW])
smax=np.zeros([nW])
sint=np.zeros([nW])
sabs=np.zeros([nW])
for j in range(nW):
    W=dx*WW[j]
    hf=np.heaviside(x-x2,0.5)
    ff.set_scales(1)
    ff['g'] = ftanh(x1-x,W)
    ff.set_scales(A)
    f2=ff['g']
    ff2=ftanh(x2-x,W)
    e=f2-ff2
    tmax[j]=max(abs(e))
    tint[j]=np.sqrt(sum(e**2)*dx)
    tabs[j]=np.sqrt(sum((ff2-hf)**2)*dx)
    
    ff.set_scales(1)
    ff['g'] = ferf(x1-x,W)
    ff.set_scales(A)
    f2=ff['g']
    ff2=ferf(x2-x,W)
    e=f2-ff2
    emax[j]=max(abs(e))
    eint[j]=np.sqrt(sum(e**2)*dx)
    eabs[j]=np.sqrt(sum((ff2-hf)**2)*dx)

    ff.set_scales(1)
    ff['g'] = fsmooth(x1-x,W)
    ff.set_scales(A)
    f2=ff['g']
    ff2 = fsmooth(x2-x,W)
    e=f2-ff2
    smax[j]=max(abs(e))
    sint[j]=np.sqrt(sum(e**2)*dx)
    sabs[j]=np.sqrt(sum((ff2-hf)**2)*dx)

     
fg, (ax1,ax2,ax3)=plt.subplots(nrows=3,ncols=1)
#ax1.plot(WW,np.log10(tmax),WW,np.log10(emax),WW,np.log10(smax))
#ax2.plot(WW,np.log10(tint),WW,np.log10(eint),WW,np.log10(sint))
ax1.plot(WW,tmax,color='r')
ax1.plot(WW,emax,color='b')
ax1.plot(WW,smax,color='k')
ax2.plot(WW,tint,color='r')
ax2.plot(WW,eint,color='b')
ax2.plot(WW,sint,color='k')
ax3.plot(tabs,tint,color='r', label='tanh')
ax3.plot(eabs,eint,color='b', label='erf')
ax3.plot(sabs,sint,color='k', label='poly')
ax3.legend(loc="upper right")

W=3*dx
ff.set_scales(1)
ff['g'] = ftanh(x1-x,W)
ff.set_scales(A)
f1t=copy.copy(ff['g'])

ff.set_scales(1)
ff['g'] = ferf(x1-x,W)
ff.set_scales(A)
f1e=copy.copy(ff['g'])

ff.set_scales(1)
ff['g'] = fsmooth(x1-x,W)
ff.set_scales(A)
f1s=copy.copy(ff['g'])

f2t = copy.copy(ftanh(  x2-x,W))
f2e = copy.copy(ferf(   x2-x,W))
f2s = copy.copy(fsmooth(x2-x,W))

fg, ((ax11,ax21),(ax12,ax22),(ax13,ax23))=plt.subplots(nrows=3,ncols=2)

ax11.plot((x2-x)/dx,f1t)
ax11.set_xlim((-5,5))
ax11.set_title('tanh')
ax12.plot((x2-x)/dx,f1e)
ax12.set_xlim((-5,5))
ax12.set_title('erf')
ax13.plot((x2-x)/dx,f1s)
ax13.set_xlim((-5,5))
ax13.set_title('poly')

ax21.plot(x2-x,f1t-f2t)
ax21.set_xlim((-0.2,0.2))
ax21.set_title('tanh')
ax22.plot(x2-x,f1e-f2e)
ax22.set_title('erf')
ax22.set_xlim((-0.2,0.2))
ax23.plot(x2-x,f1s-f2s)
ax23.set_xlim((-0.2,0.2))
ax23.set_title('poly')

#ax1.plot(x2,f2,x2,ff2)

#ax2.plot(x1,np.maximum(0,np.maximum(f1-1,-f1)))
#ax2.plot(x2,np.maximum(0,np.maximum(f2-1,-f2)))
#ax2.plot(x1,np.maximum(0,np.maximum(f1-1,-f1)))

#ax2.plot(x2-x,e)
#ax2.set_xlim((-0.5,0.5))
#ax1.set_title('max(e) {:8.5f}, int(e) {:8.5f}, W/dx={:4.1f}'.format(max(abs(e)),np.sqrt(sum(e**2)*dx),W/dx))

plt.draw()
plt.waitforbuttonpress()
