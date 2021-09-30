#!/usr/bin/env python3
import numpy as np
import scipy.special as sp
import math
import mpmath as mp
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import warnings
import jutil as ju
import jtheta

warnings.simplefilter('ignore',np.ComplexWarning)

plt.ion()


U=1
W=2
H=3
Nx=99;
Ny=100;
Nz=101;

xx=np.linspace(0,1, num=Nx, endpoint=True, dtype=np.float64)
yy=np.linspace(-W/2,W/2, num=Ny, endpoint=True, dtype=np.float64)
zz=np.linspace(-H/2,H/2, num=Nz, endpoint=True, dtype=np.float64)
dx=xx[1]-xx[0]
dy=yy[1]-yy[0]
dz=zz[1]-zz[0]
x=xx.reshape((Nx,1,1))
y=yy.reshape((1,Ny,1))
z=zz.reshape((1,1,Nz))
x,y,z=np.broadcast_arrays(x,y,z)
a=16*x**2*(1-x)**2
da=32*x*(1-x)*(1-2*x)
print('x.shape {}'.format(x.shape))
print('y.shape {}'.format(y.shape))
print('z.shape {}'.format(z.shape))

u,v,w=jtheta.periodic_gaussian(a,da,y,z,U,W,H,1)

TU=u.sum((1,2))*dy*dz
TV=v.sum((1,2))*dy*dz
TW=w.sum((1,2))*dy*dz
f,ax0=plt.subplots(nrows=3,ncols=1)
f,ax1=plt.subplots(nrows=2,ncols=2)
ax0[0].plot(xx,TU)
ax0[1].plot(xx,TV)
ax0[2].plot(xx,TW)
print('xx.shape {}'.format(xx.shape))
print('yy.shape {}'.format(yy.shape))
print('zz.shape {}'.format(zz.shape))
print('u.shape {}'.format(u.shape))
print('v.shape {}'.format(v.shape))
print('w.shape {}'.format(w.shape))
ax1[0,0].plot(xx,v[:, 0,:].reshape((Nx,Nz)))
ax1[0,0].plot(xx,v[:,-1,:].reshape((Nx,Nz)))
ax1[0,0].set_xlabel('x')
ax1[0,0].set_ylabel('v')
ax1[0,1].plot(zz,v[:, 0,:].reshape((Nx,Nz)).transpose((1,0)))
ax1[0,1].plot(zz,v[:,-1,:].reshape((Nx,Nz)).transpose((1,0)))
ax1[0,1].set_xlabel('z')
ax1[0,1].set_ylabel('v')
ax1[1,0].plot(xx,w[:,:, 0].reshape((Nx,Ny)))
ax1[1,0].plot(xx,w[:,:,-1].reshape((Nx,Ny)))
ax1[1,0].set_xlabel('x')
ax1[1,0].set_ylabel('w')
ax1[1,1].plot(yy,w[:,:, 0].reshape((Nx,Ny)).transpose((1,0)))
ax1[1,1].plot(yy,w[:,:,-1].reshape((Nx,Ny)).transpose((1,0)))
ax1[1,1].set_xlabel('y')
ax1[1,1].set_ylabel('w')



plt.draw()
plt.waitforbuttonpress()



quit()

def F1c(r,a):
    n=3 # JacobiTheta3
    z=1j*a*r
    q=np.exp(-a**2)
    f=np.zeros_like(r)
    for j in range(len(z)): f[j]=a/sqrtpi*np.exp(-(a*r[j])**2)*np.real(mp.jtheta(n, z[j], q, derivative=0))
    return(f)

def F0ar(a,r): # swapped arguments for calculating derivatives
    return(jtheta.F0(r,a))

N=5000

L=1

r=np.linspace(-L*1.1, L*1.1, num=N, endpoint=True, dtype=np.float64)    
aa=np.array([0,1e-3,jtheta.atol*.99,jtheta.atol/.99,1,10,100])
na=len(aa)
f,ax0=plt.subplots(nrows=na,ncols=1,sharex=True)
f,ax1=plt.subplots(nrows=na,ncols=3,sharex=True)
f,ax2=plt.subplots(nrows=na,ncols=3,sharex=True)
f,ax3=plt.subplots(nrows=na,ncols=3,sharex=True)

dr=1e-4
for j in range(na):
    a=aa[j]
    
    if j==0:
        ax0[j].set_title('jtheta.F0')
        ax1[j,0].set_title('jtheta.F1')
        ax2[j,0].set_title('jtheta.F2')
        ax3[j,0].set_title('jtheta.F3')
        
    ax0[j].plot(r,jtheta.F0(r,a)-r)
    ax0[j].set_ylabel('{}'.format(a))
    ax1[j,0].plot(r,jtheta.F1(r,a))
    ax1[j,1].plot(r,ju.ndif1(jtheta.F0,r,dr,a))
    ax1[j,2].plot(r,jtheta.F1(r,a)-ju.ndif1(jtheta.F0,r,dr,a))
    ax1[j,0].set_ylabel('{}'.format(a))
    
    
    da=max(1e-5,a*1e-4)
    ax2[j,0].plot(r,jtheta.F2(r,a))
    ax2[j,1].plot(r,ju.ndif1(F0ar,a,da,r))
    ax2[j,2].plot(r,jtheta.F2(r,a)-ju.ndif1(F0ar,a,da,r))
    ax2[j,0].set_ylabel('{}'.format(a))
    
    ax3[j,0].plot(r,jtheta.F3(r,a))
    ax3[j,1].plot(r,ju.ndif1(jtheta.F2,r,dr,a))
    ax3[j,2].plot(r,jtheta.F3(r,a)-ju.ndif1(jtheta.F2,r,dr,a))
    ax3[j,0].set_ylabel('{}'.format(a))
    


rr=np.array([0,0.01,0.1,0.5])
a=np.linspace(-5,5, num=N, endpoint=True, dtype=np.float64)
nr=len(rr)
f,ax0=plt.subplots(nrows=nr,ncols=1,sharex=True)
f,ax1=plt.subplots(nrows=nr,ncols=1,sharex=True)
f,ax2=plt.subplots(nrows=nr,ncols=1,sharex=True)
f,ax3=plt.subplots(nrows=nr,ncols=1,sharex=True)

for j in range(nr):
    r=rr[j]

    if j==0:
        ax0[j].set_title('F0')
        ax1[j].set_title('F1')
        ax2[j].set_title('F2')
        ax3[j].set_title('F3')
        
    ax0[j].plot(a,jtheta.F0(r,a))
    ax1[j].plot(a,jtheta.F1(r,a))
    ax2[j].plot(a,jtheta.F2(r,a))
    ax3[j].plot(a,jtheta.F3(r,a))

 
plt.draw()
plt.waitforbuttonpress()


