#!/usr/bin/env python3
import numpy as np
import scipy.special as sp
import scipy.interpolate as it
import matplotlib.pyplot as plt

def front1(e,x):
    btol=np.max(e)/20
    N=len(x)
    if e[-1]>btol: return(x[-1])
    for i in range(N-1,1,-1):
        b1=e[i-1]- btol
        b2=e[i]  - btol
        if e[i-1]>btol and e[i]<btol and e[i]>0: return((b1*x[i]-b2*x[i-1])/(b1-b2))
    return(np.nan)

# Look for an inflection point in the cubic and a strictly decreasing region
# This is only going to work well for a fully resolved simulations
def front2(e,x):
    N=len(x)
    e2=0
    e3=0
    e4=0
    for i in range(N-1,1,-1):
        e1=max(0,e[i])
        if i<N-1: e2=max(0,e[i+1])
        if i<N-2: e3=max(0,e[i+2])
        if i<N-3: e4=max(0,e[i+3])
        x1=x[min(N-1,i+1)]
        x2=x[min(N-1,i+2)]
        dd1 = e1-2*e2+e3
        dd2 = e2-2*e3+e4
        if e1>2*e2 and e2>0   and e3==0 and e4==0:   return((e2*x2+(e1-2*e2)*x1)/(e1-e2))
        if e1>=e2  and e2>=e3 and e3>e4 and e4>0 and dd1<=0 and dd2>=0 and dd1<dd2: return((dd1*x2-dd2*x1)/(dd1-dd2))
    return(np.nan)

def front3(f,x):
    e=it.pchip_interpolate(x, np.maximum(0,f), x, der=2)
    N=len(x)
    if e[-1]<0: return(x[-1])
    for i in range(N-1,1,-1):
#        if i<N-2 and i>2:
#            if f[i]<=0 and f[i-1]<=0 and f[i-2]>=0 and f[i-3]>0:
#                f1=f[i-3]
#                f2=f[i-2]
#                x1=f[i-3]
#                x2=f[i-2]
#                return((f2*x2+(f1-2*f2)*x1)/(f1-f2))
        
        e1=e[i-1]
        e2=e[i] 
        if e1<0 and e2>=0: return((e1*x[i]-e2*x[i-1])/(e1-e2))
    return(np.nan)


# Actually we want the point of maximum curvature
# To make this smooth we need higher interpolation
# Or just require that points are zero ??
# we want to find the minimum of f'' if we have clipped points easy otherwise solve a quadratic
# f''=0 is the same a finding the minimum of f' which is probably more stable

plt.ion()

Nx=501
Ny=10

x=np.linspace(-3,3, num=Nx, endpoint=True, dtype=np.float64)

y=np.linspace(-1,2, num=Ny, endpoint=True, dtype=np.float64)

f,ax=plt.subplots(nrows=4,ncols=1)
X=np.zeros_like(y)

for j in range(Ny):
    f=sp.erf(20*(y[j]-x))+sp.erf(10*(2+x))
    #f=np.minimum(1,np.maximum(0,y[j]-x))
    X[j]=front3(f,x)
    ax[0].plot(x,f)
    
ax[1].plot(y,X)
ax[2].plot(y[0:-2],X[1:-1]-X[0:-2])
ax[3].plot(y,y-X)
   

plt.draw()
plt.waitforbuttonpress()
