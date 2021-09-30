#!/usr/bin/env python3
import numpy as np
import scipy.special as sp
import math

N=25
pi=math.pi
sqrtpi=np.sqrt(math.pi)

atol=0.5
sminreal=1e-100 # sqrt of minimum real

def F0(r,a): # diff(F0(r,a),r)
    a=np.maximum(sminreal,np.abs(a))
    if np.isscalar(a):
        if a>atol: f=F0a(r,a)
        else:      f=F0b(r,a)
    elif np.isscalar(r):
        f=np.zeros_like(a)
        f[a >  atol]=F0a(r,a[a >  atol])
        f[a <= atol]=F0b(r,a[a <= atol])
    else:
        r,a=np.broadcast_arrays(r, a)
        f=np.zeros_like(r)
        f[a >  atol]=F0a(r[a >  atol],a[a >  atol])
        f[a <= atol]=F0b(r[a <= atol],a[a <= atol])
    return(f)

def F1(r,a): # diff(F0(r,a),r)
    a=np.maximum(sminreal,np.abs(a))
    if np.isscalar(a):
        if a>atol: f=F1a(r,a)
        else:      f=F1b(r,a)
    elif np.isscalar(r):
        f=np.zeros_like(a)
        f[a >  atol]=F1a(r,a[a >  atol])
        f[a <= atol]=F1b(r,a[a <= atol])
    else:
        r,a=np.broadcast_arrays(r, a)
        f=np.zeros_like(r*a)
        f[a >  atol]=F1a(r[a >  atol],a[a >  atol])
        f[a <= atol]=F1b(r[a <= atol],a[a <= atol])
    return(f)

def F2(r,a):  # diff(F0(r,a),a)
    s=np.sign(a)
    a=np.maximum(sminreal,np.abs(a))
    if np.isscalar(a):
        if a>atol: f=s*F2a(r,a)
        else:      f=s*F2b(r,a)
    elif np.isscalar(r):
        f=np.zeros_like(a)
        f[a >  atol]=s[a >  atol]*F2a(r,a[a >  atol])
        f[a <= atol]=s[a <= atol]*F2b(r,a[a <= atol])
    else:
        r,a=np.broadcast_arrays(r, a)
        f=np.zeros_like(r)
        f[a >  atol]=s[a >  atol]*F2a(r[a >  atol],a[a >  atol])
        f[a <= atol]=s[a <= atol]*F2b(r[a <= atol],a[a <= atol])
    return(f)

def F3(r,a):  # diff(F0(r,a),a)
    s=np.sign(a)
    a=np.maximum(sminreal,np.abs(a))
    if np.isscalar(a):
        if a>atol: f=s*F3a(r,a)
        else:      f=s*F3b(r,a)
    elif np.isscalar(r):
        f=np.zeros_like(a)
        f[a >  atol]=s[a >  atol]*F3a(r,a[a >  atol])
        f[a <= atol]=s[a <= atol]*F3b(r,a[a <= atol])
    else:
        r,a=np.broadcast_arrays(r, a)
        f=np.zeros_like(r*a)
        f[a >  atol]=s[a >  atol]*F3a(r[a >  atol],a[a >  atol])
        f[a <= atol]=s[a <= atol]*F3b(r[a <= atol],a[a <= atol])
    return(f)

def F0a(r,a):
    f=sp.erf(a*r)
    for j in range(1,N): f=f-sp.erfc(a*(r+j))+sp.erfc(a*(j-r)) # derf(a*r,j*a)
    return(f/2)

def F1a(r,a):
    f=a*np.exp(-(r*a)**2)/sqrtpi
    for n in range(1,N): f=f+a/sqrtpi*(np.exp(-a**2*(r-n)**2)+np.exp(-a**2*(r+n)**2))
    return(f)

def F2a(r,a):
    f=r/a*F1(r,a)
    for n in range(1,N): f=f+n/sqrtpi*(np.exp(-a**2*(r+n)**2)-np.exp(-a**2*(r-n)**2))
    return(f)

def F3a(r,a):
    f=np.zeros_like(r*a)
    for n in range(-N,N+1): f=f+1/sqrtpi*(1-2*a**2*(n+r)**2)*np.exp(-a**2*(n+r)**2)
    return(f)

def F0b(r,a):
    f=r+np.exp(-pi/a**2)*np.sin(2*pi*r)/(2*pi)
    return(f)

def F1b(r,a):
    f=1+np.exp(-pi/a**2)*np.cos(2*pi*r)
    return(f)

def F2b(r,a):
    f=np.exp(-pi/a**2)*np.sin(2*pi*r)/a**3
    return(f)

def F3b(r,a):
    f=2*pi*np.exp(-pi/a**2)*np.cos(2*pi*r)/a**3
    return(f)

def F1c(r,a):
    n=3 # JacobiTheta3
    z=1j*a*r
    q=np.exp(-a**2)
    f=np.zeros_like(r*a)
    for j in range(len(z)): f[j]=a/sqrtpi*np.exp(-(a*r[j])**2)*np.real(mp.jtheta(n, z[j], q, derivative=0))
    return(f)

def periodic_gaussian(a,da,y,z,U,W,H,r):
    a,da,y,z=np.broadcast_arrays(a,da,y,z)
    Fy=F1(y/W,a*W)
    Fz=F1(z/H,a*H)
    u = r**2*pi*U/(H*W)*Fz*Fy
    v =-r**2*pi*U*W/H*da*Fz*F2(y/W,a*W)
    w =-r**2*pi*U*H/W*da*Fy*F2(z/H,a*H)
    return(u,v,w)

#unnormalised
def periodic_gaussianu(a,y,z,W,H):
    a,y,z=np.broadcast_arrays(a,y,z)
    Fy=F1(y/W,a*W)
    Fz=F1(z/H,a*H)
    u = Fz*Fy
    return(u)



