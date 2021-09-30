#!/usr/bin/env python3
import jutil as ju
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.special as sp

plt.ion()
# Test integration functin sin jutil

N=1e4
dx=1e-3
x=np.linspace(-1, 1, num=N, endpoint=True, dtype=np.float64)    
aa=np.logspace(-3, 3, num=20, endpoint=True, dtype=np.float64)
aa=np.append([0],aa)

#aa=[1,2]
#na=len(aa)
#fg, axa=plt.subplots(nrows=na,ncols=1)
#fg, axb=plt.subplots(nrows=na,ncols=1)
#n=0
#f=ju.pgfn(x,n)
#for j in range(na):
#    a=aa[j]
#    axa[j].plot(x,ju.ndif2(ju.pgf1n,x,dx,n,a)-f*np.exp(-(a*x)**2))
#    axb[j].plot(x,np.sqrt(math.pi)*sp.erf(a*x)/(2*a)-ju.pgf1n(x,n,a))
    
a=1
f20=(-2*x*np.exp(-(a*x)**2)*a+np.sqrt(math.pi)*sp.erf(a*x))/(4*a**3)

fg, ax=plt.subplots(nrows=6,ncols=1)
for n in range(6): ax[n].plot(x,ju.pgfn(x,n))
fg, ax=plt.subplots(nrows=6,ncols=1)
for n in range(6): ax[n].plot(x,ju.pgf1n(x,n,a))
fg, ax=plt.subplots(nrows=6,ncols=1)
for n in range(6): ax[n].plot(x,ju.pgf2nb(x,n,(a*x)**2))

nn=range(6)
for n in nn:
    print("f1     n     {}".format(n))
    for j in range(len(aa)):
        a=aa[j]
        f=ju.pgfn(x,n)
        g1=f*np.exp(-(a*x)**2)
        g2=ju.ndif2(ju.pgf1n,x,dx,n,a)
        maxag=np.maximum(1,np.abs(g1))
        e=max(np.abs(g1-g2)/maxag)
        print('dx={:7.1e}, a={:7.1e} e={:10.4e}'.format(dx,a,e))
    print(' ')
nn=range(6)
for n in nn:
    print("f2     n     {}".format(n))
    for j in range(len(aa)):
        a=aa[j]
        f=ju.pgfn(x,n)
        g1=x**2*f*np.exp(-(a*x)**2)
        g2=ju.ndif2(ju.pgf2n,x,dx,n,a)
        maxag=np.maximum(1,np.abs(g1))
        e=max(np.abs(g1-g2)/maxag)
        print('dx={:7.1e}, a={:7.1e} e={:10.4e}'.format(dx,a,e))
    print(' ')

plt.draw()
plt.waitforbuttonpress()

quit()

N=1e4

x=np.linspace(-10, 10, num=N, endpoint=True, dtype=np.float64)    
dx=2e-2
print('Testing intxncos')
for n in range(21):
    g1=x**n*np.cos(x)
    g2=ju.ndif1(ju.intxncos,x,dx,n)
    e=max(np.abs(g1-g2))
    print('n={:2.0f}, abs error={:10.4e}, rel error={:10.4e}'.format(n,e,e/max(np.abs(g1))))

print('Testing intexpcosabx')
aa=[0.001,0.01,0.1,1,10,100]
aa=np.logspace(-3, 3, num=20, endpoint=True, dtype=np.float64)
bb=np.logspace(-1, 1, num=20, endpoint=True, dtype=np.float64)  


for n in range(len(aa)):
    for j in range(len(bb)):
        a=aa[n]
        b=bb[j]
        dx=1e-3/np.sqrt(a)
        g=np.exp(-(a*x)**2)*np.cos(b*x)
        #g1=ju.ndif2(ju.intexpcosabx1,x,dx,a,b)
        #g2=ju.ndif2(ju.intexpcosabx2,x,dx,a,b)
        g3=ju.ndif2(ju.intexpcosabx,x,dx,a,b)
        maxag=np.maximum(1,np.abs(g))
        #e1=max(np.abs(g1-g)/maxag)
        #e2=max(np.abs(g2-g)/maxag)
        e3=max(np.abs(g3-g)/maxag)
        #print('dx={:7.1e}, a={:7.1e}, e=[{:10.4e} {:10.4e} {:10.4e}]'.format(dx,a,e1,e2,e3))
        print('dx={:7.1e}, a={:7.1e}, b={:7.1e}, e={:10.4e}'.format(dx,a,b,e3))

    




