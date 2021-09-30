#!/usr/bin/env python3
import numpy as np
import math
#import mpmath as mp
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
#import warnings
import serf as se
import random

plt.ion()

# 5th order accurate finite difference check
def check_diff(x,d,f,df,s):
    adf=(-f(x+2*d)+f(x-2*d)+8*f(x+d)-8*f(x-d))/(12*d)
    e=df(x)-adf
    emax=np.max(np.abs(e))
    erms=np.sqrt(np.mean(e**2))
    print('{} rms: {:8.1e}, max: {:8.1e}'.format(s,emax,erms))
    return(adf)

def res(z):
    return(np.sqrt((W*(y-z)**2).sum()/H))

H=random.uniform(0, 10)
h=random.uniform(-H, 2*H)
w=0.1*H*random.uniform(0,1)**2
U=random.gauss(0,1)**2
B=random.gauss(0,1)**2
N=random.randrange(10,1000)

x=np.linspace(0, H, num=N, endpoint=True, dtype=np.float)    
dx=x[1]-x[0]
W=dx*np.ones((N,), dtype=np.float)
W[0]=dx/2
W[-1]=dx/2

u=se.fu(x,W,H,A=U,h=h,w=w)
b=se.fb(x,W,H,A=B,h=h,w=w)

bAhw = se.fb(x,W,H,None,None,None)
bhw  = se.fb(x,W,H,B,None,None)
bAw  = se.fb(x,W,H,None,h,None)
bAh  = se.fb(x,W,H,None,None,w)
bA   = se.fb(x,W,H,None,h,w)
bh   = se.fb(x,W,H,B,None,w)
bw   = se.fb(x,W,H,B,h,None)

uAhw = se.fu(x,W,H,None,None,None)
uhw  = se.fu(x,W,H,U,None,None)
uAw  = se.fu(x,W,H,None,h,None)
uAh  = se.fu(x,W,H,None,None,w)
uA   = se.fu(x,W,H,None,h,w)
uh   = se.fu(x,W,H,U,None,w)
uw   = se.fu(x,W,H,U,h,None)


tol=5e-5
check_diff(x,tol,u.f,u.df,'u du')
check_diff(x,tol,u.If,u.f,'Iu u')

check_diff(x,tol,b.f,b.df,'b db')
check_diff(x,tol,b.If,b.f,'Ib b')

fg, ((ax1,ax4),(ax2,ax5),(ax3,ax6))=plt.subplots(nrows=3,ncols=2)
ax1.plot(x,u.df(x))
ax2.plot(x,u.f(x))
ax3.plot(x,u.If(x))
ax4.plot(x,b.df(x))
ax5.plot(x,b.f(x))
ax6.plot(x,b.If(x))
ax1.set_ylabel('du')
ax2.set_ylabel('u')
ax3.set_ylabel('Iu')
ax4.set_ylabel('db')
ax5.set_ylabel('b')
ax6.set_ylabel('Ib')

y=b.f(x)+B*np.random.normal(0, .1, x.shape)

fg, ((ax1,ax5),(ax2,ax6),(ax3,ax7),(ax4,ax8))=plt.subplots(nrows=4,ncols=2)


bAhw.fit(y)
bAh.fit(y)
bAw.fit(y)
bhw.fit(y)
bA.fit(y)
bh.fit(y)
bw.fit(y)


if True:
    ax1.plot(x,y); z=bAhw.f(x);ax1.plot(x,z) ; ax1.set_ylabel('Ahw {:8.4f}'.format(res(z)))
    ax2.plot(x,y); z=bAh.f(x) ;ax2.plot(x,z) ; ax2.set_ylabel('Ah  {:8.4f}'.format(res(z)))
    ax3.plot(x,y); z=bAw.f(x) ;ax3.plot(x,z) ; ax3.set_ylabel('Aw  {:8.4f}'.format(res(z)))
    ax4.plot(x,y); z=bhw.f(x) ;ax4.plot(x,z) ; ax4.set_ylabel('hw  {:8.4f}'.format(res(z)))
    ax5.plot(x,y); z=bA.f(x)  ;ax5.plot(x,z) ; ax5.set_ylabel('A   {:8.4f}'.format(res(z)))
    ax6.plot(x,y); z=bh.f(x)  ;ax6.plot(x,z) ; ax6.set_ylabel('h   {:8.4f}'.format(res(z)))
    ax7.plot(x,y); z=bw.f(x)  ;ax7.plot(x,z) ; ax7.set_ylabel('w   {:8.4f}'.format(res(z)))
    
    
fg, ((ax1,ax5),(ax2,ax6),(ax3,ax7),(ax4,ax8))=plt.subplots(nrows=4,ncols=2)


y=u.f(x)+U*np.random.normal(0, .1, x.shape)
uAhw.fit(y)
uAh.fit(y)
uAw.fit(y)
uhw.fit(y)
uA.fit(y)
uh.fit(y)
uw.fit(y)

mu=(u.f(0)+u.f(H))/2
mb=(b.f(0)+b.f(H))/2
xu=u.solve(mu)
xb=u.solve(mb)



ax1.plot(x,y); z=uAhw.f(x);ax1.plot(x,z) ; ax1.set_ylabel('Ahw {:8.4f}'.format(res(z)))
ax2.plot(x,y); z=uAh.f(x) ;ax2.plot(x,z) ; ax2.set_ylabel('Ah  {:8.4f}'.format(res(z)))
ax3.plot(x,y); z=uAw.f(x) ;ax3.plot(x,z) ; ax3.set_ylabel('Aw  {:8.4f}'.format(res(z)))
ax4.plot(x,y); z=uhw.f(x) ;ax4.plot(x,z) ; ax4.set_ylabel('hw  {:8.4f}'.format(res(z)))
ax5.plot(x,y); z=uA.f(x)  ;ax5.plot(x,z) ; ax5.set_ylabel('A   {:8.4f}'.format(res(z)))
ax6.plot(x,y); z=uh.f(x)  ;ax6.plot(x,z) ; ax6.set_ylabel('h   {:8.4f}'.format(res(z)))
ax7.plot(x,y); z=uw.f(x)  ;ax7.plot(x,z) ; ax7.set_ylabel('w   {:8.4f}'.format(res(z)))


plt.draw()
plt.waitforbuttonpress()
