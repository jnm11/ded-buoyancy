#!/bin/python3
import numpy as np
import math
import scipy.special as sp
import scipy.optimize as opt
import matplotlib.pyplot as plt

plt.ion()

def ffta_heaviside(x,k,M):
    xx=x;
    xc = (2*M+1)/k/2
    x  = np.minimum(math.pi/2,np.maximum(-math.pi/2,k*x/(2*M+1)))
    If = np.maximum(0,(x+math.pi)*xc)+np.maximum(0,xx-xc*math.pi)
    f  = 1/2+np.zeros_like(x)
    If0=0
    df = np.zeros_like(x)
    for n in range(M+1):
        sf = k*(2*n+1)/(2*M+1)
        c=2/math.pi*np.exp(2*sp.gammaln(M+3/2)-sp.gammaln(M+2+n)-sp.gammaln(M+1-n))/(2*n+1)
        If  += -c/sf*np.cos((2*n+1)*x)
        f   +=  c   *np.sin((2*n+1)*x)
        df  +=  c*sf*np.cos((2*n+1)*x)
    df0 = k*np.exp(sp.gammaln(2*M+1)-M*np.log(4)-2*sp.gammaln(M+1))
    return(f,df,df0,If)



#w*(u-A*f)^2
#r=wuu-2Awuf+A^2wff
#A=wfu/wff
#r=wuu-2(wfu)^2/wff

def pf(x,X,p):
    tol=1e-5
    y=np.maximum(x-X,tol)
    f=np.power(y,p)
    df=np.where(x<=X,0,p*np.power(y,p-1))
    return(f,df)

def qf(u,x,w,X,p):
    f=pf(x,X,p)[0]
    wfu=(f*u*w).sum()
    wff=(f*f*w).sum()
    A=wfu/wff
    r=(u*u*w).sum()-A*wfu
    return(r,A)

def dqf(u,x,w,X,p):
    f,df=pf(x,X,p)
    wfu=(f*u*w).sum()
    wff=(f*f*w).sum()
    dwfu=-(df*u*w).sum()
    dwff=-2*(df*f*w).sum()
    dr=-2*wfu*dwfu/wff+wfu**2*dwff/wff**2
    return(dr)

def fitAX(u,x,w,p,Xmin,Xmax):
    f  = lambda X: qf(u,x,w,X,p)[0]
    df = lambda X: dqf(u,x,w,X,p)
    X=opt.fmin_bfgs(f,[(Xmin+Xmax)/2],fprime=df)
    #X = opt.brent(f,(),(Xmin,Xmax))
    A=qf(u,x,w,X,p)[1]
    return(A,X)

M=10
L=30
N=1024
x=np.linspace(0,L, num=N, endpoint=True, dtype=np.float)
dx=x[1]-x[0]
maxkx=math.pi/dx


X0 = 2*M*math.pi/maxkx # Top of fluid source region
#weight function 
w=ffta_heaviside(x-X0,maxkx,M)[0]-ffta_heaviside(x-L+X0,maxkx,M)[0]
A=np.random.normal()
X=-np.random.uniform()
p=-1/5
print('A {:6.3f}'.format(A))
print('X {:6.3f}'.format(X))
print('p {:6.3f}'.format(p))

u=A*pf(x,X,p)[0]+A*0.1*np.random.normal(size=N)
u=u/np.sqrt((u*u*w).sum())
Xmin=X-1
Xmax=X+1
B,Y=fitAX(u,x,w,p,Xmin,Xmax)
r1,B1=qf(u,x,w,Xmin,p)
r2,B2=qf(u,x,w,X,p)
r3,B3=qf(u,x,w,Xmax,p)
r4,B4=qf(u,x,w,Y,p)
v=B*pf(x,Y,p)[0]

print('X {:6.3f} B {:6.3f} r {:6.3f}'.format(Xmin,B1,r1))
print('X {:6.3f} B {:6.3f} r {:6.3f}'.format(X,   B2,r2))
print('X {:6.3f} B {:6.3f} r {:6.3f}'.format(Xmax,B3,r3))
print('Y {}  B {} r {}'.format(Y,   B ,r4))
print('Y {:6.3f} B {:6.3f} r {:6.3f}'.format(Y,   B ,r4))


fg, a0=plt.subplots(nrows=3,ncols=1)
a0[0].plot(x,w*u,x,w*v)
a0[1].plot(x,w*(u-v))
a0[2].plot(x,w)


plt.draw()
plt.waitforbuttonpress()
quit(1)
