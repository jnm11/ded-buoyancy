#!/usr/bin/env python3
# test_alias
"""
Test spectral aliasing of front 

Usage: 
   test_alias.py [options]

Options:
  -L, --Length=<>  : Length of domain       [default:  1]
  -N, --Number=<>  : Number of grid points  [default:  64]
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
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys
import jutil as ju
import logging
import math
import copy
import scipy.special as sp
import ffta
import dftf
#import QtWidgets from Qt
#import QtCore from Qt
#QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
plt.ion()

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

ju.logger=logging.getLogger(__name__)



# 5th order accurate finite difference check
def check_diff(x,d,f,df,s):
    e=df(x)-(-f(x+2*d)+f(x-2*d)+8*f(x+d)-8*f(x-d))/(12*d)
    emax=np.max(np.abs(e))
    erms=np.sqrt(np.mean(e**2))
    print('{} rms: {:8.1e}, max: {:8.1e}'.format(s,emax,erms))
        

L=np.float64(args['--Length'])
N=np.int64(  args['--Number'])
W=np.float64(args['--Width'])
M=np.int64(  args['--Order'])
x=np.float64(args['--Front'])
A=np.float64(args['--Alias'])
y1=np.float64(args['--y1'])
y2=np.float64(args['--y2'])

NN=np.int64(A*N)
nfx=10

d=list()
ff=list()
fs=list()
w=list()

LL=(-0.1,1.2)
mL=(LL[0]+LL[1])/2;



N=100+np.asarray([0,1,2,3])

for M in range(5):
    dfg, a=plt.subplots(nrows=4,ncols=3)
    if M==0: f = lambda g,d: dftf.sawtooth(g,d,'x')
    else:    f = lambda g,d: dftf.heaviside(g,d,'x',M)
    for j in range(len(N)):
        c=(de.SinCos(   'x', N[j], interval = LL, dealias=A),)
        d=de.Domain(c, np.float64)
        x = d.bases[0].grid(scale=1)
        g = d.new_field()
        if M==0: sf=1
        else:    sf=0
        g.meta['x']['parity']  = -1
        f(g,d)
        a[j,0].plot(x,g['g'],x,sf*(x-LL[0]))
        
        g.meta['x']['parity']  = 1
        f(g,d)
        a[j,1].plot(x,g['g'],x,sf*(x-mL))
        
        c=(de.Fourier('x', N[j], interval = LL, dealias=A),)
        d=de.Domain(c, np.float64)
        x = d.bases[0].grid(scale=1)
        g = d.new_field()
        f(g,d)
        a[j,2].plot(x,g['g'],x,sf*(x-mL))
    
plt.draw()

cFR=de.Fourier(   'x', N, interval = LL, dealias=A)
cSC=de.SinCos(    'x', N, interval = LL, dealias=A)
cCB=de.Chebyshev( 'x', N, interval = LL, dealias=A)
c=(cFR,cSC,cCB)

nc   = len(c)
MM=[5,7,9,11,13,15,17]
x=np.linspace(-L,L,1e4)
dx=x[1]-x[0]
tol=1e-4;
nM=len(MM)

d=de.Domain([c[1]], np.float64)
b=d.bases[0]
x1 = b.grid(scale=1)
x2 = b.grid(scale=A)
FFTA=ffta.ffta(b,5,nfx)

x1=0.1
x2=0.9
x=0.51
for n in range(1,10):
    Hn1=FFTA.Hn(x,x1,x2,n)
    Hn2=FFTA.Hn2(x,x1,x2,n)
    ju.logger.info('n {}, Hn1 {}, Hn2 {}, {}'.format(n,Hn1,Hn2,Hn1-Hn2))
quit(1)

ni=4
nd=3

g   = d.new_field()
gs  = d.new_field()
gc  = d.new_field()


dfg, da0=plt.subplots(nrows=ni,ncols=nd)

g   = d.new_field()
gs  = d.new_field()
gc  = d.new_field()
for d in range(nd):
    if np.mod(d,2)==0:
        g.meta['x']['parity']  = -1
        gs.meta['x']['parity'] =  1
        gc.meta['x']['parity'] = -1
    else:
        g.meta['x']['parity']  =  1
        gs.meta['x']['parity'] = -1
        gc.meta['x']['parity'] =  1
    for i in range(ni):
        g.set_scales(1)
        gs.set_scales(1)
        gc.set_scales (1)
        fd   = lambda x: FFTA.delta(x-mL,i+1)
        if   i==0:
            f   = lambda x: FFTA.odd_constant(x,1,0)
            df  = lambda x: FFTA.odd_constant(x,1,1)
            ddf = lambda x: FFTA.odd_constant(x,1,2)
        elif i==1: 
            f   = lambda x: FFTA.odd_constantd3(x,1,0)
            df  = lambda x: FFTA.odd_constantd3(x,1,1)
            ddf = lambda x: FFTA.odd_constantd3(x,1,2)
        elif i==2:
            f   = lambda x: FFTA.periodic_constant(x,1,0)
            df  = lambda x: FFTA.periodic_constant(x,1,1)
            ddf = lambda x: FFTA.periodic_constant(x,1,2)
            g.meta['x']['parity']  = 1
            gs.meta['x']['parity'] = 1
            gc.meta['x']['parity'] = 1
        elif i==3:
            f   = lambda x: FFTA.even_constant(x,1,0)
            df  = lambda x: FFTA.even_constant(x,1,1)
            ddf = lambda x: FFTA.even_constant(x,1,2)
            g.meta['x']['parity']  = 1
            gs.meta['x']['parity'] = 1
            gc.meta['x']['parity'] = 1


        fg, a0=plt.subplots(nrows=6)
        if d==0:
            g['g']   = f(x1)
            g2       = f(x2)
            gs['g']  = f(x1)**2
            gc['g']  = f(x1)**3
            gs2      = f(x2)**2
            gc2      = f(x2)**3
        elif d==1:
            g['g']   = df(x1)
            g2       = df(x2)
            gs['g']  = 2*f(x1)*df(x1)
            gc['g']  = 3*f(x1)**2*df(x1)
            gs2      = 2*f(x2)*df(x2)
            gc2      = 3*f(x2)**2*df(x2)
        elif d==2:
            g['g']   = ddf(x1)
            g2       = ddf(x2)
            gs['g']  = 2*df(x1)**2       + 2*f(x1)   *ddf(x1)
            gc['g']  = 6*f(x1)*df(x1)**2 + 3*f(x1)**2*ddf(x1)
            gs2      = 2*df(x2)**2       + 2*f(x2)   *ddf(x2)
            gc2      = 6*f(x2)*df(x2)**2 + 3*f(x2)**2*ddf(x2)
            
        if   i==0:ts='ffta.odd_constant      i={}, d={}'.format(i,d)
        elif i==1:ts='ffta.odd_constantd3    i={}, d={}'.format(i,d)
        elif i==2:ts='ffta.periodic_constant i={}, d={}'.format(i,d)
        elif i==3:ts='ffta.even_constant     i={}, d={}'.format(i,d)
        fg.suptitle(ts)
        check_diff(x1,tol,f,df,ts)
        check_diff(x1,tol,df,ddf,ts)
        
       
        a0[0].set_title('1')
        a0[1].set_title('2')
        a0[2].set_title('3')
        a0[0].plot(x1,g['g'],  marker='s')
        a0[1].plot(x1,gs['g'], marker='s')
        a0[2].plot(x1,gc['g'], marker='s')
        g.set_scales(A)
        gs.set_scales(A)
        gc.set_scales(A)
        a0[0].plot(x2,g['g'])
        a0[1].plot(x2,gs['g'])
        a0[2].plot(x2,gc['g'])
        a0[0].plot(x2,g2)
        a0[1].plot(x2,gs2)
        a0[2].plot(x2,gc2)
        a0[3].plot(x2,g['g']-g2)
        a0[4].plot(x2,gs['g']-gs2)
        a0[5].plot(x2,gc['g']-gc2)

        da0[i,d].plot(x1,fd(x1)**(1+d), marker='s')
        da0[i,d].plot(x2,fd(x2)**(1+d))
        S1=fd(x1).sum()*(x1[1]-x1[0])
        S2=fd(x2).sum()*(x2[1]-x2[0])
        print('{} {} S1: {:8.1e}, S2: {:8.1e}'.format(i,d,1-S1,1-S2))
        

plt.draw()
plt.waitforbuttonpress()
quit(1)
  

for i in [1,2]:
    for k in range(nc):
        d=de.Domain([c[k]], np.float64)
        b=d.bases[0]
        x1 = b.grid(scale=1)
        x2 = b.grid(scale=A)
        #print('k {}, x1[0] {} x1[-1] {}, N:, len(x1): {}'.format(k,N*(x1[0]/L),N*(1-x1[-1]/L),N,len(x1)))
        fg, a0=plt.subplots(nrows=nM,ncols=5)
        g   = d.new_field()
        dg  = d.new_field()
        gs  = d.new_field()
        dgs = d.new_field()
        if k==1:
            g.meta['x']['parity']   = -1
            dg.meta['x']['parity']  =  1
            gs.meta['x']['parity']  =  1
            dgs.meta['x']['parity'] =  1
        else:
            g.meta['x']['parity']   =  0
            dg.meta['x']['parity']  =  0
            gs.meta['x']['parity']  =  0
            dgs.meta['x']['parity'] =  0
        for j in range(nM):
            FFTA=ffta.ffta(b,MM[j],nfx)
            f1  = FFTA.f(x-tol/2,i)
            df1 = FFTA.df(x-tol/2,i)
            f2  = FFTA.f(x+tol/2,i)
            df2 = FFTA.df(x+tol/2,i)
            f   = FFTA.f(x,i)
            df  = FFTA.df(x,i)
            dfc = (f2-f1)/tol
            g.set_scales(1)
            dg.set_scales(1)
            gs.set_scales(1)
            dgs.set_scales(1)
            (aa,bb)=FFTA.dfsl(x1,i)
            if k==0:
                g['g']   = FFTA.f(x1,i)  - FFTA.f(x1-L/2,i)  + FFTA.f(x1-L,i)
                dg['g']  = FFTA.df(x1,i) - FFTA.df(x1-L/2,i) + FFTA.df(x1-L,i)
                g2       = FFTA.f(x2,i)  - FFTA.f(x2-L/2,i)  + FFTA.f(x2-L,i)
                dg2      = FFTA.df(x2,i) - FFTA.df(x2-L/2,i) + FFTA.df(x2-L,i)
            elif k==1:
                g['g']   = FFTA.f(x1,i)  - FFTA.f(x1-L,i)  - 1
                dg['g']  = FFTA.df(x1,i) - FFTA.df(x1-L,i) 
                g2       = FFTA.f(x2,i)  - FFTA.f(x2-L,i)- 1
                dg2      = FFTA.df(x2,i) - FFTA.df(x2-L,i)
            elif k==2:
                g['g']   = FFTA.f(x1,i)
                dg['g']  = FFTA.df(x1,i)
                g2       = FFTA.f(x2,i)
                dg2      = FFTA.df(x2,i)
            gs['g']  = g['g']**2
            dgs['g'] = dg['g']**2
            gs2      = g2**2
            dgs2     = dg2**2
            
            MMMM=i*np.int((MM[j]-1)/2)+1
            MMM=np.asarray(np.where(g['g'] == 1)).min()
            #if k==2: MMM=MM[j]
            #print('M {}, MMM {}'.format(MM[j],MMM))
            a0[j,0].plot(x1[0:MMM],g['g'][0:MMM], marker='s')
            a0[j,1].plot(x1[0:MMM],dg['g'][0:MMM], marker='s')
            a0[j,2].plot(x1,g['g'], marker='s')
            a0[j,3].plot(x1,dg['g'], marker='s')
            
            g.set_scales(A)
            dg.set_scales(A)
            gs.set_scales(A)
            dgs.set_scales(A)
            a0[j,2].plot(x2,g['g'])
            a0[j,3].plot(x2,dg['g'])
            a0[j,2].plot(x2,g2)
            a0[j,3].plot(x2,dg2)
            a0[j,4].plot(x,df-dfc)  

            e1=max(abs(df-dfc))
            e2=max(abs(g['g']-g2))
            e3=max(abs(dg['g']-dg2))
            e4=max(abs(gs['g']-gs2))
            e5=max(abs(dgs['g']-dgs2))

            print('i: {}, M: {:2d}, M1: {:2d}, M2: {:2d}, err: {:8.1e} {:8.1e} {:8.1e} {:8.1e} {:8.1e}'.format(i,MM[j],MMM,MMMM,e1,e2,e3,e4,e5))
  
nc=3
c=list()
c.append(de.SinCos(    'x', N, interval = LL, dealias=A))
c.append(de.Fourier(   'x', N, interval = LL, dealias=A))
c.append(de.SinCos(    'x', N, interval = LL, dealias=A))

kmax = np.ndarray((nc,),dtype=np.float64)
P    = np.ndarray((nc,),dtype=np.int64)
d=list()
for k in range(nc):
    P[k]=k-1;
    d.append(de.Domain([c[k]], np.float64))
    ff.append(d[k].new_field())
    ff[k].meta['x']['parity']=P[k]
    fs.append(d[k].new_field())
    fs[k].meta['x']['parity']=P[k]
    w.append(np.abs(d[k].bases[0].wavenumbers))
    kmax[k]=max(w[k])

MM=range(4,12)
n=len(MM)
tol=1e-5;
fg, a0=plt.subplots(nrows=n,ncols=1)
for j in range(n):
    FFTA=ffta.ffta(b,2*MM[j]-5,MM[j])
    f=FFTA.constant(x2)
    a0[j].plot(x2,f)
    a0[j].plot(x2,FFTA.odd_constant(x2))
  
plt.draw()
plt.waitforbuttonpress()
quit(1)



kmax=np.ndarray((nc,),dtype=np.float64)
P=np.ndarray((nc,),dtype=np.int64)
for k in range(nc):
    P[k]=k-1;
    d.append(de.Domain([c[k]], np.float64))
    ff.append(d[k].new_field())
    ff[k].meta['x']['parity']=P[k]
    fs.append(d[k].new_field())
    fs[k].meta['x']['parity']=P[k]
    w.append(np.abs(d[k].bases[0].wavenumbers))
    kmax[k]=max(w[k])



float_formatter = lambda x: '{:8.1e}'.format(x)
np.set_printoptions(formatter={'float_kind':float_formatter})


M=6
fg, a0=plt.subplots(nrows=M,ncols=3)
fg, a1=plt.subplots(nrows=M,ncols=2)

tol=1e-5;
for j in range(M):
    FFTA=ffta.ffta(b,5,M)
    [f2,df2,df0,If2] = FFTA.heaviside(x2+tol/2)
    [f1,df1,df0,If1] = FFTA.heaviside(x2-tol/2)
    [f,df,df0,If]=     FFTA.heaviside(x2)
    a0[j,0].plot(x,f)
    a0[j,1].plot(x,df)
    e0=max(np.maximum(0,f-1))
    e1=max(np.maximum(0,1-f))
    e2=df.sum()*dx-2
    e3=max(np.maximum(0,-df))
    print('M: {:2d}, {:8.1e} {:8.1e} {:8.1e} {:8.1e}'.format(j,e0,e1,e2,e3))
    a0[j,2].plot(x,If)

    a1[j,0].plot(x,f-(If2-If1)/tol)
    a1[j,1].plot(x,df-(f2-f1)/tol)
    
    
plt.draw()
plt.waitforbuttonpress()
quit(1)

MM=np.arange(nc,7)
e1=np.ndarray((nc,))
e2=np.ndarray((nc,))
e3=np.ndarray((nc,))
YY=[0,0.1,0.3,0.5]

def Pext(f,P):
    if   P==0: f=np.concatenate((f,f,f,f))
    elif P==1: f=np.concatenate((np.flip( f),f,np.flip( f),f))
    else:      f=np.concatenate((np.flip(-f),f,np.flip(-f),f))
    return(f)

fg, ax=plt.subplots(nrows=len(YY),ncols=6)
fg, bx=plt.subplots(nrows=len(YY),ncols=6)
fg, cx=plt.subplots(nrows=len(YY),ncols=6)
for k in range(len(YY)):
    y=YY[k]
    ax[k,0].set_ylabel('y {:4.2f}'.format(y))
    for j in range(nc):
        x1 = c[j].grid(scale=1)
        x2 = c[j].grid(scale=A)
        x1 = c[j].grid(scale=1)
        x2 = c[j].grid(scale=A)
        f=ff[j]
        h=fs[j]
        f.set_scales(1)
        h.set_scales(1)
        z=ju.ffta_delta(x1,y,P[j],kmax[j]/2,M,L)
        f['g']=z
        h['g']=z**2
        fc=copy.copy(np.abs(f['c']))
        hc=copy.copy(np.abs(h['c']))
        f.set_scales(A)
        h.set_scales(A)
        f1=copy.copy(f['g'])
        f2=copy.copy(h['g'])
        f1=Pext(f1,P[j])
        f2=Pext(f2,P[j]**2)
        x2=np.concatenate((x2-L,x2,x2+L,x2+2*L))
        
        #print('j={}, w.shape {}, fc.shape {}, hc.shape {}'.format(j,w[j].shape,fc.shape,hc.shape))
        ax[k,j  ].plot(x2,f1)
        ax[k,j+3].plot(x2,f2)
        bx[k,j  ].plot(w[j]/kmax[j],fc)
        bx[k,j+3].plot(w[j]/kmax[j],hc)

        ax[k,j  ].set_xlim((-L,3*L))
        ax[k,j+3].set_xlim((-L,3*L))
        bx[k,j  ].set_xlim((0,1))
        bx[k,j+3].set_xlim((0,1))

        ax[0,j  ].set_title('P={:2d}'.format(P[j]))
        ax[0,j+3].set_title('P={:2d}'.format(P[j]))
        bx[0,j  ].set_title('P={:2d}'.format(P[j]))
        bx[0,j+3].set_title('P={:2d}'.format(P[j]))

        f.set_scales(1)
        f['g']=ju.ffta_interval(x1,y,y+L/5,P[j],kmax[j],M,L)
        f.set_scales(A)
        cx[k,j].plot(x2,Pext(f['g'],P[j]))
        f.set_scales(1)
        f['g']=ju.ffta_interval(x1,y+L/5,y,P[j],kmax[j],M,L)
        f.set_scales(A)
        cx[k,j+3].plot(x2,Pext(f['g'],P[j]))
    
    e3[j]=f1.sum()/NN/4-1
    e1[j]=max(f1.max()-1,-f1.min())
    e2[j]=max(f2.max()-1,-f2.min())
    print('M: {:2d}, I=[{:8.1e} {:8.1e}] e1=[{:8.1e} {:8.1e}] e2=[{:8.1e} {:8.1e}]'.format(M,e3[0],e3[1],e1[0],e1[1],e2[0],e2[1]))

plt.draw()
plt.waitforbuttonpress()





