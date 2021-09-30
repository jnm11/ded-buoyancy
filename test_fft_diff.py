#!/usr/bin/env python3
#import matplotlib
#matplotlib.use("TkAgg")
# test_alias
#import matplotlib
#matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import scipy.fftpack as fft
import jutil as ju

plt.ion()




L=1

M=100
x  = L*(np.arange(M)+0.5)/M*2*math.pi
y  = L*(np.arange(M))/M*2*math.pi
dx=x[1]-x[0]
dy=y[1]-y[0];

for n in range(6):
    w1=(n+1/2)/L
    w2=(n+1)/L
    s=np.sin(w1*x)
    c=np.cos(w1*x)
    f=np.sin(w2*y+0.2)
    
    ds=+w1*np.cos((n+1/2)*x/L)
    dc=-w1*np.sin((n+1/2)*x/L)
    df=+w2*np.cos((n+1  )*y/L+0.2)

    dsa=ju.jfft_diff(s,dx,-1)
    dca=ju.jfft_diff(c,dx, 1) 
    dfa=ju.jfft_diff(f,dy, 0)

    fg, a=plt.subplots(nrows=3,ncols=2)                    
    
    a[0,0].plot(x,ds,x,dsa)
    a[1,0].plot(y,df,y,dfa)
    a[2,0].plot(x,dc,x,dca)

    a[0,1].plot(x,ds-dsa)
    a[1,1].plot(y,df-dfa)
    a[2,1].plot(x,dc-dca)

    print('sin {}'.format((ds*dsa).sum()/(dsa*dsa).sum()))
    print('cos {}'.format((dc*dca).sum()/(dca*dca).sum()))
    print('fft {}'.format((df*dfa).sum()/(dfa*dfa).sum()))
    
    #print('s0.max  {:7.3f}'.format(np.abs(s0).max()))
    #print('s0a.max {:7.3f}'.format(np.abs(s0a).max()))
    #print('s0b.max {:7.3f}'.format(np.abs(s0b).max()))
    #print('s0c.max {:7.3f}'.format(np.abs(s0c).max()))
    #print('c0.max  {:7.3f}'.format(np.abs(c0).max()))
    #print('c0a.max {:7.3f}'.format(np.abs(c0a).max()))
    #print('c0b.max {:7.3f}'.format(np.abs(c0b).max()))
    #print('c0c.max {:7.3f}'.format(np.abs(c0c).max()))
    
plt.draw()
plt.waitforbuttonpress()
#quit(1)
