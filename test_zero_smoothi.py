#!/usr/bin/env python3
import jutil as ju
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
plt.ion()


x=np.linspace(-3, 3, num=1e4, endpoint=True, dtype=np.float)    
N1=5
N2=5
for n in range(N1,1+N2):
    fg, ((ax1,ax4),(ax2,ax5),(ax3,ax6))=plt.subplots(nrows=3,ncols=2)
    f,F,df=ju.smoothi(x,n,-1,1)
    d=1e-5;
    f1,F1,df1=ju.smoothi(x+2*d,n,-1,1)
    f2,F2,df2=ju.smoothi(x+1*d,n,-1,1)
    f3,F3,df3=ju.smoothi(x-1*d,n,-1,1)
    f4,F4,df4=ju.smoothi(x-2*d,n,-1,1)
    e1=f -(-F1+8*F2-8*F3+F4)/(12*d)
    e2=df-(-f1+8*f2-8*f3+f4)/(12*d)
    ax1.plot(x,df)
    ax2.plot(x,f)
    ax3.plot(x,F)
    ax4.plot(x,e2)
    ax5.plot(x,e1)


plt.draw()
plt.waitforbuttonpress()

