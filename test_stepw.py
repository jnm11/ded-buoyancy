#!/usr/bin/env python3
import jutil as ju
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
plt.ion()

L=2
x=np.linspace(-L/2, 3*L/2, num=1e5, endpoint=True, dtype=np.float)    
x1=0
x2=1*L/10
x3=5*L/10
x4=6*L/10
n=5

fg, ((ax1,ax4),(ax2,ax5),(ax3,ax6))=plt.subplots(nrows=3,ncols=2)
f,F,df=ju.stepw(x,n,x1,x2,x3,x4,L)
d=1e-4;
f1,F1,df1=ju.stepw(x+2*d,n,x1,x2,x3,x4,L)
f2,F2,df2=ju.stepw(x+1*d,n,x1,x2,x3,x4,L)
f3,F3,df3=ju.stepw(x-1*d,n,x1,x2,x3,x4,L)
f4,F4,df4=ju.stepw(x-2*d,n,x1,x2,x3,x4,L)
e1=f -(-F1+8*F2-8*F3+F4)/(12*d)
e2=df-(-f1+8*f2-8*f3+f4)/(12*d)
ax1.plot(x,df)
ax1.set_ylabel('df')
ax2.plot(x,f)
ax2.set_ylabel('f')
ax3.plot(x,F)
ax3.set_ylabel('F')
ax4.plot(x,e2)
ax4.set_ylabel('df error')
ax5.plot(x,e1)
ax5.set_ylabel('F error')
ax6.plot(x,0.5-abs(f-0.5))
ax6.set_ylabel('range')
ax6.set_ylim([0, 0.0001])
plt.draw()
plt.waitforbuttonpress()

