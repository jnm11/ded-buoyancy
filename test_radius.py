#!/usr/bin/env python3
#import matplotlib
#matplotlib.use("TkAgg")
# test_alias
"""
Test calculation of part circle radius given centre

Usage: 
   test_alias.py [options]

Options:
  -R, --Radius=<>  : Length of domain       [default:   2.0]
  -Y, --YCentre=<> : Number of grid points  [default:   0.1]
  -Z, --ZCentre=<> : Width of front region  [default:  -0.1]

"""

from docopt import docopt
import numpy as np
import scipy.special as sp
import jutil as ju
import math
import matplotlib.pyplot as plt
#plt.ion()

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

R=np.float64(args['--Radius'])
Y=np.float64(args['--YCentre'])
Z=np.float64(args['--ZCentre'])

fg, ax=plt.subplots(nrows=1,ncols=1)

Py=0
Pz=0
AA=1/(1+Py)/(1+Pz)
for k in range(10):
    y=np.random.normal()
    z=np.random.normal()
    r=np.random.uniform(low= 0, high=1)
    r1=r/2
    r2=r+np.sqrt(y**2+z**2)
    R=ju.get_radius(r,y,z,Py,Pz)
    A=ju.circle_area(R,-y,-z,Py,Pz)
    print('A={:6.4f}, R={:6.4f}, r={:6.4f}, z={:6.3f}, y={:6.3f}, '.format(A/(math.pi*r**2),R,r,z,y))
    
    #R=np.linspace(r1, r2, 1000)
    #A=ju.circle_area(R,y,z,Py,Pz)
    #ax.plot(R/r,A/math.pi-r**2*AA)
quit(1)
plt.draw()
plt.waitforbuttonpress()

    
    
fg, (ax1,ax2)=plt.subplots(nrows=2,ncols=1)

Y=np.linspace(-1.1, 1.1, 1000)
for Z in np.linspace(-1.1, 1.1, 7):
    A=ju.circle_area(1,Y,Z,1,1)
    ax1.plot(Y,A/math.pi)
ax1.plot(Y,ju.jasin(Y)/math.pi)

f=np.sin
for k in range(10):
    x=np.random.uniform(low=-math.pi/2, high=math.pi/2)
    a=f(x)
    y=ju.secant(f,a,-math.pi/2,math.pi/2)
    print('a={:7.4f}, x={:7.4f}, y={:7.4f}, abs(x-y)={:8.2e}, abs(f(y)-a)={:8.2e}'.format(a,x,y,abs(x-y),abs(f(y)-a)))
    
s=np.linspace(-1.1, 1.1, 1000)
     
#ax1.plot(WW,np.log10(tmax),WW,np.log10(emax),WW,np.log10(smax))
#ax2.plot(WW,np.log10(tint),WW,np.log10(eint),WW,np.log10(sint))
ax2.plot(s, ju.jasin(s)/math.pi)


plt.draw()
plt.waitforbuttonpress()
