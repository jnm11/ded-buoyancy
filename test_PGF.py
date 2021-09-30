#!/usr/bin/env python3
import jutil as ju
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math

if 0:
    dd=np.array([ 0.6,  0.6,  0.6,  0.6])
    yy=np.array([+0.1, +0.1, -0.1, -0.1])
    zz=np.array([+0.1, -0.1, -0.1, +0.1])
    h=1
    W=3
    H=3
    for j in range(len(dd)):
        d=dd[j]
        y=yy[j]
        z=zz[j]
        
        (fy,fz)=ju.periodic_gaussian(d,y,z,h,W,H)
        print('d={:3.1f}, y={:4.1f}, z={:4.1f}, h={:3.1f}, w={:1.0f}, H={:1.0f}, fy={:10.7f}, fz={:10.7f}'.format(d,y,z,h,W,H,fy,fz))
        
plt.ion()



Nx=100
Ny=104
Nz=102

W=np.float64(20)
H=np.float64(21)
h=np.float64(1)

print('W={}, Ny={}'.format(W,Ny))
print('H={}, Nz={}'.format(H,Nz))
print('h={}'.format(h))

tol=np.float64(1e-4)
h=np.float64(1)
nd=10;
dd=np.linspace(0,1,nd)**4/h
y=np.linspace(-W/2, W/2, num=Ny, endpoint=True, dtype=np.float64)    
z=np.linspace(-H/2, H/2, num=Nz, endpoint=True, dtype=np.float64)    
dy=W/(Ny-1)
dz=H/(Nz-1)

yy,zz = np.meshgrid(y,z)
Q1=np.ndarray(nd)
Q2=np.ndarray(nd)
M=np.ndarray(nd)
u0=1;
A0=h**2*math.pi
Q0=A0*u0
M0=A0*u0**2

for j in range(len(dd)):
    d=dd[j]
    fy,fz,c = ju.periodic_gaussian(d,yy,zz,    h,W,H)
    fy1,fz1,c = ju.periodic_gaussian(d,yy+tol,zz,    h,W,H)
    fy2,fz2,c = ju.periodic_gaussian(d,yy-tol,zz,    h,W,H)
    fy3,fz3,c = ju.periodic_gaussian(d,yy    ,zz+tol,h,W,H)
    fy4,fz4,c = ju.periodic_gaussian(d,yy    ,zz-tol,h,W,H)
    u=u0*(fy1-fy2+fz3-fz4)/(2*tol)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(zz,yy, u, color='b')
    v=u**2
    Q1[j]  = u.sum(axis=None)*dy*dz
    M[j]   = v.sum(axis=None)*dy*dz
    print('fy.shape {}'.format(fy.shape))
    print('fz.shape {}'.format(fz.shape))
    Q2[j]  = u0*(np.sum(fz[Nz-1,:]-fz[0,:],axis=0)*dy+np.sum(fy[:,Ny-1]-fy[:,0],axis=0)*dz)
    plt.draw()
    print('h={:8.6f}, d={:8.6f}, Q1={:8.6f}, Q2={:8.6f}, M={:8.6f}'.format(h,d,Q1[j]/Q0,Q2[j]/Q0,M[j]/M0))

print('Q1 {}'.format(Q1/Q0))
print('M {}'.format(M/M0))


plt.waitforbuttonpress()


#y=np.array([0, 1])*W/4;
#z=np.array([0, 1])*H/4;
#yy,zz,xx = np.meshgrid(y,z,x)
#dd=ju.cubici(xx,0.5,1.5)/h
#d=ju.cubici(x,0.5,1.5)/h

#fy,fz = ju.periodic_gaussian(dd,yy,zz,h,W,H)

#fy1,fz1 = ju.periodic_gaussian(dd,yy+tol,0,h,W,H)
#fy2,fz2 = ju.periodic_gaussian(dd,yy-tol,0,h,W,H)
#fy3,fz3 = ju.periodic_gaussian(dd,yy,zz+tol,h,W,H)
#fy4,fz4 = ju.periodic_gaussian(dd,yy,zz-tol,h,W,H)
#u=(fy1-fy2+fz3-fz4)/(2*tol)

#fy=fy.reshape(4,Nx).transpose()
#fz=fz.reshape(4,Nx).transpose()
#xx=xx.reshape(4,Nx).transpose()
#u=u.reshape(4,Nx).transpose()
#fg, ax=plt.subplots(nrows=3,ncols=2)
#ax[0,0].plot(x,d)
#ax[1,0].plot(xx,fy)
#ax[2,0].plot(xx,fz)
#ax[0,1].plot(xx,u)
