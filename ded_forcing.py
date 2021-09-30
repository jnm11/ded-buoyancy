#!/usr/bin/env python3
import jutil_131 as ju
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
plt.ion()

L=8.
H=4.
Nx=769
Nz=128
h=1.0
Wz=0.5

x=np.linspace(0, L, Nx, endpoint=True, dtype=np.float) 
z=np.linspace(0, H, Nz, endpoint=True, dtype=np.float) 



h1=ju.cubici(z,h+Wz/2,h-Wz/2)

fgx, ax=plt.subplots(nrows=3,ncols=2,figsize=(10,5))
fgz, az=plt.subplots(nrows=3,ncols=2,figsize=(5,10))
fgf, af=plt.subplots(nrows=3,ncols=4,figsize=(8,10))

fuv, auv=plt.subplots(nrows=8,ncols=1,figsize=(8,8))


nfx=2;
nfz=np.int(np.round(Nz*h/(2*H)));

x0=0.0
x1=0.5
x2=1.5
x3=2.0

x4=1
x5=1.5
x6=2.5
x7=3

Wx=0.25
U=1
U1=0.2

# We add some extra terms to the x forcing to make sure it it is fully periodic
# These only occur where the weight is zero however
dx12=x2-x1
dx23=x3-x2
fx1,Fx1,dfx1 = ju.smoothi(x,nfx,x1,x2)
fx2,Fx2,dfx2 = ju.smoothi(x,nfx,x3,x3+dx12)
fx3,Fx3,dfx3 = ju.smoothi(x,nfx,x3+dx12+dx23,x3+2*dx12+dx23)
fz,Fz,dfz    = ju.smoothi(z,nfz, 0, 2*h)
fzH,FzH,dfzH = ju.smoothi(H,nfz, 0, 2*h)
print('nfz: {} fzH: {} FzH: {} (H-h): {} dfzH: {} '.format(nfz,fzH,(H-h),FzH,dfzH ))
dfx = dfx1 - 2*dfx2 + dfx3
fx  = fx1   -2*fx2  + fx3
Fx  = Fx1   -2*Fx2  + Fx3
wu=ju.cubici(x,x0,x0+Wx)-ju.cubici(x,x2,x3)       # Velocity weight
wb=ju.cubici(x,x0,x0+Wx)-ju.cubici(x,x6,x7)       # Density  weight

#Fz=Fz*H/FzH
#fz=fz*H/FzH
#dfz=dfz*H/FzH

zz   = np.tile(np.reshape(z, (1,Nz)),(Nx,1))
Fzz  = np.tile(np.reshape(Fz,(1,Nz)),(Nx,1))
fzz  = np.tile(np.reshape(fz,(1,Nz)),(Nx,1))
dfzz  = np.tile(np.reshape(dfz,(1,Nz)),(Nx,1))
fxx  = np.tile(np.reshape(fx,(Nx,1)),(1,Nz))
dfxx = np.tile(np.reshape(dfx,(Nx,1)),(1,Nz))
psi = -zz + (1+U1/U)*fxx* (zz - H/(H-h)*Fzz)
u   = -1  + (1+U1/U)*fxx* (1  - H/(H-h)*fzz)
v   =     - (1+U1/U)*dfxx*(zz - H/(H-h)*Fzz)
du  =     + (1+U1/U)*fxx* (   - H/(H-h)*dfzz)

print('shape(psi) {}'.format(psi.shape))

bfx = ju.cubici(x,x4,   x5)-ju.cubici(x,x7,x7+x5-x4) # Bouyancy field
bfz = ju.cubici(z,3/2*h,h)

dfx = dfx/max(abs(dfx))
fx  = fx/max(abs(fx))
Fx  = Fx/max(abs(Fx))

bfx = bfx/max(abs(bfx))
wb  = wb/max(abs(wb))
wu  = wu/max(abs(wu))

xrg1=0
xrg2=L
yrg1=-1.1
yrg2=+1.1


ax[0,0].plot(x,dfx)
ax[0,0].set_xlabel('x')
ax[0,0].set_ylabel('dfx')
ax[1,0].plot(x,fx)
ax[1,0].set_xlabel('x')
ax[1,0].set_ylabel('fx')
ax[2,0].plot(x,Fx)
ax[2,0].set_xlabel('x')
ax[2,0].set_ylabel('Fx')

ax[0,1].plot(x,bfx)
ax[0,1].set_xlabel('x')
ax[0,1].set_ylabel('bfx')
ax[1,1].plot(x,wb)
ax[1,1].set_xlabel('x')
ax[1,1].set_ylabel('wb')
ax[2,1].plot(x,wu)
ax[2,1].set_xlabel('x')
ax[2,1].set_ylabel('wu')

for j in range(0,3):
   for k in range(0,2):
        ax[j,k].set_xlim(xrg1,xrg2)
        ax[j,k].set_ylim(yrg1,yrg2)
        if k==0: ax[j,k].add_patch(patches.Rectangle((x0,yrg1),x3-x0,yrg2-yrg1,linewidth=1,edgecolor='none',facecolor=(0.8,0.7,0.7,0.5)))
        if k==1: ax[j,k].add_patch(patches.Rectangle((x0,yrg1),x7-x0,yrg2-yrg1,linewidth=1,edgecolor='none',facecolor=(0.7,0.7,0.9,0.5)))
        ax[j,k].plot([x0,x0],[yrg1,yrg2],'r')
        ax[j,k].plot([x1,x1],[yrg1,yrg2],'r')
        ax[j,k].plot([x2,x2],[yrg1,yrg2],'r')
        ax[j,k].plot([x3,x3],[yrg1,yrg2],'r')
        ax[j,k].plot([x4,x4],[yrg1,yrg2],'r')
        ax[j,k].plot([x5,x5],[yrg1,yrg2],'r')
        ax[j,k].plot([x6,x6],[yrg1,yrg2],'r')
        ax[j,k].plot([x7,x7],[yrg1,yrg2],'r')
 
        
az[0,0].plot(dfz,z)
az[0,0].set_xlabel('dfz')
az[0,0].set_ylabel('z')
az[1,0].plot(fz,z)
az[1,0].set_xlabel('fz')
az[1,0].set_ylabel('z')
az[2,0].plot(Fz,z)
az[2,0].set_xlabel('Fz')
az[2,0].set_ylabel('z')

az[0,1].plot(h1,z)
az[0,1].set_xlabel('h1')
az[0,1].set_ylabel('z')

for k in range(0,4):
    n=int(((3-k)*x1+k*x2)*(Nx-1)/L/3)
    af[0,k].plot(psi[n,:],z)
    af[0,k].set_xlabel('psi {:7.5f} '.format(x[n]))
    af[0,k].set_ylabel('z')
    af[1,k].plot(u[n,:],z)
    af[1,k].set_xlabel('u {:7.5f} '.format(x[n]))
    af[1,k].set_ylabel('z')
    af[2,k].plot(v[n,:],z)
    af[2,k].set_xlabel('v {:7.5f} '.format(x[n]))
    af[2,k].set_ylabel('z')


n1=int(1.05*h/H*(Nz-1))
n2=int(0.95*h/H*(Nz-1))


auv[0].plot(x,v[:,-1])
auv[0].set_xlabel('x')
auv[0].set_ylabel('v top')
auv[1].plot(x,v[:,n2])
auv[1].set_xlabel('x')
auv[1].set_ylabel('v 1/3')
auv[2].plot(x,v[:,n1])
auv[2].set_xlabel('x')
auv[2].set_ylabel('v 2/3')
auv[3].plot(x,v[:,0])
auv[3].set_xlabel('x')
auv[3].set_ylabel('v bottom')
auv[4].plot(x,du[:,-1])
auv[4].set_xlabel('x')
auv[4].set_ylabel('du top')
auv[5].plot(x,du[:,n2])
auv[5].set_xlabel('x')
auv[5].set_ylabel('du 1/3')
auv[6].plot(x,du[:,n1])
auv[6].set_xlabel('x')
auv[6].set_ylabel('du 2/3')
auv[7].plot(x,du[:,0])
auv[7].set_xlabel('x')
auv[7].set_ylabel('du bottom')

fig1, ax1 = plt.subplots()
ax1.contour(x,z,np.transpose(psi),50)
plt.draw()
plt.waitforbuttonpress()
