import numpy as np
import scipy.special as sp
import jutil as ju
import math
from dedalus import public as de
from dedalus import core as co

# Fourier series for basic functions with appropriate parity
# The series are chosen so as to have no ringing



# Heaviside function located as close as possible to the origin
# f(x0)=1 f^(k)(x0)=0 for k=0..2M
# f^(k)(x1)=0 k=0..2(N-M)
# int_x^x1 sin(x)^n*cos(x/2)^(2*(N-n)) dx # for even function

def heaviside(f,domain,d,M=1):
    dim     = domain.dim
    gcshape = domain.dist.coeff_layout.global_shape(scales=1)
    lgshape = domain.dist.grid_layout.local_shape(scales=1)
    ggshape = domain.dist.grid_layout.global_shape(scales=1)
    lslice  = domain.dist.coeff_layout.slices(scales=1)
    if 'parity' in f.meta[d]: parity=f.meta[d]['parity']
    else:                     parity=0
    for j in range(dim):
        if domain.bases[j].name==d: break
    basis=domain.bases[j]
    if basis.name !=d: ju.logger.info('Could not find the write direction')
    L = basis.interval[1]-basis.interval[0]
    w = basis.wavenumbers
    #ju.logger.info('sawtooth: P: {} d: {} w: {}'.format(parity,d,L*w/math.pi))
    cslice = lslice[j]
    Nc = gcshape[j]
    Ng = ggshape[j]
    #ju.logger.info('dftf.py.heaviside: Nc: {:4d}, Ng: {:4d}, Parity: {:2d}, L: {}, j: {}, d: {}, gsh: {}, csh: {}'.format(Nc,Ng,parity,L,j,d,ggshape,gcshape))
    if parity!=-1:
        if parity==0:     # periodic function
            k=np.floor(L/(2*math.pi)*w[cslice])
            N=ju.jsum(w>0);
            A = 2
        elif parity==1:   # even function
            k = np.arange(cslice.start,cslice.stop)
            N = Nc-1
            A = 4
        A *= np.exp(-N*np.log(4)+sp.gammaln(2*N+1)-sp.gammaln(N+1-k)-sp.gammaln(N+1+k))
        if   M==2: A *=          (3*N          -2*k**2        -1)/   (2*N-1)
        elif M==3: A *= (15*N**2-25*N+4*k**4-20*N*k**2+20*k**2+6)/(2*(2*N-1)*(2*N-3))
        elif M==4: A *= (105*N**3+(-210*k**2-420)*N**2+(84*k**4+630*k**2+441)*N-8*k**6-140*k**4-392*k**2-90)/(2*(2*N-1)*(2*N-3)*(2*N-5))
            
        if parity==1: A[k==0]=A[k==0]/2
  
        A=A*L/math.pi
        f['c']=0
        if dim==3:
            if   j==0 and lslice[1].start==0 and lslice[2].start==0: f['c'][:,0,0]=A
            elif j==1 and lslice[0].start==0 and lslice[2].start==0: f['c'][0,:,0]=A
            elif j==2 and lslice[0].start==0 and lslice[1].start==0: f['c'][0,0,:]=A
            
        elif dim==2:
            if   j==0 and lslice[1].start==0: f['c'][:,0]=A
            elif j==1 and lslice[0].start==0: f['c'][0,:]=A
        else:
            f['c']=A
    return(f)

# Calculate a sawtooth function centred in the domain
# if m is the  midpoint then
# f(m)=1
# f^(n)(m)=1 for all m>0
# f(x1)=f(x2)=0 at the end points
def sawtooth(f,domain,d,AA=1):
    dim     = domain.dim
    gcshape = domain.dist.coeff_layout.global_shape(scales=AA)
    lgshape = domain.dist.grid_layout.local_shape(scales=AA)
    ggshape = domain.dist.grid_layout.global_shape(scales=AA)
    lslice  = domain.dist.coeff_layout.slices(scales=AA)
    if 'parity' in f.meta[d]: parity=f.meta[d]['parity']
    else:                     parity=0
    for j in range(dim):
        if domain.bases[j].name==d: break
    basis=domain.bases[j]
    if basis.name !=d: ju.logger.info('Could not find the right direction')
    L = basis.interval[1]-basis.interval[0]
    w = basis.wavenumbers
    #ju.logger.info('sawtooth: P: {} d: {} w: {}'.format(parity,d,L*w/math.pi))
    cslice = lslice[j]
    Nc = gcshape[j]
    Ng = ggshape[j]
    #ju.logger.info('dftf.py.sawtooth: Nc: {:4d} Ng: {:4d} Parity: {:2d} L: {} j: {} d: {} gsh: {} csh: {}'.format(Nc,Ng,parity,L,j,d,ggshape,gcshape))
    if parity==0:     # periodic function
        k=np.floor(L/(2*math.pi)*w[cslice])
        N=ju.jsum(w>0);
        A = 1j*np.exp(2*sp.gammaln(N)-sp.gammaln(k+N)-sp.gammaln(N-k))/(2*k)
    elif parity==1:   # even function
        k = np.arange(cslice.start,cslice.stop)
        N = (Nc-1)/2
        k = (k+1)/2
        A = np.where(k==np.floor(k),-np.exp(2*sp.gammaln(2*N)-2*np.log(2*k-1)-sp.gammaln(N+1-k)-sp.gammaln(k+N)-(N-1)*np.log(16)-2*sp.gammaln(N)),0)
    elif parity==-1:  # odd function
        k = np.arange(cslice.start,cslice.stop)
        #N = (Nc-1)/2
        #k = k/2
        #A = np.where(k==np.floor(k),-2*sp.gammaln(N+1)-sp.gammaln(N+k+1)-sp.gammaln(1+N-k)-np.log(k)),0)
        N=Nc;
        A=         -2*(-1)**k/k*np.exp(2*sp.gammaln(N+1)-sp.gammaln(N+k+1)-sp.gammaln(1+N-k));
    A[k==0]=0
    A=A*L/math.pi
    f['c']=0
    if dim==3:
        if   j==0 and lslice[1].start==0 and lslice[2].start==0: f['c'][:,0,0]=A
        elif j==1 and lslice[0].start==0 and lslice[2].start==0: f['c'][0,:,0]=A
        elif j==2 and lslice[0].start==0 and lslice[1].start==0: f['c'][0,0,:]=A
        
    elif dim==2:
        if   j==0 and lslice[1].start==0: f['c'][:,0]=A
        elif j==1 and lslice[0].start==0: f['c'][0,:]=A
    else:
        f['c']=A

    return(f)

# Set the coefficients in the field xb to have a single hump at y of height 1
# I is the integral of the field
def ffta_single(f,domain,d,y):
    if 'parity' in f.meta[d]: parity=f.meta[d]['parity']
    else: parity=0
    dim     = domain.dim
    for k in range(dim):
        if domain.bases[k].name==d: break
    basis=domain.bases[k]
    if basis.name !=d:
        ju.logger.info('Could not find the write direction')
    cslice  = domain.dist.coeff_layout.slices(scales=1)[k]
    lcshape = domain.dist.coeff_layout.local_shape(scales=1)[k]
    N       = domain.dist.coeff_layout.global_shape(scales=1)[k]-1
    w       = basis.wavenumbers[cslice]
    i=np.arange(cslice.start,cslice.stop)
    sz=np.ones((dim,),dtype=np.int64)
    sz[k]=lcshape
    if parity==0: # periodic function
        f['c'] = np.reshape(np.exp(-1j*w*y)*A,sz)
        I=A[0]
    elif parity==1: # even function
        f['c'] = np.reshape(2*np.cos(w*y)*A/(1+np.heaviside(-i,1)),sz)
        I=A[0]
    elif parity==-1: # odd function
        f['c'] = np.reshape(4*np.sin(w*y)*A,sz)
        I=2*np.sum(f['c'][1:N+1:2]/w[1:N+1:2])
    return(f,I)




# x and y should be scales to -Pi..Pi
# This is chosen to have the first m normal derivatives vanish at the boundary where it takes a value of 1
# It is maximally "round" in the middle where it takes a value of 0
# We then raise this to a large power while satisfying the aliasing criterion 
# Point decreases linearly with n
def heaviside2(x,y,N,m=6,k=1):
    n=np.int(np.floor(N/m))
    #ju.logger.info('heaviside2: N:{} n:{}'.format(N,n))
    #ju.srange('heaviside2: x range [{:6.3f},{:6.3f}]',x/math.pi)
    #ju.srange('heaviside2: y range [{:6.3f},{:6.3f}]',y/math.pi)
    x=np.minimum(math.pi,np.abs(x))
    y=np.minimum(math.pi,np.abs(y))
    if   m==1:
        c1y=np.cos(1*y)
        c2y=np.cos(2*y)
        c3y=np.cos(3*y)
        f=( +( -924 +150*c1y +60*c2y +10*c3y)
            +( +150 +225*c1y +90*c2y +15*c3y)*np.cos(x)
            +(  +60  +90*c1y +36*c2y  +6*c3y)*np.cos(2*x)
            +(  +10  +15*c1y  +6*c2y    +c3y)*np.cos(3*x))/1024
    elif m==4:
        c1y=np.cos(1*y)
        c2y=np.cos(2*y)
        c3y=np.cos(3*y)
        c4y=np.cos(4*y)
        f=(+(-146424*c1y +13676*c2y +4920*c3y -1355*c4y +288543)
           +(-115776*c1y +36512*c2y +3648*c3y -2216*c4y -146424)*np.cos(x)
           +( +36512*c1y +20080*c2y -3488*c3y  -732*c4y  +13676)*np.cos(2*x)
           +(  +3648*c1y  -3488*c2y -2112*c3y  +104*c4y   +4920)*np.cos(3*x)
           +(  -2216*c1y   -732*c2y  +104*c3y   -25*c4y   -1355)*np.cos(4*x))/442368
    elif m==5:
        c1y=np.cos(1*y)
        c2y=np.cos(2*y)
        c3y=np.cos(3*y)
        c4y=np.cos(4*y)
        c5y=np.cos(5*y)
        f=(+(+2077164 -795900*c1y  -55440*c2y  +63930*c3y  -3420*c4y -2238*c5y)
           +( -795900 -905180*c1y   +8400*c2y +123210*c3y  -3540*c4y -9070*c5y)*np.cos(x)
           +( -55440    +8400*c1y +125120*c2y  +60360*c3y  -8240*c4y -7320*c5y)*np.cos(2*x)
           +( +63930  +123210*c1y  +60360*c2y  -12975*c3y -12930*c4y +1125*c5y)*np.cos(3*x)
           +(  -3420    -3540*c1y   -8240*c2y  -12930*c3y  -3700*c4y +1110*c5y)*np.cos(4*x)
           +(  -2238    -9070*c1y   -7320*c2y   +1125*c3y  +1110*c4y  -503*c5y)*np.cos(5*x))/(2*1376256);
    elif m==6:
        c1y=np.cos(1*y)
        c2y=np.cos(2*y)
        c3y=np.cos(3*y)
        c4y=np.cos(4*y)
        c5y=np.cos(5*y)
        c6y=np.cos(6*y)
        f=(+(+34522852 -15494160*c1y   -99390*c2y +1480600*c3y -321540*c4y  -82824*c5y +36190*c6y)
           +(-15494160 -15708480*c1y +2463480*c2y +2181600*c3y -658800*c4y -112800*c5y +49800*c6y)*np.cos(x)
           +(   -99390  +2463480*c1y +3428745*c2y  +333900*c3y -537810*c4y  +13500*c5y +19335*c6y)*np.cos(2*x)
           +( +1480600  +2181600*c1y  +333900*c2y  -640720*c3y -204120*c4y  +77040*c5y  +7540*c6y)*np.cos(3*x)
           +(  -321540   -658800*c1y  -537810*c2y  -204120*c3y  +29220*c4y  +33480*c5y   +690*c6y)*np.cos(4*x)
           +(   -82824   -112800*c1y   +13500*c2y   +77040*c3y  +33480*c4y   -1104*c5y  -1020*c6y)*np.cos(5*x)
           +(   +36190    +49800*c1y   +19335*c2y    +7540*c3y    +690*c4y   -1020*c5y   +105*c6y)*np.cos(6*x))/48234496;

    if k>1:
        f=1-(1-f)**k
        n=np.floor(n/k)
        
    f=f**n
    return(f)

    
