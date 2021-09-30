import numpy as np
import scipy.special as sp
import scipy.optimize as op
import math
import os

SERF_DEBUG_SOLVE = os.getenv('SERF_DEBUG_SOLVE') is not None
SERF_DEBUG_FIT   = os.getenv('SERF_DEBUG_FIT') is not None

# Symmetric error function class
# f = A + B*(erf((h-x)/w)+erf((h+x)/w))
# x   not none default coordinates
# W   not none default integration weights
# f0  not None f(0) is fixed
# f1  not None f(1) is fixed
# I   not None int(f(x) dx) is fixed
# A0  [false] A is zero
# h   not none h is fixed
# w   not none w is fixed
#
# There are different forms that are convenient 

Pi=math.pi
spi  = np.sqrt(Pi)
s2pi = np.sqrt(2*Pi)
eps  = np.finfo(float).eps

serf_debug=True

def estimate_yhw(W,y):
    miny=min(y)
    maxy=max(y)
    ry=maxy-miny
    h=(W*(y-miny)).sum()/ry
    w=s2pi*(W*(y-miny)*(maxy-y)).sum()/ry**2
    return(miny,maxy,h,w)

def estimate_yhw1(W,y):
    maxy=max(y)
    h=(W*y).sum()/maxy
    w=s2pi*(W*y*(maxy-y)).sum()/maxy**2
    return(maxy,h,w)

def solve_secant(f,y,H): # find x that is closest to u(x) = y
    x0=0
    x1=H
    y0  = f(0)
    y1  = f(H)
    
    toly=1e-10*abs(y1-y0)
    tolx=1e-10*H
    
    # Consider the case when there is no exact solution
    if (y0-y)*(y1-y)>0:
        if abs(y0-y)<abs(y1-y0): return(0)
        else:                    return(H)
            
    maxit=100
        
    for j in range(maxit):
        # (x1-x0)*z=(x1-x)*y0+(x-x0)*y1 linear interpolation
        # (x1-x0)*z=+(x1*y1-x0*y0) + x*(y1-y0) linear interpolation
        x2 = (x0*(y1 - y)+ x1*(y-y0))/(y1 - y0)
        y2 = f(x2)
        if SERF_DEBUG_SOLVE:
            print('j={:3d} x=[{:8.4f} {:8.4f} {:8.4f} ] y=[{:8.4f},{:8.4f},{:8.4f} ] tol=[{:5.1e}{:5.1e}]'\
                  .format(j,x0,x2,x1,y0,y2,y1,x1-x0,abs(y-y2))) 
        if abs( y-y2)<toly  : break
        if abs(x1-x0)<tolx : break
        if (y0-y)*(y2-y)<=0:
            y1=y2
            x1=x2
        else:
            y0=y2
            x0=x2
    return(x2)
    

def serfutil(x,h,w):
    if w==0:
        xn=np.sign(x-h)*np.inf
        xp=np.sign(x+h)*np.inf
        xnf=np.sign(x-h)
        xpf=np.sign(x+h)
        xnp=0*x
        xpp=0*x
    else:
        xn=(x-h)/w
        xp=(x+h)/w
        xnf=sp.erf(xn)
        xpf=sp.erf(xp)
        xnp=np.exp(-xn**2)/(w*spi)
        xpp=np.exp(-xp**2)/(w*spi)
        
    return(xp,xn,xpf,xnf,xnp,xpp)


# basic serf function 
def serff(x,h,w):
    (xp,xn,xpf,xnf,xnp,xpp)=serfutil(x,h,w)
    f=(xpf-xnf)/2
    return(f);

# Integral of basic serf function 
def serfF(x,h,w):
    if w==0:
        F = np.maximum(-h,np.minimum(h,x))
    else:
        (xp,xn,xpf,xnf,xnp,xpp)=serfutil(x,h,w)
        F= w/2*((xpp - xnp)*w - xn*xnf + xpf*xp)
    return(F)

# Derivative of basic serf function 
def serfd(x,h,w):
    (xp,xn,xpf,xnf,xnp,xpp)=serfutil(x,h,w)
    df=xpp - xnp
    return(df)

# basic serf function with derivatives
def serffd(x,h,w):
    (xp,xn,xpf,xnf,xnp,xpp)=serfutil(x,h,w)
    f   = (xpf-xnf)/2
    if w==0:
        fh  = 0*x
        fw  = fh
        fhh = fh
        fhw = fh
        fww = fh
    else:
        fh  = xpp + xnp
        fw  = xn*xnp - xp*xpp
        fhh = 2/w*(xn*xnp - xp*xpp)
        fhw = 2/w*(xn**2*xnp + xp**2*xpp - xnp - xpp)
        fww = 2/w*(xnp*xn*(xn**2 - 1) + xpp*xp*(1-xp**2))
    return(f,fh,fw,fhh,fhw,fww)

# Integral of basic serf function with derivatives
def serfFd(x,h,w):
    (xp,xn,xpf,xnf,xnp,xpp)=serfutil(x,h,w)
    if w==0:
        F = np.maximum(-h,np.minimum(h,x))
    else:
        F   = w/2*((xpp - xnp)*w - xn*xnf + xpf*xp)
    Fw  = w/2*(xpp - xnp)
    Fh  = (xnf+ xpf)/2
    Fww = xp**2*xpp  - xn**2*xnp
    Fhh = xpp - xnp
    Fhw = -xn*xnp - xp*xpp   
    return(F,Fh,Fw,Fhh,Fhw,Fww)

# serf function normalised so that f(0)=1
def serf1(x,h,w):
    f=serff(x,h,w)/serff(0,h,w)
    return(f)

def serf1d(x,h,w):
    (e,eh,ew,ehh,ehw,eww)=serffd(x,h,w)
    (g,gh,gw,ghh,ghw,gww)=serffd(0,h,w)
    f=e/g
    fh=eh/g-e*gh/g**2
    fw=ew/g-e*gw/g**2
    return(f,fh,fw)
    
def Iserf1(x,h,w):
    f=serfF(x,h,w)/serff(0,h,w)
    return(f)

def dserf1(x,h,w):
    f=serfd(x,h,w)/serff(0,h,w)
    return(f)

# serf function normalised so that int( f(x) dx) = 1
def serfH(x,h,w,H):
    f=serff(x,h,w)/serfF(H,h,w)
    return(f)

# serf function normalised so that int( f(x) dx) = 0 and f(0)=1
def serfzd(x,h,w,H):
    (a,ah,aw,ahh,ahw,aww)=serffd(x,h,w)
    (b,bh,bw,bhh,bhw,bww)=serffd(0,h,w)
    (c,ch,cw,chh,chw,cww)=serffd(H,h,w)
    (d,dh,dw,dhh,dhw,dww)=serfFd(H,h,w)

    ad=a-d/H
    bc=b-c/H
    
    f=ad/bc
    fh=(ah-dh/H)/bc-ad*(bh-ch/H)/bc**2
    fw=(aw-dw/H)/bc-ad*(bw-cw/H)/bc**2
    return(f,fh,fw)

# serf function normalised so that int( f(x) dx) = 0 and f(0)=1
def serfz(x,h,w,H):
    f=(serff(x,h,w)-serfF(H,h,w)/H)/(serff(0,h,w)-serfF(H,h,w)/H)
    return(f)

# Integral of serfz
def Iserfz(x,h,w,H):
    f=(serfF(x,h,w)-x*serfF(H,h,w)/H)/(serff(0,h,w)-serfF(H,h,w)/H)
    return(f)

# Derivative of serfz
def dserfz(x,h,w,H):
    f=serfd(x,h,w)/(serff(0,h,w)-serfF(H,h,w)/H)
    return(f)

#lsqlin best fit for y=a*f+b*g+c*h
def lsqlin3(y,W=1,f=1,g=1,h=1):
    ff = (W*f*f).sum()
    fg = (W*f*g).sum()
    fh = (W*f*h).sum()
    gg = (W*g*g).sum()
    gh = (W*g*h).sum()
    hh = (W*h*h).sum()
    yf = (W*y*f).sum()
    yg = (W*y*g).sum()
    yh = (W*y*h).sum()
    d = ff*gh**2 + hh*fg**2 + gg*fh**2-ff*gg*hh - 2*fg*fh*gh
    a = (-fg*gh*hy + fg*gy*hh + fh*gg*hy - fh*gh*gy - fy*gg*hh + fy*gh**2)/d
    b = (+ff*gh*hy - ff*gy*hh - fg*fh*hy + fg*fy*hh + fh**2*gy - fh*fy*gh)/d
    c = (-ff*gg*hy + ff*gh*gy + fg**2*hy - fg*fh*gy - fg*fy*gh + fh*fy*gg)/d
    y = a*f+b*g+c*h
    return(a,b,c,y)

#lsqlin best fit for y=a*f+b*g+c*h with constraint a*F+b*G+c*H=C
def lsqlin3c(y,W=1,f=1,g=1,h=1,F=1,G=1,H=1,C=1):
    ff = (W*f*f).sum()
    fg = (W*f*g).sum()
    fh = (W*f*h).sum()
    gg = (W*g*g).sum()
    gh = (W*g*h).sum()
    hh = (W*h*h).sum()
    yf = (W*y*f).sum()
    yg = (W*y*g).sum()
    yh = (W*y*h).sum()
    d = + F**2*(gh**2-gg*hh) + G**2*(fh**2-ff*hh) + H**2*(fg**2-ff*gg) + 2*H*G(ff*gh - fg*fh) + 2*F*G*(hh*fg - fh*gh)  + 2*F*H*(gg*fh-fg*gh) 
    a = ((+fh*hy - fy*hh)*G**2 + ((2*fy*gh - fg*hy - fh*gy)*H + (-C*fh - F*hy)*gh + hh*(C*fg + F*gy))*G + (+fg*gy - fy*gg)*H**2 + ((-C*fg - F*gy)*gh + gg*(C*fh + F*hy))*H - C*F*(gg*hh - gh**2))/d
    b = ((+gh*hy - gy*hh)*F**2 + ((2*fh*gy - fg*hy - fy*gh)*H + (-C*gh - G*hy)*fh + hh*(C*fg + G*fy))*F + (-ff*gy + fg*fy)*H**2 + ((-C*fg - G*fy)*fh + ff*(C*gh + G*hy))*H - C*G*(ff*hh - fh**2))/d
    c = ((-gg*hy + gh*gy)*F**2 + ((2*fg*hy - fh*gy - fy*gh)*G + (-C*gh - H*gy)*fg + gg*(C*fh + H*fy))*F + (-ff*hy + fh*fy)*G**2 + ((-C*fh - H*fy)*fg + ff*(C*gh + H*gy))*G - C*H*(ff*gg - fg**2))/d
    y = a*f+b*g+c*h
    return(a,b,c,y)

# lsqlin best fit for y=a*f+b*g
def lsqlin2(y,W=1,f=1,g=1):
    ff = (W*f*f).sum()
    fg = (W*f*g).sum()
    gg = (W*g*g).sum()
    yf = (W*y*f).sum()
    yg = (W*y*g).sum()
    d = ff*gg - fg**2
    a = (gg*yf-fg*yg)/d
    b = (ff*yg-fg*yf)/d
    y = a*f+b*g
    return(a,b,y)


# lsqlin best fit for y=a*f+b*g subject to linear constraint a*F+b*G=c
def lsqlin2c(y,W=1,f=1,g=1,F=1,G=1,C=1):
    ff = (W*f*f).sum()
    fg = (W*f*g).sum()
    gg = (W*g*g).sum()
    yf = (W*y*f).sum()
    yg = (W*y*g).sum()
    d = F**2*gg - 2*F*G*fg + G**2*ff
    a = (G**2*yf - F*G*yg + F*gg*C - G*fg*C)/d
    b = (F**2*yg - F*G*yf - F*fg*C + G*ff*C)/d
    y = a*f+b*g
    return(a,b,y)

def fit_serf(x,y,W,h=None,w=None):
    if h is None or w is None:
        (miny,maxy,eh,ew)=estimate_yhw(W,y)
        if h is None: h=eh
        if w is None: w=ew
         
    H  = W.sum()
    WS=np.sqrt(W)
    def f(p): return(WS*(y-lsqlin2(y,W,serff(x,p[0],p[1]))[2]))
    p0=np.array([h,w])
    lp=np.array([0,0])
    up=np.array([H,H])
    r=op.least_squares(f, np.array([h,w]), bounds=(lp,up), method='dogbox')
    h=r.x[0]
    w=r.x[1]
    (a,b,yy)=lsqlin2(y,W,serff(x,h,w))
    if SERF_DEBUG_FIT:
        y0=b+a*serff(0,h,w)
        y1=b+a*serff(H,h,w)
        minyy=min(yy)
        maxyy=max(yy)
        print('h {:6.4f} w {:6.4f} a {:6.4f} b {:6.4f} y0 {:6.4f} y1 {:6.4f} min {:6.4f} max {:6.4f}'.format(h,w,a,b,y0,y1,minyy,maxyy))
        r1=(((WS*(yy-y))**2).sum())/2
        r2=((r.fun**2).sum())/2
        print('cost {:8.6f} r1 {:8.6f} r2 {:8.6f}'.format(r.cost,r1,r2))
    #print('h {:} w {:} a {:} b {:} yy{:}'.format(h,w,a,b,yy))
    return(h,w,a,b,yy)
    
# the general form is A + C*serf(x,h,w)
#
# density  form
# b    = B*serf1(x,h,w)
# b(0) = B
# b(x) -> 0 as x-> infinity
#
# velocity form is 
# u    = (1+U)*serfz(x,h,w,H)-1
# u(0) = U
# int u = -H

class fb:
    def __init__(self,x,W,H,A=None,h=None,w=None):
        self.x=x
        self.W=W
        self.H=H
        self.fh=h # fixed value for h     
        self.fw=w # fixed value for w  
        self.fA=A # fixed value for at x=0
        self.h=h  # current best fit h
        self.w=w  # current best fit w
        self.A=A  # current best fit U
        self.dx=np.diff(x).min()   # smallest dx

    def print(self):
        print('class fb')
        if self.x  is not None : print('  x  [{:6.2f} {:6.2f}]'.format(min(self.x),max(self.x)))
        if self.H  is not None : print('  H  {:6.2f}'.format(self.H))
        if self.fh is not None : print('  fh {:6.2f}'.format(self.fh))
        if self.fw is not None : print('  fw {:6.2f}'.format(self.fw))
        if self.fA is not None : print('  fA {:6.2f}'.format(self.fA))
        if self.h  is not None : print('  h  {:6.2f}'.format(self.h))
        if self.w  is not None : print('  w  {:6.2f}'.format(self.w))
        if self.A  is not None : print('  A  {:6.2f}'.format(self.A))
        if self.dx is not None : print('  dx {:6.2e}'.format(self.dx))
        
    
    def f(self,x): # evaluate the density function
        if self.A is None or self.h is None or self.w is None: f=0*x
        else: f=self.A*serf1(x,self.h,self.w)
        return(f)
        
    def If(self,x): # evaluate the density function
        f=self.A*Iserf1(x,self.h,self.w)
        return(f)
        
    def df(self,x): # evaluate the density function
        f=self.A*dserf1(x,self.h,self.w)
        return(f)

    def fit_Ahw(self,y): # Fit the density function
        def f(x,A,h,w): return(A*serf1(x,h,w))
        def J(x,A,h,w):  (g,gh,gw)=serf1d(x,h,w) ; return(np.transpose([g,A*gh,A*gw]))
        (self.A,self.h,self.w)=op.curve_fit(f,self.x, y,[self.A,self.h,self.w], check_finite=True, method='lm', jac=J)[0]

    def fit_Ah(self,y): # Fit the density function
        def f(x,A,h): return(A*serf1(x,h,self.fw))
        def J(x,A,h):  (g,gh,gw)=serf1d(x,h,self.fw) ; return(np.transpose([g,A*gh]))
        (self.A,self.h)=op.curve_fit(f,self.x, y,[self.A,self.h], check_finite=True, method='lm', jac=J)[0]
 
    def fit_Aw(self,y): # Fit the density function
        def f(x,A,w): return(A*serf1(x,self.fh,w))
        def J(x,A,w):  (g,gh,gw)=serf1d(x,self.fh,w) ; return(np.transpose([g,A*gw]))
        (self.A,self.w)=op.curve_fit(f,self.x, y,[self.A,self.w], check_finite=True, method='lm', jac=J)[0]
 
    def fit_hw(self,y): # Fit the density function
        def f(x,h,w): return(self.fA*serf1(x,h,w))
        def J(x,h,w):  (g,gh,gw)=serf1d(x,h,w) ; return(np.transpose([self.fA*gh,self.fA*gw]))
        (self.h,self.w)=op.curve_fit(f,self.x, y,[self.h,self.w], check_finite=True, method='lm', jac=J)[0]
        
    def fit_A(self,y): # Fit the density function
        def f(x,A): return(A*serf1(x,self.fh,self.fw))
        def J(x,A):  (g,gh,gw)=serf1d(x,self.fh,self.fw) ; return(np.transpose([g]))
        self.A=op.curve_fit(f,self.x, y,[self.A], check_finite=True, method='lm', jac=J)[0][0]

    def fit_w(self,y): # Fit the density function
        def f(x,w): return(self.fA*serf1(x,self.fh,w))
        def J(x,w):  (g,gh,gw)=serf1d(x,self.fh,w) ; return(np.transpose([self.fA*gw]))
        self.w=op.curve_fit(f,self.x, y,[self.w], check_finite=True, method='lm', jac=J)[0][0]

    def fit_h(self,y): # Fit the density function
        def f(x,h): return(self.fA*serf1(x,h,self.fw))
        def J(x,h):  (g,gh,gw)=serf1d(x,h,self.fw) ; return(np.transpose([self.fA*gh]))
        self.h=op.curve_fit(f,self.x, y,[self.h], check_finite=True, method='lm', jac=J)[0][0]
       
    def fit(self,y): # Fit the density function
          
        (maxy,h,w)=estimate_yhw1(self.W,y)
        if self.A is None: self.A = maxy
        if self.h is None: self.h = h
        if self.w is None: self.w = w

        self.A = maxy
        self.h = h   
        self.w = w   
        
        if SERF_DEBUG_FIT: A=self.A;h=self.h;w=self.w;

        while True:
            if   self.fA is     None and self.fh is     None and self.fw is     None: self.fit_Ahw(y)
            elif self.fA is     None and self.fh is not None and self.fw is     None: self.fit_Aw(y)
            elif self.fA is     None and self.fh is     None and self.fw is not None: self.fit_Ah(y)
            elif self.fA is     None and self.fh is not None and self.fw is not None: self.fit_A(y)
            elif self.fA is not None and self.fh is     None and self.fw is     None: self.fit_hw(y)
            elif self.fA is not None and self.fh is not None and self.fw is     None: self.fit_w(y)
            elif self.fA is not None and self.fh is     None and self.fw is not None: self.fit_h(y)
            elif self.fA is not None and self.fh is not None and self.fw is not None: print('Nothing to fit')
            if self.w>=0: break
            self.w=-self.w
            
        if SERF_DEBUG_FIT:
            print('        {:8s} {:8s} {:8s}'.format('A','h','w'))
            print('Initial {:8.4f} {:8.4f} {:8.4f}'.format(A,h,w))
            print('Final   {:8.4f} {:8.4f} {:8.4f}'.format(self.A,self.h,self.w))
            print('Change  {:8.4f} {:8.4f} {:8.4f}'.format(self.A-A,self.h-h,self.w-w))
 
    def solve(self,b):
        return(solve_secant(self.f,b,self.H)) # find x that is closest to b(x) = y


class fu:
    def __init__(self,x,W,H,A=None,h=None,w=None):
        self.x=x
        self.W=W
        self.H=H    
        self.fh=h # fixed value for h     
        self.fw=w # fixed value for w  
        self.fA=A # fixed value for at x=0
        self.h=h  # current best fit h
        self.w=w  # current best fit w
        self.A=A  # current best fit A
        self.dx=np.diff(x).min()   # smallest dx

    def print(self):
        print('class fu');
        if self.x  is not None : print('  x  [{:6.2f} {:6.2f}]'.format(min(self.x),max(self.x)))
        if self.H  is not None : print('  H  {:6.2f}'.format(self.H))
        if self.fh is not None : print('  fh {:6.2f}'.format(self.fh))
        if self.fw is not None : print('  fw {:6.2f}'.format(self.fw))
        if self.fA is not None : print('  fA {:6.2f}'.format(self.fA))
        if self.h  is not None : print('  h  {:6.2f}'.format(self.h))
        if self.w  is not None : print('  w  {:6.2f}'.format(self.w))
        if self.A  is not None : print('  A  {:6.2f}'.format(self.A))
        if self.dx is not None : print('  dx {:6.2e}'.format(self.dx))
        
    def fAhw(self,x,A,h,w): # evaluate the velocity function
        f=(1+A)*serfz(x,h,w,self.H)-1
        return(f)

    def fAhwd(self,x,A,h,w): # evaluate the velocity function
        (g,gh,gw)=serfzd(x,h,w,self.H)
        #f=(1+A)*g-1
        return(np.transpose([g,(1+A)*gh,(1+A)*gw]))

    def f(self,x=None): # evaluate the velocity function
        if self.A is None or self.h is None or self.w is None: f=0*x-1
        elif x is not None: f=(1+self.A)*serfz(x,self.h,self.w,self.H)-1
        else:          f=(1+self.A)*serfz(self.x,self.h,self.w,self.H)-1
        return(f)
        
    def If(self,x): # evaluate the velocity function
        f=(1+self.A)*Iserfz(x,self.h,self.w,self.H)-x
        return(f)
        
    def df(self,x): # evaluate the velocity function
        f=(1+self.A)*dserfz(x,self.h,self.w,self.H)
        return(f)

    def fit_Ahw(self,y): # Fit the velocity function
        def f(x,A,h,w): return(self.fAhw(x,A,h,w))
        def J(x,A,h,w): return(self.fAhwd(x,A,h,w))
        (self.A,self.h,self.w)=op.curve_fit(f,self.x, y,[self.A,self.h,self.w], check_finite=False, method='lm', jac=J)[0]
 
    def fit_Ah(self,y): # Fit the velocity function
        def f(x,A,h): return(self.fAhw( x,A,h,self.fw))
        def J(x,A,h): return(self.fAhwd(x,A,h,self.fw)[:,(0,1)]) 
        (self.A,self.h)=op.curve_fit(f,self.x, y,[self.A,self.h], check_finite=False, method='lm', jac=J)[0]
 
    def fit_Aw(self,y): # Fit the velocity function
        def f(x,A,w): return(self.fAhw( x,A,self.fh,w))
        def J(x,A,w): return(self.fAhwd(x,A,self.fh,w)[:,(0,2)]) 
        (self.A,self.w)=op.curve_fit(f,self.x, y,[self.A,self.w], check_finite=False, method='lm', jac=J)[0]
 
    def fit_hw(self,y): # Fit the velocity function
        def f(x,h,w): return(self.fAhw( x,self.fA,h,w))
        def J(x,h,w): return(self.fAhwd(x,self.fA,h,w)[:,(1,2)]) 
        (self.h,self.w)=op.curve_fit(f,self.x, y,[self.h,self.w], check_finite=False, method='lm', jac=J)[0]

    def fit_A(self,y): # Fit the velocity function
        def f(x,A): return(self.fAhw( x,A,self.fh,self.fw))
        def J(x,A): return(self.fAhwd(x,A,self.fh,self.w)[:,(0,)]) 
        self.A=op.curve_fit(f,self.x, y,[self.A], check_finite=False, method='lm', jac=J)[0][0]

    def fit_w(self,y): # Fit the velocity function
        def f(x,w): return(self.fAhw( x,self.fA,self.fh,w))
        def J(x,w): return(self.fAhwd(x,self.fA,self.fh,w)[:,(2,)]) 
        self.w=op.curve_fit(f,self.x, y,[self.w], check_finite=False, method='lm', jac=J)[0][0]

    def fit_h(self,y): # Fit the velocity function
        def f(x,h): return(self.fAhw( x,self.fA,h,self.fw))
        def J(x,h): return(self.fAhwd(x,self.fA,h,self.fw)[:,(1,)]) 
        self.h=op.curve_fit(f,self.x, y,[self.h], check_finite=False, method='lm', jac=J)[0][0]

    def fit(self,y): # Fit the velocity function

        (miny,maxy,h,w)=estimate_yhw(self.W,y)
        if self.A is None: self.A = maxy
        if self.h is None: self.h = h
        if self.w is None: self.w = w
        
        my  = (self.W*y).sum()/self.H
        y   = y-my   # subtract the mean value

 
  
        if SERF_DEBUG_FIT: A=self.A;h=self.h;w=self.w;

        if   self.fA is     None and self.fh is     None and self.fw is     None: self.fit_Ahw(y)
        elif self.fA is     None and self.fh is not None and self.fw is     None: self.fit_Aw(y)
        elif self.fA is     None and self.fh is     None and self.fw is not None: self.fit_Ah(y)
        elif self.fA is     None and self.fh is not None and self.fw is not None: self.fit_A(y)
        elif self.fA is not None and self.fh is     None and self.fw is     None: self.fit_hw(y)
        elif self.fA is not None and self.fh is not None and self.fw is     None: self.fit_w(y)
        elif self.fA is not None and self.fh is     None and self.fw is not None: self.fit_h(y)
        elif self.fA is not None and self.fh is not None and self.fw is not None: print('Nothing to fit')
        
        if SERF_DEBUG_FIT:
            print('        {:8s} {:8s} {:8s}'.format('A','h','w'))
            print('Initial {:8.4f} {:8.4f} {:8.4f}'.format(A,h,w))
            print('Final   {:8.4f} {:8.4f} {:8.4f}'.format(self.A,self.h,self.w))
            print('Change  {:8.4f} {:8.4f} {:8.4f}'.format(self.A-A,self.h-h,self.w-w))
       
    def solve(self,u):
        return(solve_secant(self.f,u,self.H)) # find x that is closest to u(x) = y
