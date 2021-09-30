import numpy as np
import scipy.special as sp
import jutil as ju
import math
from dedalus import public as de
from dedalus import core as co

# Band limited signals everything
# is built around the odd heaviside function and its derivative
# k sets the highest wavenumber of the function
# M sets the order of the highest vanishing derivative
# Heaviside function and derivative
# k is maximum wavenumber, M is order
# x are evaluation points
# H(x) = 1 x>>0
#      = 0 x<<0

def beta_ab(h,w,H):
    c = h*(H-h)/w**2-1
    a = h/H*c
    b = (1-h/H)*c
    return(a,b)

def beta_Heaviside(x,h,w,H):
    (a,b)=beta_ab(h,w,H)
    f = sp.betainc(a,b,x/H)

def beta_Delta(x,h,w,H):
    t = x/H
    d = np.power(t,a)*np.power(1-t,b)/(H*sp.beta(a,b))


class ffta:
    def __init__(self,b,M,nfx):
        self.name = b.name
        self.x0   = b.interval[0]
        self.span = b.interval[1]-b.interval[0]
        x=b.grid(scale=1)
        self.N    = len(x)
        self.dx   = self.span/self.N
        if   isinstance(b,co.basis.Chebyshev): self.maxk=math.pi/self.dx
        elif isinstance(b,co.basis.SinCos):    self.maxk=math.pi*(self.N-1)/self.span
        elif isinstance(b,co.basis.Fourier):   self.maxk=math.pi*(self.N-2)/self.span
        if not isinstance(b,co.basis.Chebyshev): maxk=max(np.abs(b.wavenumbers))
        else: maxk=self.maxk
        #ju.logger.info('ffta.init: {} {} maxk {} dx {}'.format(self.name,type(b),self.maxk/maxk,self.dx/max(np.diff(x))))
        self.width  = (M+1/2)*self.dx
        self.iwidth = M+1
        self.M      = M
        self.nfx    = nfx

        # cos^j(x)*sin^k(x) is power, number of vanising derivatives
        # m is oversampling ratio
        # dx is the grid spacing on the unaliased grid
    def deltajks(self,x,j,k,m):
        w  = math.pi/(self.dx*(j+k)*m) # Wavenumber
        y  = np.abs(w*x)
        c  = np.where(y<=math.pi/2)[0]
        s  = slice(c[0],c[-1]+1)
        y  = y[s]
        f  = np.power(np.cos(y),j)*np.power(np.sin(y),k)
        df = 0
        if j != 0:   df += -w*j*np.power(np.cos(y),j-1)*np.power(np.sin(y),k+1)
        if k != 0:   df += +w*k*np.power(np.cos(y),j+1)*np.power(np.sin(y),k-1)
        a  = (k+1)/2
        b  = (j+1)/2
        N  = sp.beta(a,b)/(2*w)
        fI = sp.betainc(a,b,np.sin(y)**2) 
        return(s,f/N,fI*np.sign(x[s]),N,df*np.sign(x[s])/N)

     # Return a compact Heaviside function that is band limited
        # Includes the end points which are 0 and 1
        # j is order with which derivative vanish
        # m is oversampling ratio
    def Hjk(self,j,m):
        n = np.int64(np.ceil(2+2*j*m))
        t = np.sin(np.linspace(0,math.pi/2,n))**2
        a = j+1/2
        f = sp.betainc(a,a,t)
        return(f)

    def Djk(self,j,m):
        n = np.int64(np.ceil(2+2*j*m))
        f = np.sin(np.linspace(0,math.pi,n))**(2*j)
        return(f)

    def HH(self,j,m,p):
        H=self.Hjk(j,m)
        n=H.shape[0]
        a1=np.int64(np.round(n*p[0]))
        a2=np.int64(np.round(n*p[1]))
        a3=np.int64(np.round(n*p[2]))
        f=np.concatenate((np.zeros(a1,),H,np.ones(a2,),np.flip(H),np.zeros(a3,)))
        #ju.logger.info('n {} p {} a {} {} {} f {}'.format(n,p,a1,a2,a3,f))
        return(f,slice(0,f.shape[0]))

    def heavisidel(self,x,n=1): #Left decreasing heaviside function
        f=self.lf(self.xs(x,n))
        z=np.asarray(np.where(f.flatten()<1))
        sl=slice(1,z.max()+1)
        f=(1-f[sl,...])/2
        return(f,sl)

    def heavisider(self,x,n=1): #right increasing heaviside function
        f=self.lf(self.xs(x,n))
        z=np.asarray(np.where(f.flatten()>-1))
        sl=slice(z.min(),x.shape[0]+1)
        f=(1+f[sl,...])/2
        return(f,sl)

    def df0(self,n): #d=du/dx
        return(self.sf(n)*self.ldf0())

    def df(self,x,n): #d=du/dx
        return(self.sf(n)*self.ldf(self.xs(x,n)))

    def Ndf(self,x,n): # Normalised to maximum value 1
        return(self.ldf(self.xs(x,n))/self.ldf0())

    def dfsl(self,x,n): #return only non zero section and slice
        df=self.sf(n)*self.ldf(self.xs(x,n))
        z=np.asarray(np.where(df.flatten()!=0))
        sl=slice(z.min(),z.max()+1)
        return(df[sl,...],sl)

    def f(self,x,n): #u
        return(self.lf(self.xs(x,n)))

    def fsl(self,x,n): #u
        f  = self.lf(self.xs(x,n))
        z  = np.asarray(np.where(np.abs(f.flatten())!=1))
        sl = slice(z.min(),z.max()+1)
        return(f[sl,...],sl)

    def fdfsl(self,x,n): #return only non zero section and slice
        f=self.lf(self.xs(x,n))
        df=self.sf(n)*self.ldf(self.xs(x,n))
        zf  = np.asarray(np.where(np.abs(f.flatten())!=1))
        zdf = np.asarray(np.where(df.flatten()!=0))
        sl=slice(min(zf.min(),zdf.min()),max(zf.max(),zdf.max())+1)
        return(f[sl,...],df[sl,...],sl)

    
    def Ndfsl(self,x,n): #return only non zero section and slice
        df=self.ldf(self.xs(x,n))/self.ldf0()
        z=np.asarray(np.where(df.flatten()!=0))
        sl=slice(z.min(),z.max()+1)
        return(df[sl,...],sl)


    def dfw(self,x,w): # scaled by width instead of maxk
        sf= math.pi/2/w;
        return(sf*self.ldf(sf*x))

    def Ndfw(self,x,w): # scaled by width instead of maxk
        return(self.ldf(sf*x)/self.ldf0())

    def dfslw(self,x,w): # scaled by width instead of maxk
        sf= math.pi/2/w;
        df=sf*self.ldf(sf*x)
        z=np.asarray(np.where(df.flatten()!=0))
        sl=slice(z.min(),z.max()+1)
        return(df[sl,...],sl)

    def Ndfslw(self,x,w): # scaled by width instead of maxk
        df=self.ldf(sf*x)/self.ldf0()
        z=np.asarray(np.where(df.flatten()!=0))
        sl=slice(z.min(),z.max()+1)
        return(df[sl,...],sl)

    def wx(self,n): # width in units
        x=math.pi/2/self.sf(n)
        return(x)
        
    def wn(self,n): # width in gridpoints
        return(np.round(self.wx(self,n)/self.dx))

    def sf(self,n):
        return(self.maxk/self.M/n)

    # f(0)=f(L)=0
    # f(x)=1 nearly everywhere else
    # f(x+L)=f(x)
    # f is periodic and even
    # Should call with n=2 in a periodic domain
    # d is derivative
    def periodic_constant(self,x,n=1,d=0):
        m=np.int(np.floor(self.N/n))
        sf=2*math.pi/self.span
        x=sf*(x-self.x0)
        c  = 1 - sp.binom(2*m-1, m-1)/2**(2*m-1)
        d4=np.mod(d,4)
        if d==0: f = c
        else:    f = np.zeros_like(x)
        #print('m:{:2d} {:2d}: {:14.12}'.format(m,0,c))
        for n in range(1,m+1):
            if n==1: c  = - sp.binom(2*m,   m-1)/2**(2*m-1)
            else:    c *= (1+m-n)/(m+n)
            #print('{:2d}: {:14.12e}'.format(n,c))
            a = c*(sf*n)**d
            if   d4==0: f += a*np.cos(n*x)
            elif d4==1: f -= a*np.sin(n*x)
            elif d4==2: f -= a*np.cos(n*x)
            elif d4==3: f += a*np.sin(n*x)
        #print('type(f) {} f.shape {} type(d) {} d={}'.format(type(f),f.shape,type(d),d))
        return(f)
    
    # Largest wavenumber is (2*m-1)
    # odd constant
    def odd_constant(self,x,n=1,d=0):
        m=np.int(np.floor(self.N/n/2))
        sf=math.pi/self.span
        x=sf*(x-self.x0)
        f  = np.zeros_like(x)
        c   = np.float64(-1)
        d4=np.mod(d,4)
        for n in range(m-1): c*=(2*n+3)**2/(4*(n+1)*(n+2))
        for n in range(m):
            c *= (2*n-1)*(m-n)/((m+n)*(2*n+1))
            a = c*(sf*(2*n+1))**d
            if   d4==0: f += a*np.sin((2*n+1)*x)
            elif d4==1: f += a*np.cos((2*n+1)*x)
            elif d4==2: f -= a*np.sin((2*n+1)*x)
            elif d4==3: f -= a*np.cos((2*n+1)*x)
        #print('type(f) {} f.shape {} type(d) {} d={}'.format(type(f),f.shape,type(d),d))
        return(f)
    
    # odd constant with vanishing 3rd derivative at the origin
    def odd_constantd3(self,x,n=1,d=0):
        m=np.int(np.floor(self.N/n/2))
        sf=math.pi/self.span
        x=sf*(x-self.x0)
        c   = np.float64(27/28)
        d4=np.mod(d,4)
        for n in range(1,m-1): c*=(6*n+1)*(3+2*n)**2/(4*n*(n+2)*(6*n+7))
        a = c*sf**d
        if   d4==0: f = +a*np.sin(x)
        elif d4==1: f = +a*np.cos(x)
        elif d4==2: f = -a*np.sin(x)
        elif d4==3: f = -a*np.cos(x)    
        for n in range(1,m):
            if (2*n**2-2*n+3-3*m)==0: c=co*(2*n-3)*(m-n+1)*(-2*n**2+3*m-2*n-3)*(m-n)/((-2*n**2+3*m+6*n-7)*(m+n-1)*(2*n+1)*(m+n))
            else:
                co=c
                c *= (m-n)/(m+n)*(2*n-1)/(2*n+1)*(2*n**2+2*n+3-3*m)/(2*n**2-2*n+3-3*m)
            a = c*(sf*(2*n+1))**d
            if   d4==0: f += a*np.sin((2*n+1)*x)
            elif d4==1: f += a*np.cos((2*n+1)*x)
            elif d4==2: f -= a*np.sin((2*n+1)*x)
            elif d4==3: f -= a*np.cos((2*n+1)*x)    
        #print('type(f) {} f.shape {} type(d) {} d={}'.format(type(f),f.shape,type(d),d))
        return(f)

    def constant(self,x,n=1):
        f=2*self.heaviside(x,n)[0]-2*self.heaviside(x-self.span,n)[0]-1
        return(f)

    # even constant
    def even_constant(self,x,n=1,d=0):
        m=np.int(np.floor(self.N/n))-1
        sf=math.pi/self.span
        x=sf*(x-self.x0)
        f  = np.zeros_like(x)
        c   = np.float64(-1)
        d4=np.mod(d,4)
        for n in range(0,m+1):
            if n==0:   c = 1 - sp.binom(2*m-1, m-1)/2**(2*m-1)
            elif n==1: c =    -sp.binom(2*m,   m-1)/2**(2*m-1)
            else:      c *= (1+m-n)/(m+n)
            a = c*(sf*n)**d
            if   d4==0: f += a*np.cos(n*x)
            elif d4==1: f -= a*np.sin(n*x)
            elif d4==2: f -= a*np.cos(n*x)
            elif d4==3: f  = a*np.sin(n*x)
            #print('m {} n {} c {}'.format(m,n,c))
        return(f)

    def heaviside(self,x,n=1):
        M=self.nfx
        xx=x;
        k=self.maxk/n/(2*M+1)
        x  = np.minimum(math.pi/2,np.maximum(-math.pi/2,k*x))
        If = np.maximum(0,(x+math.pi)/(2*k))+np.maximum(0,xx-math.pi/(2*k))
        f  = 1/2+np.zeros_like(x)
        If0=0
        df = np.zeros_like(x)
        for n in range(M+1):
            j    = 2*n+1
            sf   = j*k
            c    = 2/(j*math.pi)*np.exp(2*sp.gammaln(M+3/2)-sp.gammaln(M+2+n)-sp.gammaln(M+1-n))
            If  += -c/sf*np.cos(j*x)
            f   +=  c   *np.sin(j*x)
            df  +=  c*sf*np.cos(j*x)
            df0 = k*np.exp(sp.gammaln(2*M+2)-M*np.log(4)-2*sp.gammaln(M+1))
        return(f,df,df0,If)
    
    def xs(self,x,n):
        #x=np.minimum(math.pi/2,np.maximum(-math.pi/2,x*self.sf(n)))
        x=x*self.sf(n)
        return(x)


    # These functions have
    # at x=0:    f=0, f''=0, f'''=0, f''''=0
    # at x=pi/2: f=1, f'=0 ... f^(n)=0 for n up to M-2 
    def lf(self,x):
        if   self.M==5:  y =(        450*np.sin(x)+    25*np.sin(3*x)-    9*np.sin(5*x))/416
        elif self.M==7:  y =(      11025*np.sin(x)+  1225*np.sin(3*x)-  147*np.sin(5*x)-   75*np.sin(7*x))/9728
        elif self.M==9:  y =(      23814*np.sin(x)+  3528*np.sin(3*x)                  -  243*np.sin(7*x)-   49*np.sin(9*x))/20480
        elif self.M==11: y =(   4802490*np.sin(x)+838530*np.sin(3*x)+68607*np.sin(5*x)-49005*np.sin(7*x)-21175*np.sin(9*x)-2835*np.sin(11*x))/4063232
        elif self.M==13: y =(  46378332*np.sin(x)+9018009*np.sin(3*x)+1288287*np.sin(5*x)-368082*np.sin(7*x)-286286*np.sin(9*x)-74529*np.sin(11*x)-7623*np.sin(13*x))/38797312
        elif self.M==15: y =( 869593725*np.sin(x)+182507325*np.sin(3*x)+34783749*np.sin(5*x)-3764475*np.sin(7*x)-6181175*np.sin(9*x)-2395575*np.sin(11*x)-467775*np.sin(13*x)-39039*np.sin(15*x))/721420288
        elif self.M==17: y =(7978177350*np.sin(x)+1772928300*np.sin(3*x)+406161756*np.sin(5*x)-57857800*np.sin(9*x)-30431700*np.sin(11*x)-8583300*np.sin(13*x)-1363791*np.sin(15*x)-96525*np.sin(17*x))/6576668672;
        y[x >=  math.pi/2]  = 1
        y[x <= -math.pi/2] = -1
        
        return(y)

    def ldf(self,x):
        if   self.M==5:  y =   15*(  30*np.cos(x)+         5*np.cos(3*x)-         3*np.cos(5*x))/416
        elif self.M==7:  y =  105*( 105*np.cos(x)+        35*np.cos(3*x)-         7*np.cos(5*x)       -5*np.cos(7*x))/9728
        elif self.M==9:  y =   63*( 378*np.cos(x)+       168*np.cos(3*x)                             -27*np.cos(7*x)        -7*np.cos(9*x))/20480
        elif self.M==11: y = 3465*(1386*np.cos(x)+       726*np.cos(3*x)+        99*np.cos(5*x)-      99*np.cos(7*x)       -55*np.cos(9*x)        -9*np.cos(11*x))/4063232
        elif self.M==13: y =(  46378332*np.cos(x)+  27054027*np.cos(3*x)+   6441435*np.cos(5*x)- 2576574*np.cos(7*x)  -2576574*np.cos(9*x)   -819819*np.cos(11*x)    -99099*np.cos(13*x))/38797312
        elif self.M==15: y =( 869593725*np.cos(x)+ 547521975*np.cos(3*x)+ 173918745*np.cos(5*x)-26351325*np.cos(7*x) -55630575*np.cos(9*x) -26351325*np.cos(11*x)  -6081075*np.cos(13*x)  -585585*np.cos(15*x))/721420288;
        elif self.M==17: y =(7978177350*np.cos(x)+5318784900*np.cos(3*x)+2030808780*np.cos(5*x)                     -520720200*np.cos(9*x)-334748700*np.cos(11*x)-111582900*np.cos(13*x)-20456865*np.cos(15*x)-1640925*np.cos(17*x))/6576668672;
        y[x>=math.pi/2]  = 0
        y[x<=-math.pi/2] = 0
        return(y)


    def ldf0(self):
        if   self.M==5:  y =    15/13       
        elif self.M==7:  y =   105/76      
        elif self.M==9:  y =    63/40       
        elif self.M==11: y =  3465/1984   
        elif self.M==13: y =  9009/4736   
        elif self.M==15: y = 45045/22016 
        elif self.M==17: y = 109395/50176
        return(y)

    # Compact delta function away from a the boundaries
    # int(delta(x),x=-infinity..infinity)
    # at x=pi/2: f=0, f'=0 ... f^(n)=0 for n up to M-1 
    def ldelta(self,x):
        M=self.M
        N=math.pi*sp.binom(M,M/2)/2**M
        y=np.where(np.abs(x)<math.pi/2,np.cos(x)**M/N,0)
        return(y)

    def delta(self,x,n):
        s=self.sf(n)
        return(s*self.ldelta(s*x))

    def deltasl(self,x,n):
        s=self.sf(n)
        f=s*self.ldelta(s*x)
        z=np.asarray(np.where(f.flatten()>0))
        sl=slice(z.min(),z.max()+1)
        return(f[sl,:],sl)
       

    # int(sin(y)^(2*n))
    # Heaviside with end points
    # f(x) = 0 x < x1
    #      = 1 x > x2
    # n derivatives vanish at x1 and x2
    
    def Hn(self,x,x1,x2,m):
        n = min(10,max(0,np.int(np.round(abs(x1-x2)/(2*self.dx*m)))))
        Pi=math.pi
        y=np.minimum(Pi,np.maximum(0,(x-x1)*Pi/(x2-x1)))
        f=y/Pi
        for j in range(1,n+1):
            if j==1 : c = -n/(n+1)/Pi
            else:     c *= (j-1)*(j-1-n)/(j*(n+j))
            f += c*np.sin(2*j*y)
        return(f)
    
