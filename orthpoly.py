import numpy as np

# Create
# These are scaled Legendre polynomials
# use np.polynomial.legendre.legfit(x, y, deg[, rcond, full, w])
#op.opoly_test(3,10)
def opoly_test(m,N):
    L=3
    x=np.linspace(0,L,N)
    s=opoly(x,m,L)
    print('m={}, len(f)={}'.format(m,len(s.f)))
    print('sum(w)={}'.format(s.w.sum()))
    print('    ', end='')
    for k in range(0,m+1): print('{:8d}'.format(k), end='')
    print('')
    for j in range(0,m+1):
        print('{:2d}  '.format(j), end='')
        for k in range(0,m+1):
    	    print('{:8.4f}'.format((s.f[j]*s.f[k]*s.w).sum()), end='')
        print('')

    for j in range(0,m+1):
        c=s.f2c(s.f[j])
        print('j={:2d}, c={}'.format(j,c))
        
    
class opoly:
    def __init__(self,x,m,L=None,w=None):
        if not L: L=1
        self.x=x
        self.L=L
        self.f=list()
        self.f.append(np.ones_like(x)/np.sqrt(L))
        if m>0: self.f.append(np.sqrt(3/L)*(2*x/L-1))
        for n in np.arange(2,m+1): self.f.append(np.sqrt((2*n+1))*(np.sqrt(2*n-1)/n*(2*x/L-1)*self.f[n-1]-(n-1)/n/np.sqrt(2*n-3)*self.f[n-2]))
        self.m=len(self.f)
        if w: self.w=w
        else: # Calculate midpoint weights these telescope so that w.sum()=L
            self.w=np.zeros_like(x)
            n=len(x)
            self.w=np.zeros((n,),dtype=np.float64)
            self.w[slice(1,n-1)]  = (x[slice(2,n)]-x[slice(0,n-2)])/2
            self.w[0]   =     (x[0]+x[1])/2
            self.w[n-1] = L - (x[n-1]+x[n-2])/2

    # Function to coefficients
    def f2c(self,f):
        c=np.zeros((self.m,),dtype=np.float64)
        for n in range(self.m): c[n]=(f*self.f[n]*self.w).sum()
        return(c)

    def c2f(self,c):
        f=np.zeros_like(self.f[0])
        for n in range(self.m): f+=c[n]*self.f[n]
        return(f)

    
