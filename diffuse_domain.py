import numpy as np
import math

pi=math.pi

def F(x,w): # x>0 is interior of region
    F = 1/(1+np.exp(-2*w*x)) #(1+tanh(w*x))/2;
    return(F)

def iF(x,w): # 1-F(x,w)
    F = 1/(1+np.exp(2*w*x)) #(1-tanh(w*x))/2;
    return(F)

def f(x,w):  # diff(F(x),x)
    f = w/(2*np.cosh(w*x)**2)
    return(f)

def lf(x,w):  # diff(log(F(x)),x)
    f = 2*w/(1+np.exp(2*w*x))
    return(f)

#These are indicator functions to use at a boundary
#The derivative is odd if n is odd and even if n is even
def F1(x,w,e=5): # Odd boundary indicator function
    F=np.where(x<=0,np.exp(-e),np.exp(-e*(2*(2*np.arctan(np.exp(-w*x))-np.exp(-w*x)))/(pi-2)))
    return(F)

def lf1(x,w,e=5): # Derivative of boundary indicator function
    f = np.where( x<=0, 0,4*e*w/(pi-2)*np.sinh(w*x)/(1+np.exp(2*w*x)))
    return(f)

def F2(x,w,e=5):
    F=np.where(x<=0,np.exp(-e),np.exp(-e*(np.exp(-2*w*x)+4*w*x-np.log(1+np.exp(4*w*x)))/(1-np.log(2))))
    return(F)
               
def lf2(x,w,e=5): # Derivative of boundary indicator function
    f = np.where( x<=0, 0,8*e*w/(1-np.log(2))*np.sinh(w*x)**2/(1+np.exp(4*w*x)))
    return(f)


