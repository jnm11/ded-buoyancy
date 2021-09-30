# coding: utf-8

import math as m
import numpy as np
from mpi4py import MPI
import time
import glob, os

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.extras import plot_tools
from dedalus.tools import post

import logging
root = logging.root
for h in root.handlers:
    h.setLevel("INFO") 
logger = logging.getLogger(__name__)

from numpy.polynomial import chebyshev as cb

def rescale(x,a,b,c,d):
    return c + (d - c)*(x - a)/(b - a)

def cheb_mode(x, m, a, b):
    x_scaled = rescale(x,a,b,-1,1)
    return cb.chebval(x_scaled, np.append(np.zeros(m),1))

def cheb_vals(x, M, interval=(-1,1)):
    a, b = interval
    return np.array([cheb_mode(x,m,a,b) for m in range(M)])

def sin_mode(x, m, a, b):
    x_scaled = rescale(x,a,b,0,2*np.pi)
    return np.sin(m*x_scaled)

def cos_mode(x, m, a, b):
    x_scaled = rescale(x,a,b,0,2*np.pi)
    return np.cos(m*x_scaled)

def fourier_mode(x, m, a, b, c):
    if m == 0: return c.real*np.ones_like(x)
    else: return 2*(c.real*cos_mode(x,m,a,b) - c.imag*sin_mode(x,m,a,b))

def fourier_cos_vals(x, M, interval=(0,2*np.pi)):
    a, b = interval
    return np.array([fourier_mode(x,m,a,b,1) for m in range(M)]).T

def fourier_sin_vals(x, M, interval=(0,2*np.pi)):
    a, b = interval
    return np.array([fourier_mode(x,m,a,b,1j) for m in range(M)]).T

def trig_cos_vals(x, M, interval=(0,2*np.pi)):
    a, b = interval
    return np.array([cos_mode(x,m/2,a,b) for m in range(M)])

def trig_sin_vals(x, M, interval=(0,2*np.pi)):
    a, b = interval
    return np.array([sin_mode(x,m/2,a,b) for m in range(M)])

def interpolate_1D(u, x, comm=None, basis_type=('Chebyshev')):
    # Get bases and shapes
    domain = u.domain
    xbasis,  = domain.bases
    gcshape = xbasis.coeff_size
    
    # Build correct interpolation functions
    if basis_type == 'Chebyshev':
        xfunc = lambda z : cheb_vals(x, gcshape, interval=xbasis.interval)
    elif basis_type == 'SinCos':
        u_parity = u.meta['z']['parity']
        if u_parity == 1:   xfunc = lambda x : trig_cos_vals(x, gcshape, interval=xbasis.interval)
        elif u_parity ==-1: xfunc = lambda x : trig_sin_vals(x, gcshape, interval=xbasis.interval)
    
    # Get grid values of the modes
    xs = xfunc(x.flatten())
    
    # Perform the transform
    A = u['c'].copy()
    G = np.dot(A, xs)
    return G

def interpolate_2D(u, x, z, comm=None, basis_types=('Fourier','Chebyshev')):
    # Get mpi communications 
    domain = u.domain
    comm = domain.dist.comm
    rank, size = comm.rank, comm.size
    
    # Get bases and shapes
    xbasis,zbasis = domain.bases
    lcshape = domain.dist.coeff_layout.local_shape(1)
    gcshape = domain.dist.coeff_layout.global_shape(1)
    
    # Build global coefficients in local processor with MPI
    if size > 1:
        # prepare to send local coefficients to rank 0
        sendbuf = u['c'].copy()
        recvbuf = None
        if rank == 0: recvbuf = np.empty([size, *lcshape], dtype=np.complex128)
        comm.Gather(sendbuf, recvbuf, root=0)
        
        # send global coefficients to every processor from rank 0
        if rank == 0: gcoeffs = np.reshape(recvbuf,gcshape)
        else: gcoeffs = np.empty(gcshape, dtype=np.complex128)
        comm.Bcast(gcoeffs, root=0)    
    else: gcoeffs = u['c']    
    
    # Build correct z basis interpolation functions
    if basis_types[1] == 'SinCos':
        u_parity = u.meta['z']['parity']
        if u_parity == 1:   zfunc = lambda z : trig_cos_vals(z, gcshape[1], interval=zbasis.interval)
        elif u_parity ==-1: zfunc = lambda z : trig_sin_vals(z, gcshape[1], interval=zbasis.interval)
    elif basis_types[1] == 'Chebyshev':
        zfunc = lambda z : cheb_vals(z, gcshape[1], interval=zbasis.interval)
    
    # Get the grid values of all the modes (split for Fourier)
    xsc = fourier_cos_vals(x.flatten(), gcshape[0], interval=xbasis.interval)
    xss = fourier_sin_vals(x.flatten(), gcshape[0], interval=xbasis.interval)
    zs = zfunc(z.flatten())

    # Perform the transform
    B, C = gcoeffs.real, gcoeffs.imag
    F = np.dot(xsc,B) + np.dot(xss,C)
    G = np.dot(F,zs)
    return G

def togrid(A, xsc, xss, zs):
    B, C = A.real, A.imag
    F = np.dot(xsc,B) + np.dot(xss,C)
    return np.dot(F,zs)

def interperic(u,x,z):
    xsc = xvalsc(x)
    xss = xvalss(x)
    zs = zvals(z)
    A = u['c'].copy()
    return togrid(A, xsc, xss, zs)

# # Testing communication of numpy arrays
# comm = MPI.COMM_WORLD
# rank, size = comm.rank, comm.size

# # Create bases and domain
# resx,resz = 1024,1024
# a, b = -1, 4
# xbasis = de.Fourier('x', resx, interval=(a,b), dealias=1)
# zbasis = de.SinCos('z', resz, interval=(a,b), dealias=1)
# domain = de.Domain([xbasis,zbasis], grid_dtype=np.float64)

# x,z = domain.grids(scales=domain.dealias)
# xx, zz = np.meshgrid(x,z,indexing='ij')
# kx, kz = xbasis.wavenumbers, np.arange(zbasis.coeff_size)
# kxx,kzz = np.meshgrid(kx,kz,indexing='ij')

# xg,zg = xbasis.grid(), zbasis.grid()
# xxg,zzg = np.meshgrid(xg,zg,indexing='ij')
# gslice = domain.dist.grid_layout.slices(1)

# u = domain.new_field(scales=domain.dealias)
# u.meta['z']['parity'] = -1
# nd = 4
# δ = 1
# u['g'] = np.exp(-(zz)**2/δ**2)*((xx-a)**nd)*(b-xx)**nd

# x1, z1 = np.linspace(a,b,10000), np.linspace(a,b,10000)
# xx1, zz1 = np.meshgrid(x1,z1,indexing='ij')

# uinterp = interpolate(u, x, z, basis_types=('Fourier','SinCos'))

# print(np.sum((u['g'] - uinterp))/(resx*resz))
















