#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Find fourier smoothing numbers

Usage: 
   ded_find_smoothing.py [options] NX NZ

Options:
  -L, --Length=<>      : Length [default: 1]
  -H, --Height=<>      : Height [default: 1]
"""

#ded_find_smoothing.py 256 64
import numpy as np
from docopt import docopt
import jutil as ju
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools  import post
from polytrope.tools.checkpointing import Checkpoint

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

Nx =  np.int(args['NX'])
Nz =  np.int(args['NZ'])
H  =  np.float(args['--Height'])
L  =  np.float(args['--Length'])


dx=L/Nx;
dz=L/Nz;

Sx=100/dx
Sz=100/dz

x_basis=de.Fourier(  'x', Nx, interval=(0, L), dealias=3/2)
z_basis=de.Chebyshev('z', Nz, interval=(0, H), dealias=3/2)
domx=de.Domain([x_basis], grid_dtype=np.float64)
domz=de.Domain([x_basis], grid_dtype=np.float64)

# Get pointers to the grid coordinates
x = domx.grid(0)
z = domz.grid(0)

fx = domx.new_field(name='fx')
fz = domz.new_field(name='fz')


x1=1*L/2
x2=3*L/4
z1=1*H/2
z2=3*H/4


gx=ju.indicatorp(x,x1,x2,Sx,L)
gz=ju.indicatorp(z,z1,z2,Sz,H)

fx['g']  = gx   # This region force uniform velocity profile
fz['g']  = gz   # This region force uniform velocity profile

hx=fx['g']
hz=fz['g']

ex=max(abs(hx-gx))
ez=max(abs(hz-gz))

mingx=-min(gx)
maxgx= max(gx)-1
minfx=-min(hx)
maxfx= max(hz)-1
print("Sx={:8.4f}".format(Sx))
print("gx=[{:11.2f},{:11.2f}]".format(mingx,maxgx))
print("fx=[{:11.2f},{:11.2f}]".format(minfx,maxfx))
print("ex={:11.2f}, ez={:11.2f}]".format(ex,ez))



