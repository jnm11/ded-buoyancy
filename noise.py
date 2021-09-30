#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test evolution of noise field

"""
Test evolution of noise field in 3d

Usage: 
   noise.py [options] NAME

Options:
  -L, --Length=<>       : Box length                      [default:  2]
  -W, --Width=<>        : Box width                       [default:  4]
  -H, --Height=<>       : Box height                      [default:  8]
  -N, --Nx=<>           : x grid points                   [default: 16]
      --Ny=<>           : y grid points                   [default: 32]
      --Nz=<>           : z grid points                   [default: 64]
      --Tz=<>           : type of z grid                  [default: Fourier]
      --dt=<>           : type of z grid                  [default: 0.01]
   -n, --noise=<>       : amplitude    for noise forcing  [default: 0.1]
      --noiseL=<>       : length scale for noise forcing  [default: 0.0]
      --noiseT=<>       : time scale for noise forcing    [default: 1.0]
      --noiseR=<>       : noise recording cadence         [default: 0.01]
      --AA=<>           : anti-aliasing default is 3/2    [default: 1.5]
"""

import time
import numpy as np
from docopt import docopt
from mpi4py import MPI
import os
import h5dict
from dedalus import public as de
from dedalus.core   import operators
from pathlib import Path

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

import logging
logger=logging.getLogger(__name__)

comm      = MPI.COMM_WORLD
mpirank   = comm.Get_rank()
mpisize   = comm.Get_size()

name =  np.str(args['NAME'])
dt        =  np.float(args['--dt'])
Nx        =    np.int(args['--Nx'])
Ny        =    np.int(args['--Ny'])
Nz        =    np.int(args['--Nz'])
AA        =  np.float(args['--AA'])
H         =  np.float(args['--Height'])
W         =  np.float(args['--Width'])
L         =  np.float(args['--Length'])
noise     =  np.float(args['--noise'])
noiseT    =  np.float(args['--noiseT'])
noiseL    =  np.float(args['--noiseL'])
noiseR    =  np.float(args['--noiseR'])
Tz        =    np.str(args['--Tz'])

sType='noise'
bdir = Path(os.environ['DEDALUS_DATA']) / sType
ddir = bdir / name
if mpirank==0:
    if not ddir.is_dir():
        ddir.mkdir()
comm.Barrier()
 
dx=L/Nx;
dy=W/Ny;
dz=H/Nz;
dV=dx*dy*dz

param={ 'Nx':Nx,'Ny':Ny,'Nz':Nz,
        'L':L,'W':W,'H':H,
        'dx':dx,'dy':dy,'dz':dz,
        'dt':dt,
        'AA':AA,
        'sType':sType,
        'name':name,
        'noise':noise,
        'noiseL':noiseL,
        'noiseR':noiseR,
        'noiseT':noiseT
}
if mpirank==0:
    fnp=ddir / 'param.h5'
    logger.info("Writing parameters to {}".format(fnp))
    h5dict.save_dict_to_hdf5(param,fnp)


rand = np.random.RandomState()#(seed=42)

# return a random field scaled by sqrt of time step in grid space
# the mean vaue doesn't matter since we are treating as  stream function and we filter out high frequencies
# see internal_heat.py

#nfx = domain.new_field(name='nfx')
#nfx.set_scales(3/2)
#nfx['c']
#fL=2

#kx = domain.bases[0].wavenumbers[:,np.newaxis,np.newaxis]/np.float(L/fL)
#ky = domain.bases[1].wavenumbers[np.newaxis,:,np.newaxis]/np.float(W/fL)
#kz = domain.bases[2].wavenumbers[np.newaxis,np.newaxis,:]/np.float(H/fL)
#kf = np.exp(-(kx**2+ky**2+kz**2))
#kgshape = domain.dist.coeff_layout.global_shape(scales=AA)
#kslices = domain.dist.coeff_layout.slices(scales=AA)
#def noise_zmean1(deltaT):
#    #    print(dir(domain.dist))
#    nfx['c'] = (rand.standard_normal(kgshape)[kslices]+1J*rand.standard_normal(kgshape)[kslices])*kf/np.sqrt(max(1e-12,deltaT))
#    #    nfx['g']
#    #    nfx.set_scales(3/2, keep_data=True)
#    x=nfx['g']
#    #    nfx.set_scales(3/2)
#    #    nfx['c']
#    return(x)


logger.info("Running noise.py: Noise evolution")
logger.info("MPI rank: {:d}, size: {:d}".format(mpirank,mpisize))
logger.info("name:            {}".format(name))
logger.info("Stype:           {}".format(sType))
logger.info("Tz:              {}".format(Tz))
logger.info("dt:              {}".format(dt))
logger.info("dx:    {:7.5f},    dy: {:7.5f},    dz: {:7.5f}".format(dx, dy,  dz))
logger.info("noise: {:7.5f} noiseL: {:7.5f} noiseT: {:7.5f} noiseR: {:7.5f}".format(noise,noiseL,noiseT,noiseR))

x_basis=de.Fourier(  'x', Nx, interval=(0, L), dealias=AA)
y_basis=de.Fourier(  'y', Ny, interval=(0, W), dealias=AA)
if   Tz=="Cheb":    z_basis=de.Chebyshev('z', Nz, interval=(0, H), dealias=AA)
elif Tz=="Fourier": z_basis=de.Fourier(  'z', Nz, interval=(0, H), dealias=AA)
else:   logger.error("Unknown z grid {}".format(Tz))

domain=de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)

# Gaussian noise scaled by 1/sqrt(dt*dV)*noise
# Adjust the scaling so that it sholdn't depend to much on the grid but it will depend on noiseL and noiseT
def anoise(sc):
    gshape = domain.dist.grid_layout.global_shape(scales=sc)
    slices = domain.dist.grid_layout.slices(scales=sc)
    x = rand.standard_normal(gshape)[slices]
    return(x)

def gnoise(deltaT):
    x = anoise(AA)*(noise*noiseL*np.sqrt(noiseL*noiseT/max(1e-12,deltaT*dV)))
    return(x)

zmn = operators.GeneralFunction(domain,'g',gnoise,args=[])
zmn.original_args = [dt]

if Tz=="Cheb":
    problem = de.IVP(domain, variables=['noisef','noisefz'], time='t')
    problem.substitutions['LP(X,Y)'] = '(dx(dx(X)) + dy(dy(X)) + dz(Y))'
    problem.substitutions['CD(X,Y)'] = '( u*dx(X)   + v*dy(X)  +  w*Y)'
    problem.add_equation("noisefz - dz(noisef) = 0")
    problem.add_bc("right(noisefz)  = 0")
    problem.add_bc("left(noisefz)   = 0")
elif Tz=="Fourier":
    problem = de.IVP(domain, variables=['noisef'], time='t')
    problem.substitutions['noisefz'] = '0'
    problem.substitutions['LP(X,Y)'] = '(dx(dx(X)) + dy(dy(X)) + dz(dz(X)))'
    problem.substitutions['CD(X,Y)'] = '( u*dx(X)  +  v*dy(X)  +  w*dz(X))'

problem.parameters['fnoise']  = zmn
problem.parameters['L']       = L
problem.parameters['W']       = W
problem.parameters['H']       = H
problem.parameters['noiseT']  = noiseT
problem.parameters['noiseL']  = noiseL


problem.add_equation("noiseT*dt(noisef) + noisef - LP(noisef,noisefz)*noiseL**2 = fnoise")
 
# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Integration parameters
solver.stop_sim_time  = 10
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

nf=solver.state['noisef']
nf['g']=noise*anoise(1)

analysis_tasks = []
s = solver.evaluator.add_file_handler(ddir / 'noise', sim_dt=noiseR, parallel=True,max_writes=np.inf,mode='overwrite')
s.add_task("noisef", scales=1, name='noisef')
analysis_tasks.append(s)

sim_time_start=solver.sim_time
t1= np.float(time.time())
while solver.ok: 
    zmn.args = [dt]
    solver.step(dt)
    t2= np.float(time.time())
    if t2>t1+30:
        logger.info('t: {:7.3f}'.format(solver.sim_time))
        t1=t2
        
