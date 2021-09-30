import numpy as np
import scipy.special as sp
from mpi4py import MPI
import time
import os
from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger=logging.getLogger(__name__)

# Parameters
Lx=8.      # Domain length
Ly=0.125   # Domain height
Lz=1.      # Domain height


#Re=1000. give a Reynolds number based on the front velocity and current height of around 420 

Re=500.
Sc=500.
Nx = 256 # x fourier modes
Ny=round(2*Ly/Lx*Nx)
Nz=round(  Lz/Lx*Nx)

# Create bases and domain
start_init_time=time.time()
x_basis=de.SinCos(   'x', Nx, interval=(0, Lx), dealias=3/2)
y_basis=de.Fourier(  'y', Ny, interval=(0, Ly), dealias=3/2)
z_basis=de.Chebyshev('z', Nz, interval=(0, Lz), dealias=3/2)
domain=de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)

# 3D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','b','u','v','w','bz','uz','vz','wz'], time='t')
problem.meta['p','b','bz','v','vz','w','wz']['x']['parity'] = 1
problem.meta['u','uz']['x']['parity'] = -1

problem.parameters['nu']    = 1/Re
problem.parameters['kappa'] = 1/Sc
problem.parameters['Lx']    = Lx
problem.parameters['Ly']    = Ly

print("Runing ded_gc_primtive.py: Boussinesq lock-exchange gravity current in primitive formulation with odd symmetry")
print("Lx={}, Ly={}, Lz={}, Nx={}, Ny={}, Nz={}, Re={}, Sc={}".format(Lx,Ly,Lz,Nx,Ny,Nz,Re,Sc))

problem.add_equation("dx(u) + dy(v) + wz = 0")
problem.add_equation("dt(b) - kappa*(dx(dx(b)) + dy(dy(b)) + dz(bz))             = - u*dx(b) - v*dy(b) - w*bz")
problem.add_equation("dt(u) -    nu*(dx(dx(u)) + dy(dy(u)) + dz(uz)) + dx(p)     = - u*dx(u) - v*dy(u) - w*uz")
problem.add_equation("dt(v) -    nu*(dx(dx(v)) + dy(dy(v)) + dz(vz)) + dy(p)     = - u*dx(v) - v*dy(v) - w*vz")
problem.add_equation("dt(w) -    nu*(dx(dx(w)) + dy(dy(w)) + dz(wz)) + dz(p) - b = - u*dx(w) - v*dy(w) - w*wz")
problem.add_equation("bz - dz(b) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("vz - dz(v) = 0")
problem.add_equation("wz - dz(w) = 0")


problem.add_bc("left(bz)   = 0")
problem.add_bc("right(bz)  = 0")
problem.add_bc("left(u)    = 0")
problem.add_bc("right(u)   = 0")
problem.add_bc("left(v)    = 0")
problem.add_bc("right(v)   = 0")
problem.add_bc("integ_z(p) = 0", condition="(nx == 0)")
problem.add_bc("left(w)    = 0", condition="(nx != 0)")
problem.add_bc("right(w)   = 0")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Get pointers to the grid coordinates
x = domain.grid(0)
y = domain.grid(1)

# Get pointers to the simulation variables and set them to zero
b=solver.state['b']
u=solver.state['u']
v=solver.state['v']
w=solver.state['w']
p=solver.state['p']
bz=solver.state['bz']
uz=solver.state['uz']
vz=solver.state['vz']
wz=solver.state['wz']


# Set the bouyancy field to lock exchange with a slightly smoothed edge to prevent ringing
b['g'] = sp.erf((x-Lx/2)*Nx/Lx)
b.differentiate('z', out=bz)

# Integration parameters
solver.stop_sim_time = 6.05
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf


# Read the output directory from the environment variable DEDALUS_DATA
data_dir=os.environ['DEDALUS_DATA']+'/gc'
print("Writing data to directory '{}'".format(data_dir))

# Set the output data to writing buoyancy every 0.1 time units
# scale increases the resolution of the output datafiles
snap = solver.evaluator.add_file_handler(data_dir, sim_dt=0.1)
snap.add_task("integ(b, 'y')", name='b', scales=2)
snap.add_task("integ(u, 'y')", name='u', scales=2)
snap.add_task("integ(w, 'y')", name='w', scales=2)

# CFL
CFL = flow_tools.CFL(solver, initial_dt=1e-2, cadence=5, safety=0.8, max_change=1.5, min_change=0.5, max_dt=0.1)
CFL.add_velocities(('u','v','w')) # maybe calculate with w*sin(y/Ly) since in y direction the grid-spacing is non uniform

# Flow properties Print som information about the simulation every cadence iterations
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("sqrt(u*u + v*v + w*w)*Ly/nu/2", name='Re')

end_init_time = time.time()
logger.info('Initialisation time: %f' %(end_init_time-start_init_time))

try:
    logger.info('Starting loop')
    start_run_time = time.time()
    while solver.ok:
        dt = CFL.compute_dt()
        solver.step(dt)
        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
        logger.info('Max Re = %f' %flow.max('Re'))
        
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*domain.dist.comm_world.size))



