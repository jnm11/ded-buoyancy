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
Lx=16. # Domain length
Ly=1.  # Domain height


#Re=1000. give a Reynolds number based on the front velocity and current height of around 420 

Re=2000.
Sc=2000.
if Re<=2000:
    Nx = 512 # x fourier modes
else:
    if Re<=5000:
        Nx = 1024 # x fourier modes

Ny=round(2*Ly/Lx*Nx)

# Create bases and domain
start_init_time=time.time()
x_basis=de.SinCos(   'x', Nx, interval=(0, Lx), dealias=3/2)
y_basis=de.Chebyshev('y', Ny, interval=(0, Ly), dealias=3/2)
domain=de.Domain([x_basis, y_basis], grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','b','u','v','by','uy','vy'], time='t')
problem.meta['p','b','v','by','vy']['x']['parity'] = 1
problem.meta['u','uy']['x']['parity'] = -1

problem.parameters['nu']    = 1/Re
problem.parameters['kappa'] = 1/Sc
problem.parameters['Lx']    = Lx
problem.parameters['Ly']    = Ly

print("Runing ded_gc_primtive.py: Boussinesq lock-exchange gravity current in primitive formulation with odd symmetry")
print("Lx={}, Ly={}, Nx={}, Ny={}, Re={}, Sc={}".format(Lx,Ly,Nx,Ny,Re,Sc))

problem.add_equation("dx(u) + vy = 0")
problem.add_equation("dt(b) - kappa*(dx(dx(b)) + dy(by))             = - u*dx(b) - v*by")
problem.add_equation("dt(u) -    nu*(dx(dx(u)) + dy(uy)) + dx(p)     = - u*dx(u) - v*uy")
problem.add_equation("dt(v) -    nu*(dx(dx(v)) + dy(vy)) + dy(p) + b = - u*dx(v) - v*vy")
problem.add_equation("by - dy(b) = 0")
problem.add_equation("uy - dy(u) = 0")
problem.add_equation("vy - dy(v) = 0")

problem.add_bc("left(by)   = 0")
problem.add_bc("right(by)  = 0")
problem.add_bc("left(uy)   = 0")
problem.add_bc("right(uy)  = 0")
problem.add_bc("left(v)    = 0", condition="(nx != 0)")
problem.add_bc("right(v)   = 0")
problem.add_bc("integ_y(p) = 0", condition="(nx == 0)")

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
p=solver.state['p']
by=solver.state['by']
uy=solver.state['uy']
vy=solver.state['vy']


# Set the bouyancy field to lock exchange with a slightly smoothed edge to prevent ringing
b['g'] = sp.erf((x-Lx/2)*Nx/Lx)
b.differentiate('y', out=by)

# Integration parameters
solver.stop_sim_time = 13.05
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf


# Read the output directory from the environment variable DEDALUS_DATA
data_dir=os.environ['DEDALUS_DATA']+'/gc'
print("Writing data to directory '{}'".format(data_dir))

# Set the output data to writing buoyancy every 0.1 time units
# scale increases the resolution of the output datafiles
snap = solver.evaluator.add_file_handler(data_dir, sim_dt=0.1)#, max_writes=10)
snap.add_task("b" , scales=2, name='b')

# CFL
CFL = flow_tools.CFL(solver, initial_dt=1e-2, cadence=5, safety=0.8, max_change=1.5, min_change=0.5, max_dt=0.1)
CFL.add_velocities(('u','v')) # maybe calculate with w*sin(y/Ly) since in y direction the grid-spacing is non uniform

# Flow properties Print som information about the simulation every cadence iterations
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("sqrt(u*u + v*v)*Ly/nu/2", name='Re')

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



