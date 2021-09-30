#2d incompressible gravity current simualtion using streamfunction vorticity formulation
#u =  f_y       x velocity
#v = -f_x       y velocity
#w = -f_xx-f_yy vorticity

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

Re=4000.
Sc=4000.
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
problem = de.IVP(domain, variables=['b','f','w','u','by','wy'], time='t')
problem.meta['b','by']['x']['parity'] = 1
problem.meta['u','w','wy','f']['x']['parity'] = -1

Re=2000.
Sc=2000.

problem.parameters['nu']    = 1/Re
problem.parameters['kappa'] = 1/Sc
problem.parameters['Lx']    = Lx
problem.parameters['Ly']    = Ly

print("Runing gc5.py: stream-function vorticity formulation with odd symmetry")
print("Lx={}, Ly={}, Nx={}, Ny={}, Re={}, Sc={}".format(Lx,Ly,Nx,Ny,Re,Sc))

# Section 2D Navier-Stokes Vorticity Stream function formulation in vorticity-2d.mw
# Boundary conditions for this formulation are not straight forward
# Normally somthing complicated for the vorticity w must be specified but this is not necessary with this scheme
problem.add_equation("dt(b) - kappa*(dx(dx(b)) + dy(by))         = - u*dx(b) + dx(f)*by")
problem.add_equation("dt(w) -    nu*(dx(dx(w)) + dy(wy)) + dx(b) = - u*dx(w) + dx(f)*wy")
problem.add_equation("w  + dy(u) + dx(dx(f)) = 0")
problem.add_equation("u  - dy(f) = 0")
#problem.add_equation("v  + dx(f) = 0")
problem.substitutions['v'] = "-dx(f)"
problem.add_equation("wy - dy(w) = 0")
problem.add_equation("by - dy(b) = 0")

problem.add_bc("left(by)   = 0")
problem.add_bc("right(by)  = 0")
problem.add_bc("left(f)    = 0") 
problem.add_bc("right(f)   = 0") 
problem.add_bc("left(u)    = 0")
problem.add_bc("right(u)   = 0")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Initial conditions

x = domain.grid(0)
y = domain.grid(1)

b=solver.state['b']
u=solver.state['u']
#v=solver.state['v']
w=solver.state['w']
f=solver.state['f']
by=solver.state['by']
wy=solver.state['wy']

b['g'] = sp.erf((x-Lx/2)*Nx/Lx/2)
b.differentiate('y', out=by)

# Integration parameters
solver.stop_sim_time = 12.
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

# Read the output directory from the environment variable DEDALUS_DATA
data_dir=os.environ['DEDALUS_DATA']+'/gc'
# Analysis -- output files
# scale increases the resolution of the output datafiles
print("Writing data to directory '{}'".format(data_dir))
snap = solver.evaluator.add_file_handler(data_dir, sim_dt=0.1)#, max_writes=10)
snap.add_task("b" , scales=2, name='b')
snap.add_task("f",  scales=2, name='f')


# CFL
CFL = flow_tools.CFL(solver, initial_dt=1e-2, cadence=5, safety=0.8,max_change=1.5, min_change=0.5, max_dt=0.1)
CFL.add_velocities(('u','v')) # maybe calculate with w*sin(y/Ly) since in y direction the grid-spacing is non uniform

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("sqrt(u*u + v*v)*Ly/nu/2", name='Re')

# Main loop
end_init_time = time.time()
logger.info('Initialization time: %f' %(end_init_time-start_init_time))
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



