import numpy as np
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger=logging.getLogger(__name__)

#import matplotlib.pyplot as plt

Nz = 64
Nx = 512

# Parameters
Lx, Lz = (16.,1.)

# Create bases and domain
start_init_time=time.time()
x_basis=de.SinCos(   'x', Nx, interval=(0, Lx), dealias=3/2)
z_basis=de.Chebyshev('z', Nz, interval=(0, Lz), dealias=3/2)
domain=de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','b','u','w','bz','uz','wz'], time='t')
problem.meta['p','b','w','bz','wz']['x']['parity'] = 1
problem.meta['u','uz']['x']['parity'] = -1
problem.parameters['nu']    = 1/2000.
problem.parameters['kappa'] =x 1/2000.
problem.parameters['gx']    = 0
problem.parameters['gz']    = -1
problem.parameters['Lx']    = Lx
problem.parameters['Lz']    = Lz


problem.add_equation("dx(u) + wz = 0")
problem.add_equation("dt(b) - kappa*(dx(dx(b)) + dz(bz))                 +  pwz = - dx(pu)  ")
problem.add_equation("dt(pu) -    nu*(dx(dx(u)) + dz(uz)) + dx(p)        + puwz = - dx(u*pu)")
problem.add_equation("dt(pw) -    nu*(dx(dx(w)) + dz(wz)) + dz(p) - gz*b + pwwz = - dx(w*pw)")
problem.add_equation("pu = u*b")
problem.add_equation("pw = w*b")
problem.add_equation("pwz - dz(pw) = 0")
problem.add_equation("puwz = dz(pu*w+pw*u)/2")
problem.add_equation("pwwz = dz(pw*w)")

problem.add_bc("left(bz)   = 0")
problem.add_bc("right(bz)  = 0")
#problem.add_bc("left(dz(u)) = 0")
#problem.add_bc("right(u)= 0")
problem.add_bc("left(u) = 0")
problem.add_bc("right(dz(u))= 0")
problem.add_bc("left(w)    = 0")
problem.add_bc("right(w)   = 0", condition="(nx != 0)")
problem.add_bc("integ_z(p) = 0", condition="(nx == 0)")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Initial conditions

x = domain.grid(0)
z = domain.grid(1)

b=solver.state['b']
u=solver.state['u']
w=solver.state['w']
p=solver.state['p']
bz=solver.state['bz']
uz=solver.state['uz']
wz=solver.state['wz']

b['g'] = -np.tanh((x-Lx/2)*Nx/Lx)
b.differentiate('z', out=bz)
                                
# Integration parameters
solver.stop_sim_time = 20.
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

# Analysis -- output files
# scale increases the resolution of the output datafiles
z = domain.grid(1)
snap = solver.evaluator.add_file_handler('/Users/jnm/data/Dedalus/gc', sim_dt=0.1)#, max_writes=10)
snap.add_task("b" , scales=2, name='b')
#snap.add_task("u", scales=2, name='u')
#snap.add_task("w", scales=2, name='w')

# CFL
CFL = flow_tools.CFL(solver, initial_dt=1e-2, cadence=5, safety=0.8,max_change=1.5, min_change=0.5, max_dt=0.1)
CFL.add_velocities(('u','w')) # maybe calculate with w*sin(z/Lz) since in z direction the grid-spacing is non uniform

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u*u + w*w)*Lz/nu/2", name='Re')


#x = domain.grid(0)
#y = domain.grid(1)
#print("len(x)",len(x))
#print("len(y)",len(y))
#xm, ym = np.meshgrid(x,y)
#print("len(xm)",len(xm))
#print("len(ym)",len(ym))
#plt.ion()
#fig, axis = plt.subplots(figsize=(10,5))
#axis.pcolormesh(xm, ym, b['g'].T, cmap='RdBu_r');
#axis.set_xlim([-Lx/2,Lx/2])
#axis.set_ylim([0,Lz])
#plt.draw()
#plt.pause(0.05)
 
#x = domain.grid(0,scales=domain.dealias)
#y = domain.grid(1,scales=domain.dealias)
#xm, ym = np.meshgrid(x,y)

# Main loop
end_init_time = time.time()
logger.info('Initialization time: %f' %(end_init_time-start_init_time))
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    while solver.ok:
        dt = CFL.compute_dt()
        solver.step(dt)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max Re = %f' %flow.max('Re'))
  #          axis.pcolormesh(xm, ym, b['g'], cmap='RdBu_r');
  #          plt.title('time %f' % solver.sim_time)
  #          plt.draw()
  #          fig.canvas.draw()
  #          fig.canvas.flush_events()
        
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*domain.dist.comm_world.size))



