"""
1D KS equation

This script should be ran serially (because it is 1D), and creates a space-time
plot of the computed solution.

"""

import numpy as np
import matplotlib.pyplot as plt

from dedalus import public as de
from dedalus.extras.plot_tools import quad_mesh, pad_limits

import logging
logger = logging.getLogger(__name__)


# Bases and domain
x_basis = de.Fourier('x', 4*2048, interval=(-100, 100), dealias=3/2)
domain = de.Domain([x_basis], np.float64)

# Problem
problem = de.IVP(domain, variables=['u'])
problem.parameters['a'] = 1
problem.parameters['b'] = 1
problem.add_equation("dt(u) + a*d(u,x=2) + b*d(u,x=4) = -u*dx(u)")


# Build solver
solver = problem.build_solver(de.timesteppers.SBDF2)
solver.stop_wall_time = np.inf
solver.stop_iteration = 200000
solver.stop_sim_time  = np.inf

# Random initial conditions

x = domain.grid(0)
u = solver.state['u']

u['g'] = np.random.normal(0,10,len(x))

#Fouier filter high wavenumbers
u.set_scales(1/4, keep_data=True)
u['c']
u['g']
u.set_scales(1, keep_data=True)


# Store data for final plot

u.set_scales(1, keep_data=True)

plot_freq = 10
nt        = 1 + solver.stop_iteration//plot_freq
u_plot    = np.zeros((nt,len(u['g'])))
t_plot    = np.zeros(nt)
u_plot[0] = np.copy(u['g'])

# Main loop
dt = 2e-3
while solver.ok:
    solver.step(dt)
    
    if solver.iteration % plot_freq == 0:
        
        plot_index = solver.iteration//plot_freq
        
        u.set_scales(1, keep_data=True)
        u_plot[plot_index] = np.copy(u['g'])
        t_plot[plot_index] = solver.sim_time
    
    if solver.iteration % 100 == 0:
        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))

# Create space-time plot

xmesh, ymesh = quad_mesh(x=x, y=t_plot)
plt.figure()
plt.pcolormesh(xmesh, ymesh, u_plot, cmap='RdBu_r')
plt.clim((-4,4))
plt.axis(pad_limits(xmesh, ymesh))
plt.colorbar()
plt.xlabel('space (x)')
plt.ylabel('time (t)')
plt.title('$\partial_{t} u \ + \ \partial_{x}^{2} u \ + \ \partial_{x}^{4} u \ = \ - u \partial_{x} u$')
plt.savefig('space_time_diagram.png')

