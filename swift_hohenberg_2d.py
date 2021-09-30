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
x_basis = de.Fourier('x', 256, interval=(-10, 10), dealias=2)
y_basis = de.Fourier('y', 256, interval=(-10, 10), dealias=2)
domain = de.Domain([x_basis,y_basis], np.float64)

# Problem
problem = de.IVP(domain, variables=['u'])
problem.parameters['a'] = 1
problem.substitutions['L(var)'] = "d(var,x=2) + d(var,y=2) + a*var"
problem.add_equation("dt(u) + L(L(u)) = -u**3")


# Build solver
solver = problem.build_solver(de.timesteppers.SBDF2)
solver.stop_wall_time = np.inf
solver.stop_iteration = 1000
solver.stop_sim_time  = np.inf

# Random initial conditions

x, y = domain.grid(0), domain.grid(1)
u    = solver.state['u']

amplitude = 20
u['g']    = np.random.normal(0,amplitude,(x.shape[0], y.shape[1]))

#Fouier filter high wavenumbers
u.set_scales(1/4, keep_data=True)
u['c']
u['g']
u.set_scales(1, keep_data=True)


# Main loop
dt = 2e-3
while solver.ok:
    solver.step(dt)
    
    if solver.iteration % 10 == 0:
        logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))

# Create 2D plot

u.set_scales(1, keep_data=True)
x.shape, y.shape = len(x), len(x)

xmesh, ymesh = quad_mesh(x=x, y=y)
plt.figure()
plt.pcolormesh(xmesh, ymesh, u['g'] , cmap='RdBu_r')
plt.clim((-0.5,0.5))
plt.axis(pad_limits(xmesh, ymesh))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('$\partial_{t} u \ + \ ( \\nabla^{2} + 1)^{2} u \ = \  -\,u^{3}$')
plt.savefig('2d_diagram.png')

