#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""beta plane barotropic vorticity model.

Solve the barotropic vorticity equation in two dimensions

    D/Dt[ω] = 0                                                             (1)

where ω = ζ + f is absolute vorticity.  ζ is local vorticity ∇ × u and
f is global rotation.

Assuming an incompressible two-dimensional flow u = (u, v),
the streamfunction ψ = ∇ × (ψ êz) can be used to give (u,v)

    u = -∂/∂y[ψ]         v = ∂/∂x[ψ]                                        (2)

and therefore local vorticity is given by the Poisson equation

    ζ = ∆ψ                                                                  (3)

Since ∂/∂t[f] = 0 equation (1) can be written in terms of the local vorticity

        D/Dt[ζ] + u·∇f = 0
    =>  D/Dt[ζ] = -vβ                                                       (4)

using the beta-plane approximation f = f0 + βy.  This can be written entirely
in terms of the streamfunction and this is the form that will be solved
numerically.

    D/Dt[∆ψ] = -β ∂/∂x[ψ]                                                   (5)

"""
import logging

import numpy as np
import matplotlib.pyplot as plt

from dedalus import public as de
from dedalus.extras import flow_tools

root = logging.root
for h in root.handlers:
    h.setLevel("INFO")

logger = logging.getLogger(__name__)

N = 256
Lx, Ly = (1., 1.)
nx, ny = (N, N)
beta = 8.0
U = 0.0

# setup the domain
x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', ny, interval=(0, Ly), dealias=3/2)
domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)

problem = de.IVP(domain, variables=['psi'])


# solve the problem from the equations
# ζ = Δψ
# ∂/∂t[∆ψ] + β ∂/∂x[ψ] = -J(ζ, ψ)

# Everytime you ask for one of the expression on the left, you will get the expression on the right.
problem.substitutions['zeta'] = "  d(psi,x=2) + d(psi,y=2) "
problem.substitutions['u']    = " -dy(psi) "
problem.substitutions['v']    = "  dx(psi) "

# This pattern matches for the 'thing' arguements. They don't have to be called 'thing'.
problem.substitutions['L(thing_1)']         = "  d(thing_1,x=2) + d(thing_1,y=2) "
problem.substitutions['J(thing_1,thing_2)'] = "  dx(thing_1)*dy(thing_2) - dy(thing_1)*dx(thing_2) "

# You can combine things if you want
problem.substitutions['HD(thing_1)']         = "  -D*L(L(thing_1)) "

problem.parameters['beta'] = beta
problem.parameters['U']    = U
problem.parameters['D']   = 0.01 # hyperdiffusion coefficient

problem.add_equation("dt(zeta) + beta*v  = J(psi,zeta) ", condition="(nx != 0) or  (ny != 0)")
problem.add_equation("psi = 0",                                         condition="(nx == 0) and (ny == 0)")


solver = problem.build_solver(de.timesteppers.CNAB2)
solver.stop_sim_time  = np.inf
solver.stop_wall_time = np.inf
solver.stop_iteration = 1000

# vorticity & velocity are no longer states of the system. They are true diagnostic variables.
# But you still might want to set initial condisitons based on vorticity (for example).
# To do this you'll have to solve for the streamfunction.

# This will solve for an inital psi, given a vorticity field.
init = de.LBVP(domain, variables=['init_psi'])

gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=42)
noise = rand.standard_normal(gshape)[slices]

init_vorticity = domain.new_field()
init_vorticity.set_scales(1)

x,y = domain.grids(scales=1)

init_vorticity['g'] =  (0.5)*noise +  3*np.exp( - 80* ( (x-Lx/2)**2 + (y-Ly/4)**2 ) )

init.parameters['init_vorticity'] = init_vorticity

init.add_equation(" d(init_psi,x=2) + d(init_psi,y=2) = init_vorticity ", condition="(nx != 0) or  (ny != 0)")
init.add_equation(" init_psi = 0",                                        condition="(nx == 0) and (ny == 0)")

init_solver = init.build_solver()
init_solver.solve()

psi = solver.state['psi']
psi['g'] = init_solver.state['init_psi']['g']

# Now you are ready to go.
# Anytime you ask for zeta, u, or v they will be non-zero because psy is non-zero.

dt = 1e-3 #Lx/nx
# You can set parameters to limit the size of the timestep.
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=2,
                     max_change=1.5, min_change=0.5, max_dt=10*dt)
CFL.add_velocities(('u','v'))

# I don't really know what this is showing. But it makes somethign that looks abotu right.
#plt.ion()
#fig, axis = plt.subplots(figsize=(10,5))
#p = axis.imshow(init_vorticity['g'].T, cmap=plt.cm.YlGnBu)
#plt.pause(1)

logger.info('Starting loop')
while solver.ok:
    # dt = cfl.compute_dt()   # this is returning inf after the first timestep
    # print(dt)
    solver.step(dt)
#    if solver.iteration % 1000 == 0:
        # This won't work any more.
        
            # Update plot of scalar field
            #p.set_data(zeta['g'].T)
            #p.set_clim(np.min(zeta['g']), np.max(zeta['g']))
        
        # There are several ways to see the output as you go.
        # I recommend creating a file handler and saving the output.
        # If you are worried about speed,
        # then plotting in real time is not going to get you to your goals.
        
#logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
#plt.pause(0.001)

# Print statistics
#logger.info('Iterations: %i' %solver.iteration)
logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))

