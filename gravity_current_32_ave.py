import numpy as np
import os
import matplotlib.pyplot as plt
import time
import shelve
from dedalus2.public import *

# Set domain
nx = 2048
nz = 128
Lx = 30.
Lz = 1.
x_basis = Fourier(nx,   interval=[-Lx, Lx])
z_basis = Chebyshev(nz, interval=[-Lz, Lz])
domain = Domain([x_basis, z_basis],grid_dtype = np.float64)

Gr = 1.25e6
Sc = 1.

nu = np.sqrt(Lz**2/Gr)
D = nu/Sc

boussinesq = ParsedProblem(axis_names=['x', 'z'],
                                    field_names=['p','u','w','s','uz','sz'],
                                    param_names=['D','nu'])

boussinesq.add_equation("dt(u) + dx(p) - nu*dx(dx(u)) - nu*dz(uz) = -u*dx(u) - w*uz")
boussinesq.add_equation("dt(w) + dz(p) + s - nu*dx(dx(w)) + nu*dx(uz) = -u*dx(w) + w*dx(u)")
boussinesq.add_equation("dx(u) + dz(w) = 0")
boussinesq.add_equation("dt(s) - D*dx(dx(s)) - D*dz(sz) = -u*dx(s) - w*sz")
boussinesq.add_equation("uz - dz(u) = 0")
boussinesq.add_equation("sz - dz(s) = 0")

boussinesq.add_right_bc("w = 0")
boussinesq.add_right_bc("u = 0")
boussinesq.add_right_bc("sz = 0")
boussinesq.add_left_bc("sz = 0")
boussinesq.add_left_bc("w = 0", condition="dx != 0")
boussinesq.add_left_bc("p = 0", condition="dx == 0")
boussinesq.add_left_bc("u = 0")

boussinesq.parameters['nu'] = nu
boussinesq.parameters['D'] = D

boussinesq.expand(domain, order=1)

ts = timesteppers.MCNAB2

# Build solver
solver = solvers.IVP(boussinesq, domain, ts)

# initial conditions
x  = domain.grid(0)
z  = domain.grid(1)
u  = solver.state['u']
w  = solver.state['w']
s  = solver.state['s']
sz = solver.state['sz']
p  = solver.state['p']

L = 15
Lw = (Gr*Sc**2)**(-1/4.)

import scipy.special as sp

s['g']  = 0.5 * (np.sign( x) + 1)*(1.-sp.erf((x-L)/Lw))/2.
s['g'] += 0.5 * (np.sign(-x) + 1)*(1.+sp.erf((x+L)/Lw))/2.

# start with temperature perturbations only
#A0 = 1e-6
#T['g'] += A0 * np.sin(np.pi * z/Lz) * np.random.randn(*T['g'].shape)

def grid_spacing(grid):
    diff = np.diff(grid)
    dg = np.empty_like(grid)
    dg[0] = diff[0]
    dg[-1] = diff[-1]
    for i in range(1, grid.size-1):
        dg[i] = min(diff[i], diff[i-1])

    return dg

dx = grid_spacing(x[:,0]).reshape((x.size, 1))
dz = grid_spacing(z[0,:]).reshape((1, z.size))

def cfl_dt(safety=1.):
    minut = np.min(np.abs(dx / solver.state['u']['g'].real))
    minwt = np.min(np.abs(dz / solver.state['w']['g'].real))
    dt = safety * min(minut, minwt)

    if domain.distributor.comm_world:
        dt = domain.distributor.comm_world.gather(dt, root=0)
        if domain.distributor.rank == 0:
            dt = [min(dt)] * domain.distributor.size
        dt = domain.distributor.comm_world.scatter(dt, root=0)

    return dt

CFL = 0.2

dt_min = 1e-2

solver.start_time=0.
solver.iteration=0
solver.stop_iteration=400
solver.sim_stop_time = np.inf
solver.wall_stop_time = np.inf
solver.dt = dt_min
dt_cadence = 5
copy_cadence = 20

analysis1 = solver.evaluator.add_file_handler("KE_integral_32", iter=10)
analysis1.add_task("Integrate(0.5*u*u+0.5*w*w)")
#analysis1.add_task("Integrate(0.5*u*u+0.5*w*w)")


#analysis1 = solver.evaluator.add_file_handler("restart", sim_dt=5.0)
#analysis1.add_system(solver.state)

#analysis2 = solver.evaluator.add_file_handler("density", sim_dt=0.05)
#analysis2.add_task("Integrate(s, dz)", name="density profile")
#analysis2.add_task(s, name="density field")

# Main loop
start_time = time.time()
while solver.ok:

  # adaptive timestep
  if solver.iteration % dt_cadence == 0:
    dt = min(dt_min,cfl_dt(CFL))

  # advance
  solver.step(dt)

  # update lists
  if solver.iteration % copy_cadence == 0:
    logger.info('Iteration: %i, Time: %e' %(solver.iteration, solver.time))

end_time = time.time()

# Print statistics
logger.info('-' * 20)
logger.info('Total time: %e' %(end_time - start_time))
logger.info('Iterations: %e' %(solver.iteration))
logger.info('Average timestep: %e' %(solver.time / solver.iteration))
logger.info('-' * 20)


