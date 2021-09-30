import numpy as np
import scipy.special as sp
from mpi4py import MPI
import time
import os
import math
from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger=logging.getLogger(__name__)

# Parameters
L=16. # Domain length
H=1.  # Domain height
h=H/2

#Re=1000. give a Reynolds number based on the front velocity and current height of around 420 

Re=8000.
Nx = 512 # x fourier modes

#Re=2000.
#Nx = 256# x fourier modes
Ny=round(2*H/L*Nx)
Sc=Re
q=0*math.pi/180

#U=-0.5  # 13
#U=-0.7  # 14 delayed velocity. Steady nose at 6.5 t~200
#U=-0.6  # 16 delayed velocity. Nose slowly moving past 10 at t~200
#U=-0.45 # 23
#U=-0.49 # 23
#U=-0.5  # 22

#Earlier had inconsistent velocity enforcement
#U=-0.50  # 31   w2=2      
#U=-0.55 # 30   w2=1
#U=-0.60 # 29   w2=1
#U=-0.65 # 28   w2=1
#U=-0.70 # 27   w2=1
#U=-0.75 # 32   w2=2 and upper noundary in b fixed
#U=-0.80  # 33   w2=2 and upper noundary in b fixed
#U=-0.80  # 34   Fixed boundaries so as to be smooth
#U=-0.75  # 35 
#U=-0.7   # 36 
#U=-0.85  # 37
U=-0.65   # 38 
U=-0.65   # 39 
U=-0.60   # 40
U=-0.70   # 41
U=-0.75   # 42
U=-0.60   # 43
# New conditions with 0 bouyancy enforced
U=-0.40   # 47 # Head stable at 16
U=-0.45   # 46 # Head stable at 16
U=-0.50   # 45 # Head stable at 16
U=-0.51   # 51 # head at 14 still advancing
U=-0.52   # 50 # head at 12 still advancing
U=-0.53   # 49 # stable head at 10 
U=-0.54   # 48 # stable head at 8.3
U=-0.55   # 44 # stable head at 7.5 but then collapses
U=-0.515  # 52 # incomplete

U=-0.515  # 53 linear initial velocity stable head at 15
U=-0.51   # 54 linear initial velocity stable head at 15
U=-0.52   # 55 linear initial velocity stable head at 14
U=-0.525  # 56 linear initial velocity stable head at 13
U=-0.53   # 57 # stable head at 10 
U=-0.54   # 58 # stable head at 8.3
U=-0.55   # 59 # stable head at 7.5 but then collapses

# Changed forcng velocity regind to 1.5 width
# changed forcing density region chnged width to 0.5

U=-0.51   # 60 
U=-0.52   # 61 
U=-0.53   # 62 
U=-0.54   # 63 
U=-0.55   # 64 

#No-slip bottom boundary
U=-0.50   # 65
U=-0.45   # 66
U=-0.40   # 67
U=-0.35   # 68

param.bc.lb='noslip'
param.bc.ub='slip'

U1 = 0; #h*U1+(H-h)*U2=H*U
U2 = (H*U-h*U1)/(H-h);

#if Re>1000:
#    Nx = 512 # x fourier modes
#if Re>2000:
#    Nx = 1024 # x fourier modes

gy=-math.cos(q)
gx= math.sin(q)

# Create bases and domain
start_init_time=time.time()
x_basis=de.Fourier(   'x', Nx, interval=(0, L), dealias=3/2)
y_basis=de.Chebyshev('y', Ny, interval=(0, H), dealias=3/2)
domain=de.Domain([x_basis, y_basis], grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','b','u','v','by','uy','vy'], time='t')

Dx=Nx/L
Dy=Ny/H


problem.parameters['nu']    = 1/Re
problem.parameters['kappa'] = 1/Sc
problem.parameters['L']    = L
problem.parameters['H']    = H
problem.parameters['gx']    = gx
problem.parameters['gy']    = gy
problem.parameters['U']     = U


# Get pointers to the grid coordinates
x = domain.grid(0)
y = domain.grid(1)

mdt=0.01
F=1/mdt

Wx=4/Dx;
Wy=4/Dx;

# Make the weight/time constants smooth
# The underlying target function can be sharp
def indicator(x,x1,x2,D):
    if D==0:
        y = (np.sign(x-x1)-np.sign(x-x2))/2
        y = np.round(y)
    else:
        y = (sp.erf((x-x1)*D)-sp.erf((x-x2)*D))/2
    return(y)

def heaviside(x,D):
    if D==0:
        y = (1+np.sign(x))/2
        y = np.round(y)
    else:
        y = (1+sp.erf(x*D))/2
    return(y)

x0=Wx;
x1=1.5+Wx;

x2=3
x3=4
x4=4.5
x5=x4+5 # initial condition 
w1=indicator(x,x0,   x1,  Dx)   # This region force uniform velocity profile
W1=indicator(x,x0-Wx,x1+Wx,0)  

w2=indicator(x,x2,   x3,  Dx)   # Force zero velocity and high and low velocity
W2=indicator(x,x2-Wx,x3+Wx,0)

w3=indicator(x,x3,   x4,  Dx)   # Buoyancy region
W3=indicator(x,x3-Wx,x4+Wx,0) 

H1=heaviside(h-y,0)
h1=heaviside(h-y,Dy)
H2=heaviside(y-(H-h),0)
h2=heaviside(y-(H-h),Dy)

nfb = domain.new_field(name='nfb')  # Target b
nfu = domain.new_field(name='nfu')  # Target u
nwb = domain.new_field(name='nwb')  # Weight for b
nwu = domain.new_field(name='nwu')  # Weight for u
nwv = domain.new_field(name='nwv')  # Weight for v
# Target v is zero

# Buoyancy 0 in regions 1 and 2 1 in region 3
I=np.maximum(0,np.minimum(1,(indicator(x,x0,x3,Dx)+indicator(x,x0,x4,Dx)*H1)))

#nwb['g'] = F*I
#nfb['g'] = indicator(x,x3,x4+Wx,0)*heaviside(Wy+h-y,0)
#nfb['g'] = (sp.erf((x-x3)*Dx)-np.sign(x-x4-Wx))/2*heaviside(Wy+h-y,0)
nwb['g']  = F*indicator(x,x0,   x4,  Dx)
nfb['g']  = (sp.erf((x-x3)*Dx)-np.sign(x-x4-Wx))/2*h1


nwv['g']  = F*(w1 + w2)
nwu['g']  = F*(w1 + w2)
nfu['g']  = U*(W1 + W2*H/(H-h)*H2)

problem.parameters['fb'] = nfb
problem.parameters['fu'] = nfu
problem.parameters['wb'] = nwb
problem.parameters['wu'] = nwu
problem.parameters['wv'] = nwv

print("Runing ded_gc_primtive.py: Boussinesq lock-exchange gravity current in primitive formulation with odd symmetry")
print("L={}, H={}, Nx={}, Ny={}, Re={}, Sc={}, U={}".format(L,H,Nx,Ny,Re,Sc,U))

problem.add_equation("dx(u) + vy = 0")
problem.add_equation("dt(b) - kappa*(dx(dx(b)) + dy(by))                =  wb*(fb-b) - u*dx(b) - v*by")
problem.add_equation("dt(u) -    nu*(dx(dx(u)) + dy(uy)) + dx(p) - gx*b =  wu*(fu-u) - u*dx(u) - v*uy")
problem.add_equation("dt(v) -    nu*(dx(dx(v)) + dy(vy)) + dy(p) - gy*b = -wv*v      - u*dx(v) - v*vy")
problem.add_equation("by - dy(b) = 0")
problem.add_equation("uy - dy(u) = 0")
problem.add_equation("vy - dy(v) = 0")

problem.add_bc("left(by)   = 0")
problem.add_bc("right(by)  = 0")

if param.bc.lb == "noslip":
    problem.add_bc("left(u)  = 0")
elif param.bc.lb == "slip":
    problem.add_bc("left(u)  = U")
else:
    logger.error('Unknown param.bc.lb "' + param.bc.lb + '"')
    
param.bc.ub='slip'

problem.add_bc("left(u)    = U")
problem.add_bc("right(uy)  = 0")

problem.add_bc("left(v)    = 0", condition="(nx != 0)")
problem.add_bc("right(v)   = 0")
problem.add_bc("integ_y(p) = 0", condition="(nx == 0)")


# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Get pointers to the simulation variables and set them to zero
b=solver.state['b']
u=solver.state['u']
v=solver.state['v']
p=solver.state['p']
by=solver.state['by']
uy=solver.state['uy']
vy=solver.state['vy']

# Set the buoyancy field to lock exchange with a slightly smoothed edge to prevent ringing
b['g'] = indicator(x,x3,x5,Dx)*h1
b.differentiate('y', out=by)

u['g'] = indicator(x,x3,x5,Dx)*U*H/(H-h)*np.minimum(1,2*y/h-1)
u.differentiate('y', out=uy)

# Integration parameters
solver.stop_sim_time = 300
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

# Read the output directory from the environment variable DEDALUS_DATA
data_dir=os.environ['DEDALUS_DATA']+'/gc'
print("Writing data to directory '{}'".format(data_dir))

# Set the output data to writing buoyancy every 0.1 time units
# scale increases the resolution of the output datafiles
snap = solver.evaluator.add_file_handler(data_dir, sim_dt=1)
snap.add_task("b",               scales=1, name='b')
snap.add_task("u",               scales=1, name='u')
snap.add_task("v",               scales=1, name='v')
snap.add_task("p",               scales=1, name='p')
snap.add_task("integ(b,   'y')", scales=1, name='B')
snap.add_task("integ(u*b, 'y')", scales=1, name='Q')
snap.add_task("integ(p,   'y')", scales=1, name='P')
snap.add_task("integ(u*u, 'y')", scales=1, name='UU')
snap.add_task("integ(u*v, 'y')", scales=1, name='UV')
snap.add_task("integ(v*v, 'y')", scales=1, name='VV')
snap.add_task("integ(b, 'x','y')",     scales=1, name='B00')
snap.add_task("integ(b*x, 'x','y')",   scales=1, name='B10')
snap.add_task("integ(b*y, 'x','y')",   scales=1, name='B01')
snap.add_task("integ(b*x*x, 'x','y')", scales=1, name='B20')
snap.add_task("integ(b*x*y, 'x','y')", scales=1, name='B11')
snap.add_task("integ(b*y*y, 'x','y')", scales=1, name='B02')

# CFL
CFL = flow_tools.CFL(solver, initial_dt=1e-2, cadence=5, safety=0.8, max_change=1.5, min_change=0.5, max_dt=mdt)
CFL.add_velocities(('u','v')) # maybe calculate with w*sin(y/H) since in y direction the grid-spacing is non uniform

# Flow properties Print some information about the simulation every cadence iterations
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("sqrt(u*u + v*v)*H/nu/2", name='Re')

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



