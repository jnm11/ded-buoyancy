#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Burgers' equation
"""
Burgers' equation

Usage: 
   ded_gc.py [options] NAME

Options:
  -h, --help            : Show this help message
      --sType=<>        : Type of simulation
      --preset=<>       : Run a simulation with preset parameters
      --scheme=<>       : Integration scheme
  -R, --Re=<>           : Reynolds number  
  -F, --force=<>        : Strength of forcing
      --clipu=<>        : Clip velocities to less than this absolute value
      --clipB=<>        : Clip bouyany between 0 and B
      --maxdt=<>        : Maximum time step
      --mindt=<>        : Minimum time step
      --Sc=<>           : Schmidt number   
  -H, --height=<>       : Box height       
  -W, --width=<>        : Box width        
  -L, --length=<>       : Box length       
      --radius=<>       : Inlet radius      
      --hb=<>           : Depth of buoyancy forcing region
      --hu=<>           : Depth of velocity forcing region
      --alpha=<>        : Entrainment coefficient
  -U, --velocity=<>     : Flow velocity    
      --U1=<>           : Current velocity 
  -N, --Nx=<>           : x grid points    
      --Ny=<>           : y grid points   
      --Nz=<>           : z grid points   
      --Tx=<>           : type of x grid
      --Ty=<>           : type of y grid
      --Tz=<>           : type of z grid
  -T, --Time=<>         : Simulation time  
  -f, --forcing=<>      : Forcing type     
      --xa=<>           : Start of region to fill with material
      --xb=<>           : End of region to fill with material
      --x0=<>           : Start of velocity uniform region
      --x1=<>           : End of velocity uniform region
      --x2=<>           : Start of density region
      --x3=<>           : end of density region
      --x4=<>           : 
      --x5=<>           : 
      --x6=<>           : 
      --x7=<>           : 
      --xn1=<>          : start noise forcing region
      --xn2=<>          : 
      --xn3=<>          : 
      --xn4=<>          : end noise forcing region
      --Wx=<>           : x transition width
      --Wy=<>           : y transition width
      --Wz=<>           : z transition width
  -c, --ck_time=<>      : Checkpoint time minutes       
  -A, --angle=<>        : Slope angle (degrees)    
  -g, --gravity=<>      : Gravity
  -l, --lbc=<>          : Lower boundary condition 
  -u, --ubc=<>          : Upper boundary condition 
  -t, --time=<>         : Max runtime 0 means infinite 
      --signals=<>      : Catch signals      
      --parallel=<>     : Use parallel HDF   
  -s, --start_new_files : Start new files while checkpointing 
  -r, --restart         : Restart  
      --rfn=<>          : Restart file name
  -p, --param           : Load from parameter file
      --pfn=<>          : Parameter file name              
      --reset           : reset solve sim time and integrated variables
  -B, --buoyancy        : Buoyancy
  -S, --sediment        : Sediment
      --SV              : Sediment velocity
  -n, --noise=<>        : amplitude    for noise forcing
      --noiseL=<>       : length scale for noise forcing
      --noised=<>       : xy, xz, or yz plane to impose noise forcing from stream function
      --noiseT=<>       : time scale for noise forcing
      --inoise=<>       : initial noise
  -i, --inlet=<>        : inlet type gaussian, circle, square, 
      --AA=<>           : anti-aliasing default is 3/2
      --AB=<>           : scale factor for 3d fields
      --VMEM=<>         : set maximum VMEM in GB (=0 unlimited)  
      --dtb=<>          : time frequency for writing out 3d density data
      --dtu=<>          : time frequency for writing out 3d u velocity data
      --dtv=<>          : time frequency for writing out 3d v velocity data
      --dtw=<>          : time frequency for writing out 3d w velocity data
      --dtp=<>          : time frequency for writing out 3d p velocity data
      --dtstats=<>      : time frequency for writing out statistics
      --dtnoise=<>      : time frequency for writing out 3d noise data
      --dtx=<>          : time frequency for writing out x integrated data
      --dty=<>          : time frequency for writing out y integrated data
      --dtz=<>          : time frequency for writing out z integrated data
      --dtxy=<>         : time frequency for writing out xy integrated data
      --dtxz=<>         : time frequency for writing out xz integrated data
      --dtyz=<>         : time frequency for writing out yz integrated data
      --dtxyz=<>        : time frequency for writing out xyz integrated data
      --dtmomb=<>       : time frequency for writing out moments of b
      --dtavrg=<>       : time frequency for writing out time averaged data
      --dtleft=<>       : time frequency for writing out left boundary data
      --dtright=<>      : time frequency for writing out right boundary data
      --dtpm=<>         : time frequency for writing out plume data
      --dtgc=<>         : time frequency for writing out gravity current data
      --dtslice=<>      : time frequency for outputting slice data
      --dtforce=<>      : time frequency for outputting forcing data
      --dtfpid=<>       : time frequency for outputting pid forcing data
      --fpidu=<>        : integral error control for u velocity
      --fpidv=<>        : integral error control for v velocity
      --fpidw=<>        : integral error control for w velocity
      --fpidb=<>        : integral error control for b
      --fpids=<>        : integral error control for s
      --fpid=<>         : integral error control scale factor
      --slicex=<>       : comma seperated list of x locations to output data
      --slicey=<>       : comma seperated list of y locations to output data
      --slicez=<>       : comma seperated list of z locations to output data
      --PIDG=<>         : PID controller sets gravity not velocity
      --PIDD=<>         : PID Set Derivative Coefficient
      --PIDI=<>         : PID Set Integral Coefficient 
      --PIDP=<>         : PID Set Proportional Coefficient           
      --PIDT=<>         : PID Computation Time Interval
      --PIDX=<>         : PID set target X value
      --PIDIT=<>        : PID integration error term
      --PIDDD=<>        : PID derivative error term
      --PIDST=<>        : PID square wave time interval
      --PIDS1=<>        : PID square X position 1 
      --PIDS2=<>        : PID square X position 2
"""

import time
import numpy as np
timea = [np.float(time.time())]
ver=np.float(1.32)
import warnings
warnings.filterwarnings('ignore',category=FutureWarning) 
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
from docopt import docopt
import scipy.special as sp
from mpi4py import MPI
import resource
import math
import h5dict
import jutil as ju
import ded_analysis as da
import ded_types as dp
import jpid
import h5py
import copy
import jtheta
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools  import post
from dedalus.core   import operators
from polytrope.tools.checkpointing import Checkpoint
from pathlib import Path

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

import logging
ju.logger=logging.getLogger(__name__)

np.seterr(divide='ignore', invalid='ignore') # Ignore /0 and 0/0

p,name = dp.name_type(args)
sType=p['sType']

bdir,ddir=ju.set_base_dir(sType,name)
ju.logger.info("name:            {}".format(name))
ju.logger.info("sType:           {}".format(sType))
ju.logger.info("bdir:            {}".format(bdir))
ju.logger.info("ddir:            {}".format(ddir))

ju.check_version(ver)

if args['--pfn']: pfn = Path(args['--pfn'])
else:             pfn = None

if args['--rfn']: rfn = Path(args['--rfn'])
else:             rfn = None
    
if args['--param'] or pfn:
    if not pfn: pfn=Path('final')
    fpfn=ju.find_param_file(bdir,name,pfn)
    param=ju.read_param(fpfn)
    p=dp.merge_param(p,param)

if args['--preset']: p=dp.preset(p,args['--preset'])
p=dp.read_args(p,args)

p=dp.validate(p)
ju.Barrier()
ju.logger.info('Parameters validated')

wall_time       = '99:00:00'
start_new_files = False
restart         = False

if args['--time']:                     wall_time =  np.str(args['--time'])
if args['--start_new_files']:          start_new_files = True
if args['--restart'] or args['--rfn']: restart = True
else:                                  restart = False
if args['--reset']:                    reset = True
else:                                  reset = False

ju.write_status("Running")
ju.write_files()

scheme   = p['scheme']
Tx       = p['Tx']
Nx       = p['Nx']
L        = p['L']
hu       = p['hu']
hb       = p['hb']
alpha    = p['alpha']
Wx       = p['Wx']
AA       = p['AA']
q        = p['q']
g        = p['g']
Re       = p['Re']
U        = p['U']
T        = p['T']
U1       = p['U1']

noise    = p['noise']
noiseL   = p['noiseL']
noised   = p['noised']
noiseT   = p['noiseT']
xn1      = p['xn1']
xn2      = p['xn2']
xn3      = p['xn3']
xn4      = p['xn4']

x0       = p['x0']
x1       = p['x1']
x2       = p['x2']
x7       = p['x7']
x3       = p['x3']
x4       = p['x4']
x5       = p['x5']
x6       = p['x6']
x7       = p['x7']

# Don't save parameters that setup initial state
inoise   = p.pop('inoise', None)
xa       = p.pop('xa',     None)
xb       = p.pop('xb',     None)

VMEM     = p['VMEM']
inlet    = p['inlet']
sType    = p['sType']
signals  = p['signals']
parallel = p['parallel']
forcing  = p['forcing']
ck_time  = p['ck_time']
lbc      = p['lbc']
ubc      = p['ubc']

dtstats  = p['dtstats']

PIDG     = p['PIDG']
PIDD     = p['PIDD']
PIDI     = p['PIDI']
PIDP     = p['PIDP']
PIDT     = p['PIDT']
PIDX     = p['PIDX']
PIDIT    = p['PIDIT']
PIDDD    = p['PIDDD']
PIDST    = p['PIDST']
PIDS1    = p['PIDS1']
PIDS2    = p['PIDS2']
maxdt    = p['maxdt']
mindt    = p['mindt']
F        = p['force']
clipu    = p['clipu']
clipB    = p['clipB']

Nxx       = np.int(np.round(p['Nx']*p['AA']))

dt = None

if maxdt*F>2:
    maxdt=2/F
    ju.logger.info('Limiting maxdt due to forcing to {:8.3e}'.format(maxdt))
    
ju.make_lock()
if signals: ju.init_signals()

wt=ju.get_sec(wall_time)

p=dp.remove_None(p)
ju.logger.info("Writing parameters")
ju.write_param(p)  

dx=L/Nx
dV=dx
    

# Should add a test to see if h5py supports parallel
ck_parallel=parallel
fh_parallel=parallel
max_writes=5
#https://stackoverflow.com/questions/47072859/how-to-append-data-to-one-specific-dataset-in-a-hdf5-file-with-h5py

ju.open_stats('stats.hdf5',dtstats,reset)

     
    # 3D Boussinesq hydrodynamics

#print("Tz: {}, Nx: {}, Ny: {}, Nz: {}".format(Tz,Nx,Ny,Nz))
#print("L: {}, W: {}, H: {}".format(L,W,H))

domain,problem,xx,yy,zz=ju.boussinseq_init(p)


if (noise>0) and ( (sType=='pm' and not forcing==1) or (sType=='gc') ):
    if noised=='xy':
        problem.substitutions['noisex'] = ' wnoise*dy(noisef2)'
        problem.substitutions['noisey'] = '-wnoise*dx(noisef2)'
        problem.substitutions['noisez'] = '0'
    elif noised=='xz':
        problem.substitutions['noisex'] = ' wnoise*dz(noisef2)'
        problem.substitutions['noisey'] = '0'
        problem.substitutions['noisez'] = '-wnoise*dx(noisef2)'
    else:         
        problem.substitutions['noisex'] = '0'
        problem.substitutions['noisey'] = ' wnoise*dz(noisef2)'
        problem.substitutions['noisez'] = '-wnoise*dy(noisef2)'

else:
    problem.substitutions['noisex'] = '0'
    problem.substitutions['noisey'] = '0'
    problem.substitutions['noisez'] = '0'
     
problem.substitutions['gy']    = '0'
problem.substitutions['SVY']   = '0'
if   q==  0:
    problem.substitutions['gx']    = '0'
    problem.substitutions['gz']    = '-1'
    problem.substitutions['SVX']   = '0'
    problem.parameters['SVZ']      = -SV
elif q==  90: 
    problem.substitutions['gx']    = '1'
    problem.substitutions['gz']    = '0'
    problem.parameters['SVX']      = SV
    problem.substitutions['SVZ']   = '0'
elif q==-180: 
    problem.substitutions['gx']    = '0'
    problem.substitutions['gz']    = '1'
    problem.substitutions['SVX']   = '0'
    problem.parameters['SVZ']      = SV
else:
    gx = math.sin(math.pi*q/180)
    gz =-math.cos(math.pi*q/180)
    problem.parameters['gx']       = gx
    problem.parameters['gz']       = gz
    problem.parameters['SVX']      = gx*SV
    problem.parameters['SVZ']      = gz*SV

    
problem.parameters['nu']    = 1/Re
problem.parameters['fnu']   = 0/Re/1000
problem.parameters['kappa'] = 1/(Re*Sc)
problem.parameters['L']     = L
problem.parameters['W']     = W
problem.parameters['H']     = H
problem.parameters['U']     = U
problem.parameters['B']     = B
problem.parameters['S']     = S
problem.parameters['noise']  = noise
problem.parameters['noiseT'] = noiseT
problem.parameters['noiseL'] = noiseL
problem.parameters['Volume'] = W*H*L

#problem.meta[:][‘z’][‘dirichlet’] = True
#if PIDX==0 or PIDG: problem.meta['U']['x','y']['constant']

# Set the length for Reynolds number calculation
if sType=="gc": RL = hu
else:           RL = 2*R 

# Get pointers to the grid coordinates
#print("domain.grid(0) {} ".format(dir(domain.grid(0))))
#print("domain.dist.grid_layout.local_shape {} ".format(dir(domain.dist.grid_layout)))
#print("domain.dist.grid_layout.local_shape {} ".format(domain.dist.grid_layout.local_shape(scales=1)))
#domain.dist.grid_layout.local_shape(scales=AA)
#print("domain.dist.grid_layout.local_shape {}".format(domain.dist.grid_layout.local_shape(scales=AA)))
#print('xx.shape {}'.format(xx.shape));
#domain.dist.grid_layout(scales=1)
#print("domain.dist.grid_layout".format(domain.dist.grid_layout))
x = domain.grid(0)  # These will contain only the local grid coordinates
if Ny>1:
    y = domain.grid(1)
    z = domain.grid(2)
else:
    y = 0
    z = domain.grid(1)

#rs=(             np.reshape(np.tile((xx-L/2)**2,len(yy)),[len(xx),len(yy)])
#   +np.transpose(np.reshape(np.tile((yy-W/2)**2,len(xx)),[len(xx), len(yy)])))/r**2

#print("x.shape: {}".format(x.shape));
#print("y.shape: {}".format(y.shape));
#print("z.shape: {}".format(z.shape));
#print("Tz:      {}".format(Tz));
#print("y: {}".format(y))
#print("z: {}".format(z))

# Make the weight/time constants smooth
# The underlying target function can be sharp

Tdiv=None
rand = np.random.RandomState()#(seed=42)

# return a random field scaled by sqrt of time step in grid space
# the mean vaue doesn't matter since we are treating as  stream function and we filter out high frequencies
# see internal_heat.py

#nfx = domain.new_field()
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
if Tz=='Cheb':
    #N=100;dz=sin(pi*(1+2*(0:N-1))/(2*N))*sin(pi/(2*N));plot(max(dz,max(dz)/10));
    dz1=H*np.sin(math.pi*(1+2*np.arange(Nz ,dtype=float))/(2*Nz ))*np.sin(math.pi/(2*Nz ))
    dz2=H*np.sin(math.pi*(1+2*np.arange(Nzz,dtype=float))/(2*Nzz))*np.sin(math.pi/(2*Nzz))
    #z1=np.maximum(dz1,max(dz1)/10)
    #dz2=np.maximum(dz2,max(dz2)/10) # Don't make the noise too strong near the walls
    #    print('sum(dz1) {} '.format(sum(dz1)))
    #    print('sum(dz2) {} '.format(sum(dz2)))
else:
    dz1=H/Nz
    dz2=H/Nzz



if Ny>1: dA=dz2*L/Nx*W/Ny
else:    dA=dz2*L/Nx

Edy=dy/W
Edz=dz2/H

def anoise(sc):
    gshape = domain.dist.grid_layout.global_shape(scales=sc)
    slices = domain.dist.grid_layout.slices(scales=sc)
    x = rand.standard_normal(gshape)[slices]
    return(x)

# Code never divides up along Cheb directions so don't need to do anything special for dz1

#print('anoise.shape {}'.format(anoise(AA).shape))
#print('dz2.shape {}'.format(dz2.shape))
#print('gshape {}'.format(gshape))
#print('slices {}'.format(slices))
#print('slices[1] {}'.format(slices[1]))
def gnoise(deltaT):
    if Ny>1:
        x = np.divide(anoise(AA),np.sqrt(dA))
#        print('max x {}'.format(a.max(dim=)))
    else:
        slices = domain.dist.grid_layout.slices(scales=AA)
        x = np.divide(anoise(AA),np.sqrt(*dA[slices[1]]))
    return(x)

zmn = operators.GeneralFunction(domain,'g',gnoise,args=[])
zmn.original_args = [0.0001]
#class ClassForce(operators.GeneralFunction):
#    def meta_parity(self, axis):
#        #parity0 = self.args[0].meta[axis]['parity']
#        return(1)

##znm = ClassForce(domain,'g', gnoise, args=[])
#zmn.original_args = [0.0001]

#        if axis == 'z':
#                            return (-1) * parity0
#            else:
#                return parity0

fb = None
fs = None
fu = None
fv = None
fw = None
fd = None

if PIDX>0 and PIDG:
    nfg = domain.new_field()  # gravity function
    nfg['g'] = g
    if   p['Tx'] == "SinCos": nfg.meta['x']['parity'] = -1
    else:                     nfg.meta['x']['parity'] = 1
    if Ny>1: nfg.meta['y']['parity'] = 1
    nfg.meta['z']['parity'] = 1
    problem.parameters['g'] = nfg
else: problem.parameters['g'] = g

U2 = None
wnoise = None
nwb = None
nws = None
nwu = None
nwv = None
nww = None

nfb = None
nfs = None
nfu = None
nfv = None
nfw = None
nfd = None

ifb = None
ifs = None
ifu = None
ifv = None
ifw = None


WN=None
WB=None
WS=None
WU=None
WV=None
WW=None

EN=None
 
if sType=="pm":
    r = np.sqrt(y*y+z*z)
    IA = math.pi*R**2 # Target area
    MEANU = U*IA/(H*W)
    problem.parameters['meanu']      = MEANU
    problem.substitutions['meanv']   = '0'
    problem.substitutions['meanw']   = '0'

    if inlet == "gaussian":
        wrmax=2
        wr = 2*np.exp(-2*(r/R)**2)             # Make maximum and area pi*R^2 and int(f^2)=int(f)
    elif inlet == "circle":
        wrmax = 1
        wr = ju.cubici(r,R+Wy/2,R-Wy/2)
    QV=ju.jsum(dy*dz*U*wr)/IA
    QM=ju.jsum(dy*dz*U*U*wr*wr)/IA
    QB=ju.jsum(dy*dz*B*wr*wr)/IA
    QQ=ju.jsum(dy*dz*U*B*wr*wr)/IA

    if noise>0:
        wnoise = domain.new_field()  # weighting for noise
        wnoise['g']=(ju.cubici( x,  xn1, xn2) - ju.cubici( x, xn3, xn4))*wr/wrmax

    if forcing==5:
        # x0 = 0.5 [0.0 0.5] negative div(u)
        # x1 = 1.0 [0.5 1.0] zero velocity
        # x2 = 1.5 [1.0 1.5] positive div(u)
        # x3 = 2.0 [1.5 2.0] enforced u,v,w        
        
        nfb = domain.new_field() # Target b
        nfu = domain.new_field() # Target u
        nfd = domain.new_field() # Divergence
        nwb = domain.new_field() # Weight for b
        nwu = domain.new_field() # Velocity shaping
    
        # This code is only correct for MPI parallelisation in the y & z directions not x
        
        
        #  0 x1 B=0
        # x1 x3 B=1
        #  0 x0 div < 0 
        # x0 x1 u=v=w=0
        # x1 x2 div > 0 v==0 u = int div  
        # x2 x3 u=1,v=0,w=0
        #
        x4=2*x3-x2
        dx1 = ju.dcubici(x, x1, x2)
        dx2 = ju.dcubici(x,  0, x0)

        dx1 = dx1/(dx*dx1.sum(0))
        dx2 = dx2/(dx*dx2.sum(0))

        fxb = ju.cubici( x, x0, x2) - ju.cubici( x, x3, x4)  
        fx1 = ju.cubici( x, x1, x2) - ju.cubici( x, x3, x4)

        fd1=wr*dx1
        fd2=np.ones(wr.shape)*dx2

        TQ1 = dV*ju.jsum(fd1)
        TQ2 = dV*ju.jsum(fd2)
        
        Q=  TQ1
        U2= -Q/TQ2
        
        nwb['g'] = ju.cubici( x,  0, x0) - ju.cubici( x, x2, x3)   
        nwu['g'] = ju.cubici( x, x0, x1) - ju.cubici( x, x2, x3)   
        nfb['g'] = B*fxb*wr
        nfu['g'] = fx1*wr 

        nfd['g'] = fd1     + U2*fd2
        Tdiv = ju.jint(nfd)

        TQ1 = dV*ju.jsum(fd1)/IA
        TQ2 = dV*U2*ju.jsum(fd2)/IA
        TQ3 = dV*ju.jsum(nfd['g'])/IA
        
        if Ny>1: problem.substitutions['wv']   = 'wu'
        problem.substitutions['ww']   = 'wu'
    elif forcing==1:
        problem.substitutions['fd']     = '0'
        problem.substitutions['fb']     = '0'
        problem.substitutions['fu']     = '0'
        if Ny>1: problem.substitutions['fv']     = '0'
        problem.substitutions['fw']     = '0'
        problem.substitutions['wb']     = '0'
        problem.substitutions['wu']     = '0'
        if Ny>1: problem.substitutions['wv']     = '0'
        problem.substitutions['ww']     = '0'
        Tdiv = 0
        TB   = B*IA
        Q    = IA
        ifw  = IA/(L*W)
    elif forcing==2:  # stream function forcing between two gaussians
        # periodic-gaussian.mw
        # x0 = 0.5 [0.0 0.5] gaussian 1
        # x1 = 1.0 [0.5 1.0] switch region
        # x2 = 1.5 [1.0 1.5] switch region
        # x3 = 2.0 [1.5 2.0] gaussian 2         

        nfb = domain.new_field() # Target b
        nfu = domain.new_field() # Target u
        if Ny>1: nfv = domain.new_field() # Target v
        nfw = domain.new_field() # Target w
        nwb = domain.new_field() # Weight for b

        #      0
        # x0 = 0.5
        # x1 = 1.0
        # x2 = 1.5
        # x3 = 2.0
        nfx=2;
        nfz=np.int(np.round(Nz*R/(2*H)));
        fx1,Ffx1,dfx1=ju.smoothi(x,nfx,x0,x2)
        fx2,Ffx2,dfx2=ju.smoothi(x,nfx,x3,L)
        fx=fx1-fx2
        bb=0.5
        dfx=(dfx1-dfx2)/R*(1+bb)/(1+bb*fx)**2
        #ib=(1+bb)*fx/(1+bb*fx)/R
        ib=fx
 
 
        U1=1;
        U2=math.pi*R**2/(H*W)
        
        nfb['g'] = B*fx*wr
        nwb['g'] = (ju.cubici( x, 0, x0) - ju.cubici( x, x2, x3))

   
        #psi = domain.new_field() # Stream function
        WS=np.array(ib.data,dtype=float)*W/np.sqrt(2.0)
        HS=np.array(ib.data,dtype=float)*H/np.sqrt(2.0)
        
        E1=sp.erfc(WS)
        E2=sp.erfc(HS)
        EE=math.pi*(E1*E2-E1-E2)/(H*W)
        wb =   np.exp(-2*r**2*ib**2)
        #psi['g'] = -U1*R**2*(wb+EE*r**2)/2;
        
        fux      =  U1*R**2*    ( 2*ib**2*wb-EE)
        fur      =  U1*R**2*dfx*(-2*ib*wb + np.sqrt(math.pi/2)*(sp.erf(WS)*np.exp(-HS**2)/W + sp.erf(HS)*np.exp(-WS**2)/H))
        
        nfu['g'] = fux
        nfv['g'] = fur*y
        nfw['g'] = fur*z
        
        problem.substitutions['wu'] = 'wb'
        problem.substitutions['wv'] = 'wb'
        problem.substitutions['ww'] = 'wb'

        del(fux,fur,fx)
 
    elif forcing==3 or forcing==4:  # stream function forcing between two gaussians with weighting
        # periodic-gaussian.mw
        # Velocity is forced netween x0 and x3
        # Velocity transition occurs between x1 and x2
        # Buoyancy is forced between x4 and x7
        # Buoyancy transition occures beteen x5 and x6

        nfb = domain.new_field() # Target b
        nfu = domain.new_field() # Target u
        nfv = domain.new_field() # Target v
        nfw = domain.new_field() # Target w
        nwb = domain.new_field() # Weight for b
        nwu = domain.new_field() # Weight for u

        #      0
        # [x0,x1] velocity transition region [x2,x3] detransition to ensure periodicity
        # [0 x2 x3 x4] foring region
        
        # Velocity is forced between x0 and x3
        # Transition occur between x1 and x2
        nwu['g'] = ju.cubici( x,x0,x0+Wx)-ju.cubici( x,x3-Wx,x3)
        U2 = math.pi*R**2/(H*W)
        c0 = np.sqrt(2)/R

        if x2>x1:
            fx       = ju.cubici( x,x1,x2)-ju.cubici( x,x3,x3+(x2-x1))
            dfx      = ju.dcubici(x,x1,x2)-ju.dcubici(x,x3,x3+(x2-x1))
            if alpha>0:   # Estimate the width of the Gaussian at the top
                R1 = R + alpha*(L-x3)  
                c1 = np.sqrt(2)/R1
            else: c1=0    # Make a uniform outflow condition
            c  =  fx*(c0-c1)+c1
            dc = dfx*(c0-c1)
        else: # velocity forcing is indepedendent of x
            c=c0
            dc=0
        
        #Buoyancy is forced between 0 and x7
        #B   0     0<x<x5
        #    0-1  x5<x<x6
        #    1    x6<x<x7
        nfb['g'] = B*wr*(ju.cubici( x, x5, x6)    - ju.cubici( x, x7,    x7+(x6-x5)))
        nwb['g'] =      (ju.cubici( x, x4, x4+Wx) - ju.cubici( x, x7-Wx, x7))

        nn=5
        if forcing==3:   nfu['g'],nfv['g'],nfw['g'] = ju.weighted_gaussian(c,dc,y,z,h,1,W,H,nn)
        else:            nfu['g'],nfv['g'],nfw['g'] = jtheta.periodic_gaussian(c,dc,y,z,1,W,H,R)
        
        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'

        if forcing==3: fu,fv,fw = ju.weighted_gaussian(c0,0,y,z,h,1,W,H,nn) # Forcing function at radius
        else:          fu,fv,fw = jtheta.periodic_gaussian(c0,0,y,z,1,W,H,R)
 
        ifu=nfu['g']
        ifv=nfv['g']
        ifw=nfw['g']


        del(fu,fv,fw,fx,dfx,c,c0,dc,nn)
        
  
if sType=="gc": 
    if forcing>0:
        EN=Nxx-1
        Eu=-1
        Ev=0
        Ew=0
        EN=None

        U0=-1
        if hu<H: U2 = (H*U0-hu*U1)/(H-hu);   # h*U1+(H-h)*U2=H*U0
        else:    U2 = 0.0
    else:
        U  = None
        U0 = None
        U1 = None
        U2 = None
    if noise>0:
        if xn4>xn1:
            wnoise = domain.new_field()  # weighting for streamfunction
            if xn1==xn2 and xn3==xn4: 
                TWN = ju.cubich(x,xn1,xn4)
                ju.logger.info('wnoise: cubich {:5.2f}-{:5.2f}'.format(xn1,xn4))
            else:
                TWN = ju.cubici(x,xn1,xn2)-ju.cubici(x,xn3,xn4)
                if xn4>L:     TWN = TWN + ju.cubici(x+L,xn4,xn3)
                if xn1<0:     TWN = TWN - ju.cubici(x-L,xn2,xn1)
                ju.logger.info('wnoise: cubici {:5.2f} {:5.2f} {:5.2f} {:5.2f}'.format(xn1,xn2,xn3,xn4))
                
            wnoise['g']=TWN
            if Tx=='SinCos': wnoise.meta['x']['parity'] = 1
        else:
            wnoise=1
            ju.logger.info('wnoise: noise everywhere')
            WN=1
 
    # All the buoyancy forcing is the same controlled by x4,x5,x6,x7 and x0 and Wx
    # weight region is x0,x0+Wx to x6,x7
    # b=B    region is x4,x5 past x7 and 0 to hb+Wz/2,hb-Wz/2
    if hb>=H: hbz=np.ones_like(z)
    else:     hbz=ju.cubici(z,hb+Wz/2,hb-Wz/2)
    if lbc=='slip':
        if hu>=H: hu1=np.ones_like(z)
        else:     hu1=ju.cubici(z,hu+Wz/2,hu-Wz/2)
    else:
        if hu>=H: hu1=ju.cubici(z,-Wz/2,Wz/2)
        else:     hu1=ju.cubici(z,-Wz/2,Wz/2)-ju.cubici(z,hu-Wz/2,hu+Wz/2)
    hu2=1-hu1
    if forcing>0:
        nfb      = domain.new_field()  # Buoyancy field
        nwb      = domain.new_field()  # Buoyancy weight 
        bfx      = ju.cubici(x,x4,   x5)-ju.cubici(x,x7,x7+x5-x4) # Buoyancy field
        nwb['g'] = ju.cubici(x,x0,x0+Wx)-ju.cubici(x,x6,x7)       # Weight field
        nfb['g'] = B*bfx*hbz
        if xb is not None:
            if xb>x5: ifb = hbz*(ju.cubici(x,x4,x5)-ju.cubici(x,xb-Wx/2,xb+Wx/2))
        del(bfx)
    else:
        if xa is not None and xb is not None: ifb = hbz*(ju.cubici(x,xa-Wx/2,xa+Wx/2)-ju.cubici(x,xb-Wx/2,xb+Wx/2))
        if xa is     None and xb is not None: ifb = hbz*(1-ju.cubici(x,xb-Wx/2,xb+Wx/2))
        if xa is not None and xb is     None: ifb = hbz*(ju.cubici(x,xa-Wx/2,xa+Wx/2)-1)
    del(hbz)

    if forcing==1:
        nwu = domain.new_field()  # Weight for u
        if Ny>1: nwv = domain.new_field()  # Weight for v
        nww = domain.new_field()  # Weight for w
        
        Dx=1/Wx
        Dy=1/Wy
        Dz=1/Wz   # Need extra vertical smoothness to prevent ringing
        
        x0=Wx      # start of uniform velocity forcing
        x1=1.5+Wx # end   of uniform velocity forcing
        
        x2=3       # start of profile velocity forcing
        x3=4       # end   of profile velocity forcing start of positive buoyancy
        
        x4=4.5     # end of buoyancy forcing
        x5=x4+6    # initial condition 
        
        w1=ju.indicatorp(x,x0,   x1,   Dx,L)   # This region force uniform velocity profile
        W1=ju.indicatorp(x,x0-Wx,x1+Wx,Dx,L) 
        
        w2=ju.indicatorp(x,x2,   x3,   Dx,L)   # Force zero velocity and high and low velocity
        W2=ju.indicatorp(x,x2-Wx,x3+Wx,Dx,L)
    
        w4=ju.indicatorp(x,x0,   x3,   Dx,L)   # v velocity forcing
        W4=ju.indicatorp(x,x0-Wx,x3+Wx,Dx,L) 
        
                
        # Buoyancy 0 in regions 1 and 2 1 in region 3
        #I=np.maximum(0,np.minimum(1,(ju.indicator(x,x0,x3,Dx)+ju.indicator(x,x0,x4,Dx)*H1)))
        
    
        if Ny>1: nwv['g']  = w4
        nww['g']  = (w1 + w2)
        nwu['g']  = (w1 + w2)
        nfu['g']  = U0*(1-W2) + W2*(U1*hu1+U2*hu2)
        if Ny>1: problem.substitutions['fv'] = '0'
        problem.substitutions['fw'] = '0'
        problem.substitutions['fd'] = '0'
    elif forcing==2:
        nwu = domain.new_field()  # Weight for u
        fu        = W1
        
        if Ny>1: nwv['g']  = w1
        
        if Ny>1: problem.substitutions['fv'] = '0'
        problem.substitutions['fw'] = '0'
        problem.substitutions['fd'] = '0'
        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'

        # Impose constant velocity between 0 horizontal velocity profile over a region and the necessary
        # divergence to be compatible with zero vertical velocities
        # Vertical and transverse velocities are forced to zero
        # between x0 and x1 u=-U
        # between x2 and x3 y= U1 h<hu, U2, h>hu
        # v and w are zero everywhere in the forcing region
        # divergence is non-zero between x1 and x2
    elif forcing==3:
        nwu = domain.new_field()  # Weight for u
        nfd = domain.new_field()  # div(u)
        nfu = domain.new_field()  # u velocity 

        # x0 x3  forcing buoyancy and velocity
        # x2 x3  constant velocity and region over which velocity weight adjusts      
        # x2 x3  positive and negative divergence
        xwu = ju.cubici( x,x0,x0+Wx) - ju.cubici( x,x2,x3)        # velocity weight
        xfu = ju.cubici( x,x1,   x2) - ju.cubici( x,x3,x3+x2-x1)  # velocity value
        
        
        nfu['g'] = (xfu*(U1*hu1+U2*hu2)-(1-xfu))
        nfu.differentiate('x', out=nfd)
        nfd['g']=np.where((xx-x0)*(x3-xx)>0,nfd['g'],0)
        nwu['g'] = xwu
        if Ny>1: problem.substitutions['fv'] = '0'
        problem.substitutions['fw'] = '0'
        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'
        if xb is not None:
            if xb>x2: xfu=(ju.cubici( x,x1,   x2) - ju.cubici( x,xb-Wx,xb+Wx))
        ifu=(xfu*(U1*hu1+U2*hu2)-(1-xfu))
    
        Tdiv=ju.jint(nfd)
        del(xwu,xfu)
    elif forcing==4:
        #velocity forcing only in x0 x2 
        #buoyancy forcing only in x0 x3
        
        nwu = domain.new_field()  # Velocity forcing weight 
        nfd = domain.new_field()  # Target divergence
        
        # Velocity forcing weight is from x0,x0+wx and x1,x2 forcing weight should rise sharply but decline smoothly
        # Buoyancy forcing weight is from x0  to x5 
        # Buoyancy field is 0 from x0 to x4 and then 1 from x4 to x5
        # Velocity field is constant
        # x0 orgin of all forcing
        # x1 & x2 velocity
        # x3 & x4 divergence
        # x5 & x6 buoyancy
        n1=np.round((x1-x0)/dx); # Choose order of horizontal smoothing
        n2=np.round((x2-x1)/dx); # Choose order of horizontal smoothing
        n1 = 3
        n2 = 3
        n3 = np.int(np.round(Nz/2));       # Choose order of vertical   smoothing
    
        nwu['g'] =  (ju.cubici(x,x0,x0+Wx)-ju.cubici(x,x1,x2))     # Velocity   weight
 

        if 0:
            wd =     ju.cubici(x,x3-Wx/2,x3+Wx/2)-ju.cubici(x,x4-Wx/2,x4+Wx/2)              # Divergence field
            if hu >=H: hb1=ju.cubici(z,hu,hu)
            else:
                z1=0
                z2=H
                bt=0.9999
                while z2-z1>1e-10:
                    z3 = (z1+z2)/2
                b3,bb3,bbb3 = ju.smoothi(z3,n3,0,H)
                #print('z1 {:6.4f}, z2 {:6.4f}, z3 {:6.4f} b {:8.6f}'.format(z1,z2,z3,b3))
                if b3>bt: z2=z3
                else:     z1=z3
            hb1,hb2,hb3 = ju.smoothi(hu+z3-z,n3,0,H)
             
            Df0,DF0,Ddf0=ju.smoothi(0,n3,0,H)
            Dfh,DFh,Ddfh=ju.smoothi(hu,n3,0,H)
            DfH,DFH,DdfH=ju.smoothi(H,n3,0,H)
            DIVA=(((-U1+U2)*hu-U2*H)*DFh+U2*(H-hu)*DF0+DFH*U1*hu)/((H-hu)*DF0-H*DFh+DFH*hu)
            DIVB=-hu*(-U2+U1)*(H-hu)/((H-hu)*DF0-H*DFh+DFH*hu)
            #print('DF0 {}  DFh {}  DFH {} DIVA {} DIVB {}'.format(DF0,DFh,DFH,DIVA,DIVB))    
            
            Df,DF,Ddf=ju.smoothi(z,n3,0,H)
            
            nfd['g'] = 0*wd*(DIVA+DIVB*Df)/(x4-x3)

        #hu*U1+(H-hu)*U2=-H*U
        U2=-(H+hu*U1)/(H-hu)
        dw=ju.dcubici(x,x4,x3)  
        nfd['g'] = -dw*(1+U1*hu1+U2*hu2)
  
        if PIDX>0 and not PIDG:
            nfu = domain.new_field()  # Velocity function
            nfu['g'] = -1
        else:
            problem.substitutions['fu'] = '-U'
 
        if Ny>1: problem.substitutions['fv'] = '0'
        problem.substitutions['fw'] = '0'
        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'
        
        Tdiv=ju.jint(nfd)

        # This forcing approach forces the horizontal velocity with a constant value, -U,
        # The Buoyancy must be large enough to overcome this
    elif forcing==5:
        #velocity forcing only in x0 x2 
        
        nwu = domain.new_field()  # Velocity forcing weight 
 
        # Constant Velocity forcing  from x0,x1 to x2,x3 amd is constant
        # x0 0.0
        # x1 0.5
        # x2 1.0
        # x3 1.5
          
        xfu =ju.cubici(x,x0,x1)-ju.cubici(x,x2,x3) # Velocity forcing weight
   
        nwu['g'] = xfu

        if Ny>1: problem.substitutions['fv'] = '0'
        if Ny>1: problem.substitutions['wv'] = 'wu'
        problem.substitutions['fw'] = '0'
        problem.substitutions['ww'] = 'wu'
        problem.substitutions['fd'] = '0'
        if PIDX>0 and not PIDG:
            nfu = domain.new_field()  # Velocity function
            nfu['g']=-1
        else: problem.substitutions['fu'] = '-U'
        
        if xb is not None:
            if xb>x1: xfu=(ju.cubici( x,x0,   x1) - ju.cubici( x,xb-Wx,xb+Wx))
        ifu=xfu*(U1*hu1+U2*hu2)+(1-xfu)*U0
         
        Tdiv=None
        
        
    elif forcing==6:
        # Streamfunction forcing with velocity controlled by x0,x1,x2 and x3
        # buoyancy forcing controlled by x4, x5, x6 and x7
         
        nwu  = domain.new_field()  # Velocity weight 
        nfu  = domain.new_field()  # u velocity field
        nfw  = domain.new_field()  # v velocity field
        psi  = domain.new_field()  # stream function
        if Ty=="SinCos":            psi.meta['y']['parity'] = 1
            
        nfx=2;
        nfz=np.int(np.round(Nz*hu/(2*H)));
        nfz=2
        
        
        fx1,Fx1,dfx1    = ju.smoothi(x,nfx,x1,x3)
        fx2,Fx2,dfx2    = ju.smoothi(x,nfx,x3,x3+x3-x1)
        fx  = fx1-fx2
        Fx  = Fx1-Fx2
        dfx = dfx1-dfx2
        xwu = ju.cubici(x,x0,x0+Wx)-ju.cubici(x,x2,x3)

        nwu['g'] = xwu   # Velocity weight
    
        if lbc=='slip':
            fz1=1
            Fz1=z
            Fz1H=H
            fz2,Fz2,dfz2    = ju.smoothi(z,nfz, hu-Wz/2, hu+Wz/2)
            fz2H,Fz2H,dfz2H = ju.smoothi(H,nfz, hu-Wz/2, hu+Wz/2)
        else: 
            fz1,Fz1,dfz1    = ju.smoothi(z,2,   -Wz/2, Wz/2)
            fz1H,Fz1H,dfz1H = ju.smoothi(H,2,   -Wz/2, Wz/2)
            Fz1  = 2*Fz1-z
            Fz1H = 2*Fz1H-H
            fz2,Fz2,dfz2    = ju.smoothi(z,nfz, hu-Wz/2, hu+Wz/2)
            fz2H,Fz2H,dfz2H = ju.smoothi(H,nfz, hu-Wz/2, hu+Wz/2)

   
        FF = z+(1+U1)*(Fz1H/Fz2H*Fz2-Fz1)
        # FF = +(1+U1)*H/(H-hu)*Fz2 -z U1 % Has the correct limit as Wz=0

        psi['g'] = -z*(1-fx) - fx*FF
        psi.differentiate('z', out=nfu)  
        psi.differentiate('x', out=nfw)
        nfw['g'] = -nfw['g']

        
        if Ny>1: problem.substitutions['fv'] = '0'
        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'
        problem.substitutions['fd'] = '0'
        
         
        if xb is not None:
            if xb>x3: xfu=(ju.cubici( x,x1,   x3) - ju.cubici( x,xb-Wx,xb+Wx))
        ifu=xwu*(U1*hu1+U2*hu2)-(1-xwu)
        
        Tdiv=None
        del(psi,fx1,Fx1,dfx1,fx2,Fx2,dfx2,fx,Fx)
    elif forcing>6:
        ju.logger.error("Unknown forcing {}".format(forcing))
        ju.set_error()

ju.check_abort()

if noise>0: problem.parameters['fnoise'] = zmn
else:    problem.substitutions['fnoise'] = '0'


rfb=None
rfs=None
rfu=None
rfv=None
rfw=None
  
if rfb is None: problem.substitutions['rfb']='0'
else:              problem.parameters['rfb']=rfb
if rfs is None: problem.substitutions['rfs']='0'
else:              problem.parameters['rfs']=rfb
if rfu is None: problem.substitutions['rfu']='0'
else:              problem.parameters['rfu']=rfu
if rfv is None: problem.substitutions['rfv']='0'
else:              problem.parameters['rfv']=rfv
if rfw is None: problem.substitutions['rfw']='0'
else:              problem.parameters['rfw']=rfw

Uflg = PIDX==0 or PIDG or sType=='pm'

WN = ju.set_wfun(problem, wnoise, 'wnoise', L*H*W, Tx, Ty, AA)
WB = ju.set_wfun(problem, nwb,    'wb'    , L*H*W, Tx, Ty, AA)
WS = ju.set_wfun(problem, nws,    'ws'    , L*H*W, Tx, Ty, AA)
WU = ju.set_wfun(problem, nwu,    'wu'    , L*H*W, Tx, Ty, AA)
WV = ju.set_wfun(problem, nwv,    'wv'    , L*H*W, Tx, Ty, AA)
WW = ju.set_wfun(problem, nww,    'ww'    , L*H*W, Tx, Ty, AA)

#if p['fpidu'] and nfu is not None:
#    ju.logger.info('Setting fpidu to xwu*(-(udxx+udzz)/Re + nfu*udx+nfw*udz)')
#    udx   = nfu.differentiate('x')
#    udxx  = udx.differentiate('x')
#    udz   = nfu.differentiate('z')
#    udzz  = udz.differentiate('z')
#    rfu   = U*xwu*(-(udxx+udzz)/Re + nfu*udx+nfw*udz)
#    del(udx,udxx,udz,udzz)
rfu=0
rfw=0
rfv=0
if p['fpid']>0: problem.parameters['fpid']=p['fpid']
else:        problem.substitutions['fpid']='0'
#if p['fpidw'] and nfw is not None:
#    ju.logger.info('Setting fpidw to nwu*(-(wdxx+wdzz)/Re + nfu*wdx+nfw*wdz)')
#    wdx   = nfw.differentiate('x')
#    wdxx  = wdx.differentiate('x')
#    wdz   = nfw.differentiate('z')
#    wdzz  = wdz.differentiate('z')
#    rfw   = U*xww*(-(wdxx+wdzz)/Re + nfu*wdx+nfw*wdz)
#    del(wdx,wdxx,wdz,wdzz)
  
fb = ju.set_ffun(problem,nfb,'fb',1,Tx,Ty, 1, 1,False,AA)
fd = ju.set_ffun(problem,nfd,'fd',U,Tx,Ty,-1, 1,Uflg,AA)
fu = ju.set_ffun(problem,nfu,'fu',U,Tx,Ty,-1, 1,Uflg,AA)
fv = ju.set_ffun(problem,nfv,'fv',U,Tx,Ty, 1,-1,Uflg,AA)
fw = ju.set_ffun(problem,nfw,'fw',U,Tx,Ty, 1, 1,Uflg,AA)

if sType=='gc':
    if Ny>1: WV=WU
    if forcing==6: WW=WU
if sType=='pm':
    if Ny>1: WV=WU
    WW=WU
    
if restart:
    if rfn is not None: fnrst=ju.find_start_file(bdir,rfn)
    else:               fnrst=ju.find_start_file(bdir,name)
    ju.logger.info('Restart file :{}'.format(fnrst[0]))

if (not restart) or start_new_files:
    mode = "overwrite"
else:
    mode = "append"

if sType=='gc':
    problem.parameters['A']  = W
    fluxd='yz'
elif sType=='pm':
    problem.parameters['A']  = 1
    if   forcing==0: fluxd='yz'
    elif forcing==1: fluxd='xy'
    elif forcing==2: fluxd='yz'
    elif forcing==3: fluxd='yz'
    elif forcing==4: fluxd='yz'
    else:
        ju.logger.error('forcing should be 0, 1, 2, 3 or 4 for plumes'.format(forcing))
        ju.set_error()
else:
    ju.logger.error('sType should be gc or pm'.format(sType))
    ju.set_error()
ju.check_abort()



#https://docs.python.org/3/library/resource.html
#resource.RLIMIT_AS   The maximum area (in bytes) of address space which may be taken by the process.
#resource.RLIMIT_VMEM The largest area of mapped memory which the process may occupy.
if VMEM>0:
    memlim=np.int64(VMEM)*1024**3
    resource.setrlimit(resource.RLIMIT_AS,  (memlim,memlim))
    resource.setrlimit(resource.RLIMIT_RSS, (memlim,memlim))
    ju.logger.info('Set maximum memory: {:6.3f} GB (AS & RSS)'.format(VMEM))
 
smem=ju.MemStats()
ju.logger.info(smem)

ju.logger.info("Running ded_gc.py: Gravity current with forcing")
ju.run_info()
ju.logger.info("signals:            {}".format(signals))
ju.logger.info("parallel:           {}".format(parallel))
ju.logger.info("forcing:            {},    force: {:4.0f},".format(forcing,F))
ju.logger.info("clipB:              {},    clipu: {:4.0f},".format(clipB,clipu))
ju.logger.info("restart:            {}".format(restart))
ju.logger.info("reset:              {}".format(reset))
ju.logger.info("start_new_files:    {}".format(start_new_files))
ju.logger.info("ck_time:            {:7.4f}".format(ck_time))
ju.logger.info("wall_time:          {}".format(wall_time))
ju.logger.info("Integration scheme: {}".format(scheme))
ju.logger.info("mindt:  {:8.6f}, maxdt: {:8.6f}".format(mindt, maxdt))
ju.logger.info("dx:     {:7.5f},    dy: {:7.5f},    dz: {:7.5f}".format(dx, dy,  dz))
ju.logger.info("Wx:     {:7.5f},    Wy: {:7.5f},    Wz: {:7.5f}".format(Wx, Wy, Wz))
if WB is not None: ju.logger.info("WB:     {:7.5f}".format(WB))
if WS is not None: ju.logger.info("WS:     {:7.5f}".format(WS))
if WU is not None: ju.logger.info("WU:     {:7.5f}".format(WU))
if WV is not None: ju.logger.info("WV:     {:7.5f}".format(WV))
if WW is not None: ju.logger.info("WW:     {:7.5f}".format(WW))
if WN is not None: ju.logger.info("WN:     {:7.5f}".format(WN))
ju.logger.info("Re:     {:7.0f},    Sc: {:7.4f},    T:  {:7.1f}".format(Re, Sc,  T))
if Tdiv: ju.logger.info("Tdiv:  {:7.5e}".format(Tdiv))
ju.logger.info("B:      {:7.4f}".format(B))
ju.logger.info("noise:  {:7.4f}, noiseL: {:7.5f}, noiseT: {:7.5f}, noised {}".format(noise,noiseL,noiseT,noised))
ju.logger.info("g:      {:7.4f},     q: {:7.4f}".format(g,  q))
#ju.logger.info("gx:     {:7.4f},    gy: {:7.4f},    gz: {:7.4f}".format(gx, gy, gz))
if U1 is not None: ju.logger.info("U1: {:7.4f}".format(U1))
if U2 is not None: ju.logger.info("U2: {:7.4f}".format(U2))
#ju.logger.info("problem.variables     {}".format(problem.variables))
#ju.logger.info("problem.substitutions {}".format(problem.substitutions))
#ju.logger.info("problem.parameters    {}".format(problem.parameters))

ju.logger.info("xn1:    {:7.4f}, xn2:     {:7.4f}, xn3:    {:7.4f},  xn4: {:7.4f}".format(xn1,xn2,xn3,xn4))
ju.logger.info("x0:     {:7.4f}, x1:      {:7.4f}, x2:     {:7.4f},  x3:  {:7.4f}".format(x0,x1,x2,x3))
ju.logger.info("x4:     {:7.4f}, x5:      {:7.4f}, x6:     {:7.4f},  x7:  {:7.4f}".format(x4,x5,x6,x7))
if inoise is not None: ju.logger.info("inoise  {:7.4f}".format(inoise))
if xa     is not None: ju.logger.info("xa:     {:7.4f}".format(xa))
if xb     is not None: ju.logger.info("xb:     {:7.4f}".format(xb))

if sType=="gc":
    if hu is not None: ju.logger.info("hu:     {:7.4f},".format(hu))
    if hb is not None: ju.logger.info("hb:     {:7.4f},".format(hb))
    if U  is not None: ju.logger.info("U:      {:7.4f},".format(U))
    ju.logger.info("lbc:    {:>7s},  ubc:     {:>7s}".format(lbc,ubc))
if sType=="pm":
    ju.logger.info("Inlet shape:             {:}".format(inlet))
    ju.logger.info("Inlet radius          R: {:7.4f}".format(R))
    ju.logger.info("Inlet velocity        U: {:7.4f}".format(U))
    ju.logger.info("Inlet density         B: {:7.4f}".format(B))
    ju.logger.info("Inlet area           IA: {:7.4f}".format(IA))
    ju.logger.info("            sqrt(IA/pi): {:7.4f}".format(np.sqrt(IA/np.pi)))
    ju.logger.info("Inlet volume flux   V/A: {:7.4f}".format(QV))
    ju.logger.info("Inlet momentum flux M/A: {:7.4f}".format(QM))
    ju.logger.info("Inlet Buoyancy      B/A: {:7.4f}".format(QB))
    ju.logger.info("Inlet Buoyancy flux Q/A: {:7.4f}".format(QQ))
    if forcing==0:
        ju.logger.info("TQ1:   {:7.4f} TQ2: {:7.4f} TQ3: {:7.4e}".format(TQ1,TQ2,TQ3))
        ju.logger.info("rho=0           [{:7.4f} {:7.4f}]".format( 0,x1))
        ju.logger.info("rho=1           [{:7.4f} {:7.4f}]".format(x1,x3))
        ju.logger.info("negative div(u) [{:7.4f} {:7.4f}]".format( 0,x0))
        ju.logger.info("zero vel        [{:7.4f} {:7.4f}]".format(x0,x1))
        ju.logger.info("positive div(u) [{:7.4f} {:7.4f}]".format(x1,x2))
        ju.logger.info("constant u + N  [{:7.4f} {:7.4f}]".format(x2,x3))
        ju.logger.info("outlet velocity U2: {:7.4f}".format(U2))
 

#for s in problem.substitutions: ju.logger.info('substitutions: {:7s}: {}'.format(s,problem.substitutions[s]))
#for s in problem.parameters:    ju.logger.info('parameters:    {:7s}: {}'.format(s,problem.parameters[s]))

ju.boussinseq_equations(domain,problem,p)


#problem.meta[:]['z']['dirichlet'] = True
if Tz=='Cheb':
    if noise>0: # If we are using this as a stream function must be constant at the ends 
        problem.add_bc("right(noisef) = 0")
        problem.add_bc("left(noisef)  = 0")
    if sType=="pm":
        problem.add_bc("left(b)    = left(BF(t,x,y,h,B))*(1+left(noisef))")
        problem.add_bc("right(bz)  = 0")
        problem.add_bc("left(v)    = 0")
        problem.add_bc("right(vz)  = 0")
        problem.add_bc("left(u)    = 0")
        problem.add_bc("right(uz)  = 0")
        problem.add_bc("left(w)    = left(BF(t,x,y,h,U))*(1+left(noisef))")
        problem.add_bc("right(wz)  = 0", condition="(nx != 0) or (ny != 0)")
        problem.add_bc("integ_z(p) = 0", condition="(nx == 0) and (ny == 0)")
    else:
        problem.meta['w']['z']['dirichlet'] = True
        #problem.meta['u']['z']['dirichlet'] = True
        #problem.meta['uz']['z']['dirichlet'] = True
        #problem.meta['w']['z']['dirichlet'] = True
        #problem.meta['wz']['z']['dirichlet'] = True
        #problem.meta['b']['z']['dirichlet'] = True
        #problem.meta['bz']['z']['dirichlet'] = True
        if   lbc=='slip':
            problem.add_bc(         "left(uz)   = 0")
            if Ny>1: problem.add_bc("left(vz)   = 0")
        elif lbc=='noslip':
            if PIDX>0 and not PIDG: problem.add_bc("left(u) = left(fu)")
            elif U is not None:     problem.add_bc("left(u) = -U")
            else:                   problem.add_bc("left(u) = 0")
            if Ny>1: problem.add_bc("left(v)    = 0")
        if   ubc=='slip':
            problem.add_bc(         "right(uz)  = 0")
            if Ny>1: problem.add_bc("right(vz)  = 0")
        elif ubc=='noslip':
            if PIDX>0 and not PIDG: problem.add_bc("right(u) = right(fu)")
            elif U is not None:     problem.add_bc("right(u) = -U")
            else:                   problem.add_bc("right(u) = 0")
            if Ny>1: problem.add_bc("right(v)   = 0")

        problem.add_bc("right(w)   = 0")
        if Ny>1: 
            problem.add_bc("left(w)    = 0", condition="(nx != 0) or (ny != 0)")
            problem.add_bc("integ_z(p) = 0", condition="(nx == 0) and (ny == 0)")
        else:
            problem.add_bc("left(w)    = 0", condition="(nx != 0)")
            problem.add_bc("integ_z(p) = 0", condition="(nx == 0)")
                
        vz=list(set(problem.variables) & set(('bz','sz','fpidbz','fpidsz','fpiduz','fpiduz','fpidwz')))
        for v in vz:
            problem.add_bc("right(" + v + ")  = 0")
            problem.add_bc("left("  + v + ")  = 0")
            problem.meta[v]['z']['dirichlet'] = True
             

ju.print_problem(problem)
# Build solver
# U. M. Ascher, S. J. Ruuth, and R. J. Spiteri, Applied Numerical Mathematics (1997).
if   scheme=='RK111': solver = problem.build_solver(de.timesteppers.RK111) # 1st-order 1-stage DIRK+ERK scheme [Ascher 1997 sec 2.1]
elif scheme=='RK222': solver = problem.build_solver(de.timesteppers.RK222) # 2nd-order 2-stage DIRK+ERK scheme [Ascher 1997 sec 2.6]
elif scheme=='RK443': solver = problem.build_solver(de.timesteppers.RK443) # 3rd-order 4-stage DIRK+ERK scheme [Ascher 1997 sec 2.8]
elif scheme=='RKSMR': solver = problem.build_solver(de.timesteppers.RKSMR) # (3-ε)-order 3rd-stage DIRK+ERK scheme [Spalart 1991 Appendix]
else:
    ju.logger.error('Unknown integration scheme {}'.format(scheme))
    ju.set_error()
ju.check_abort()
ju.logger.info('Solver built for initial value problem with integration scheme {}'.format(scheme)) 

ju.logger.info('checkpointing in {}'.format(ddir))
checkpoint = Checkpoint(ddir)
checkpoint.set_checkpoint(solver, wall_dt=ck_time*60, mode=mode,parallel=ck_parallel)

# Initial conditions
u=None
v=None
w=None
b=None
s=None
uz=None
vz=None
wz=None
bz=None
sz=None
fpidu=None
fpidv=None
fpidw=None
fpidb=None
fpids=None
nf=None
nf2=None

if 'u'       in problem.variables:  u     = solver.state['u']
if 'v'       in problem.variables:  v     = solver.state['v']
if 'w'       in problem.variables:  w     = solver.state['w']
if 's'       in problem.variables:  s     = solver.state['s']
if 'b'       in problem.variables:  b     = solver.state['b']
if 'uz'      in problem.variables:  uz    = solver.state['uz']
if 'vz'      in problem.variables:  vz    = solver.state['vz']
if 'wz'      in problem.variables:  wz    = solver.state['wz']
if 'sz'      in problem.variables:  sz    = solver.state['sz']
if 'bz'      in problem.variables:  bz    = solver.state['bz']
if 'fpidu'   in problem.variables:  fpidu = solver.state['fpidu']
if 'fpidv'   in problem.variables:  fpidv = solver.state['fpidv']
if 'fpidw'   in problem.variables:  fpidw = solver.state['fpidw']
if 'fpidb'   in problem.variables:  fpidb = solver.state['fpidb']
if 'fpids'   in problem.variables:  fpids = solver.state['fpids']
if 'noisef'  in problem.variables:  nf1   = solver.state['noisef']
if 'noisef2' in problem.variables:  nf2   = solver.state['noisef2']

if restart:
    dt = ju.restart(fnrst, solver)
else:
    ju.logger.info('Initialising fields')


    # Get pointers to the simulation variables and set them to zero
    if b is not None: b['g'] = 0
    if s is not None: s['g'] = 0
    if u is not None: u['g'] = 0
    if v is not None: v['g'] = 0
    if w is not None: w['g'] = 0
    
    if ifb is None and 'fb' in problem.parameters and 'wb' in problem.parameters:
        ifb = problem.parameters['fb']['g']*problem.parameters['wb']['g']
        b.set_scales(AA)
        b['g']   = ifb
    else:
        b['g']   = ifb
        
    if ifs is not None: s['g']   = ifs
    if ifu is not None: u['g']   = U*ifu
    if ifv is not None: v['g']   = U*ifv
    if ifw is not None: w['g']   = U*ifw
    if nf2 is not None: nf2['g'] = anoise(1)
            
    if bz is not None: b.differentiate('z', out=bz)
    if sz is not None: s.differentiate('z', out=sz)
    if uz is not None: u.differentiate('z', out=uz)
    if vz is not None: v.differentiate('z', out=vz)
    if wz is not None: w.differentiate('z', out=wz)

    if fpidu is not None: fpidu = rfu 
    if fpidv is not None: fpidv = rfv 
    if fpidw is not None: fpidw = rfw
    if fpidb is not None: fpidb = rfb 
    if fpids is not None: fpids = rfs
 
del(ifb,ifs,ifu,ifv,ifw)
            
if inoise is not None:
    ju.logger.info('Adding noise to velocities with amplitude {:8.3f}'.format(inoise))
    u.set_scales(AA)
    if Ny>1: v.set_scales(AA)
    w.set_scales(AA)
    zn=domain.new_field()
    if   p['Tx'] == "SinCos": zn.meta['x']['parity'] = -1
    else:                     zn.meta['x']['parity'] = 1
    if Ny>1 and p['Ty'] == "SinCos": zn.meta['y']['parity'] = 1
    zn.set_scales(AA)
    zn['g']=inoise*anoise(AA)
    zn['g']=zn['g']-ju.jint(zn)/problem.parameters['Volume']
    u['g']=u['g']+zn['g']
    
    zn['g']=inoise*anoise(AA)
    zn['g']=zn['g']-ju.jint(zn)/problem.parameters['Volume']
    w['g']=w['g']+zn['g']

    if Ny>1:
        if  p['Ty'] == "SinCos":  zn.meta['y']['parity'] = -1
        zn['g']=inoise*anoise(AA)
        zn['g']=zn['g']-ju.jint(zn)/problem.parameters['Volume']
        v['g']=v['g']+zn['g']


    del(zn)
    
maxu=ju.jmax(np.abs(u['g']))
if Ny>1: maxv=ju.jmax(np.abs(v['g']))
else:    maxv=0
maxw=ju.jmax(np.abs(w['g']))
ju.logger.info('Maximum velocities: [{:6.3f},{:6.3f},{:5.2f}]'.format(maxu,maxv,maxw))

if maxu==0 and maxv==0 and maxw==0: dt =np.sqrt(mindt*maxdt)
else:                               dt = 0.5/max(maxu/dx,maxv/dy,maxw/dz)
dt=max(mindt,min(maxdt,dt))
ju.logger.info('Initial timestep:   [{:6.4f},{:6.4f},{:6.4f}]'.format(mindt,dt,maxdt))

if clipB: b['g']=np.maximum(0,np.minimum(B,b['g']))
if clipu>0:
    u['g']=np.maximum(-clipu,np.minimum(clipu,u['g']))
    if Ny>1: v['g']=np.maximum(-clipu,np.minimum(clipu,v['g']))
    w['g']=np.maximum(-clipu,np.minimum(clipu,w['g']))

zmn.args = [dt]

if reset:
    ju.logger.info('Reset')
    solver.sim_time=0
    solver.iteration=0
    ju.logger.info("setting solver.sim_time=0")
    ju.logger.info("setting solver.iteration=0")
    for field in solver.state.fields:
        if field.name[0]=='a':
            ju.logger.info('zeroing: {}'.format(field.name))
            xxx=solver.state[field.name]
            xxx['g']=0

 
timea.append(np.float(time.time()))
wt=np.maximum(ju.get_sec(wall_time)-200-timea[0]+timea[1],60)
# Integration parameters
solver.stop_sim_time  = T
CFL=flow_tools.CFL(solver, initial_dt=dt, cadence=5, safety=0.8, max_change=1.5, min_change=0.5, max_dt=maxdt, min_dt=mindt)
if Ny>1: CFL.add_velocities(('u','v','w'))
else:    CFL.add_velocities(('u','w'))

solver.stop_wall_time = wt
solver.stop_iteration = np.inf

analysis_tasks=da.set_analysis_tasks(ddir,solver,problem,max_writes,fh_parallel,p)

# Flow properties Print some information about the simulation every cadence iterations
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
#if Ny>1: flow.add_property("sqrt(u*u + v*v + w*w)*RL/nu", name='Re')
#else:    flow.add_property("sqrt(u*u + w*w)*RL/nu", name='Re')
if WB is not None: flow.add_property("(b-fb)**2*abs(wb)", name='FB')
if WS is not None: flow.add_property("(s-fs)**2*abs(ws)", name='FS')
if WU is not None: flow.add_property("(u-fu)**2*abs(wu)", name='FU')
if WV is not None: flow.add_property("(v-fv)**2*abs(wv)", name='FV')
if WW is not None: flow.add_property("(w-fw)**2*abs(ww)", name='FW')


if noise>0: flow.add_property("wnoise*(noisex*noisex + noisey*noisey + noisez*noisez)",   name='NS')

timea.append(np.float(time.time()))
ju.logger.info('Initialisation time: {:f}'.format(timea[2]-timea[0]))

ju.check_abort()

smem=ju.MemStats()
ju.logger.info(smem)

sim_iter_start=solver.iteration
sim_time_start=solver.sim_time
ju.write_time(sim_time_start)
ju.logger.info('Starting simulation loop at sim time {:8.3f}, wall time: {}'.format(sim_time_start,time.time()))
st0=sim_time_start
if PIDT>0: nt0=np.floor(st0/PIDT)
else:      nt0=0
    
gshape = domain.dist.grid_layout.global_shape(scales=AA)
slices = domain.dist.grid_layout.slices(scales=AA)
#print('ju.mpirank {}, slices {}, gshape {}'.format(ju.mpirank,slices,gshape))


tp= np.float(time.time())

if PIDST>0:
    if np.int(np.floor(st0/PIDST) % 2)==0: PIDX=PIDS1
    else: PIDX=PIDS2
    ju.logger.info("Stepping PID input to calibrate controller")
    ju.logger.info("PIDST:   {:7.3f},   PIDS1: {:7.4f},   PIDS2: {:7.4f}".format(PIDST, PIDS1, PIDS2))


# We need to make sure that it is fine if the front position is initially a NaN    
b=solver.state['b']
tX0=None
X=None

if sType=='gc':
    X=ju.calc_front(b)
    if not np.isfinite(X): X=0
    tX0=solver.sim_time
    ju.logger.info("Front position X: {:7.4f}".format(X))

X0=X
tX1=tX0


if PIDX>0 and ju.mpirank==0:
    if PIDIT==0:
        if PIDG: PIDIT=g
        else:    PIDIT=U
    if PIDG: ju.logger.info("PID controller for front position is active on g: {:7.4f}".format(g))
    else:
        ju.logger.info("PID controller for front position is active on U: {:7.4f}".format(U))
    ju.logger.info("PIDX: {:7.4f},  PIDIT: {:7.4f}, PIDDD: {:7.4f}".format(PIDX, PIDIT, PIDDD))
    ju.logger.info("PIDT: {:7.4f},   PIDP: {:7.4f}   PIDI: {:7.5f},  PIDD: {:7.4f}".format(PIDT, PIDP, PIDI, PIDD))
    if PIDG: pid = jpid.pid(PIDP, PIDI, PIDD, PIDX, PIDIT, PIDDD, X,  g, st0, 0.05,20,5,-1)
    else:    pid = jpid.pid(PIDP, PIDI, PIDD, PIDX, PIDIT, PIDDD, X,  U, st0, 0.05, 1,0.5,1)

NS=None
FB=None
FS=None
FU=None
FV=None
FW=None
UF=None
cRe=None

dforce=True

if PIDG or PIDX==0: del(fd,fs,fb,fu,fv,fw)
del(wnoise,nwb,nws,nwu,nwv,nww)
del(nfb,nfd,nfu,nfv,nfw)
del(xx,yy,zz)
del(x,y,z)

ju.logger.info('Starting main loop');
t1=timea[2]-30

if False:
    maxu=ju.jmax(np.abs(u['g']))
    if Ny>1: maxv=ju.jmax(np.abs(v['g']))
    else:    maxv=0
    maxw=ju.jmax(np.abs(w['g']))
    ju.logger.info('Maximum velocities: [{:6.3f},{:6.3f},{:6.3f}]'.format(maxu,maxv,maxw))

check_finite=True
fflag=False
st0=solver.sim_time
while solver.ok and not ju.abort:
    if check_finite:
        if b is not None: fflag =fflag or not np.all(np.isfinite(b['g']))
        if s is not None: fflag =fflag or not np.all(np.isfinite(w['g']))
        if u is not None: fflag =fflag or not np.all(np.isfinite(u['g']))
        if v is not None: fflag =fflag or not np.all(np.isfinite(v['g']))
        if w is not None: fflag =fflag or not np.all(np.isfinite(w['g']))
        if fflag:
            ju.set_error('Not finite')
            ju.logger.info('Not finite')
        ju.check_abort()

    #mu=ju.jint(u)/(L*H*W)
    #ju.logger.info('mean u velocity {:}'.format(mu))
    if tX1 is not None: tX0=copy.copy(tX1)
    if sType=='gc' and np.isfinite(X): XO=copy.copy(X)
    
    dt = CFL.compute_dt()
    zmn.args = [dt]
    solver.step(dt)
    ju.jdot()
    
    if clipB: b['g']=np.maximum(0,np.minimum(B,b['g']))
    if clipu>0:
        u['g']=np.maximum(-clipu,np.minimum(clipu,u['g']))
        if Ny>1: v['g']=np.maximum(-clipu,np.minimum(clipu,v['g']))
        w['g']=np.maximum(-clipu,np.minimum(clipu,w['g']))
    st=solver.sim_time
    if PIDT>0: nt=np.floor(st/PIDT)
    else:      nt=nt0+1
    if sType=='gc':
        XX=ju.calc_front(b)
        if np.isfinite(XX): X=copy.copy(XX)
        if np.isfinite(X) and np.isfinite(X0) and tX1>tX0:
            tX1=st
            UC=(X-XO)/(tX1-tX0)
            UF = UC + np.exp(-dt)*(UF-UC) # Calculate smoothed front velocity
            cRe=RL*np.abs(U+UF)*Re
    else: cRe=RL*U*Re
    if EN is not None:
        EU=ju.calc_uerror(domain,EN,u,v,w,Eu,Ev,Ew,Edy,Edz)
        ju.logger.info('EU = {}'.format(EU))
        
    if WB is not None: FB=np.sqrt(flow.volume_average('FB')/WB)
    if WS is not None: FS=np.sqrt(flow.volume_average('FS')/WS)
    if WU is not None: FU=np.sqrt(flow.volume_average('FU')/WU)
    if WV is not None: FV=np.sqrt(flow.volume_average('FV')/WV)
    if WW is not None: FW=np.sqrt(flow.volume_average('FW')/WW)
    if WN is not None: NS=np.sqrt(flow.volume_average('NS')/WN)
 
    if cRe is not None:
        if not np.isfinite(cRe):
            ju.logger.error('Reynolds number is infinite or NaN')
            ju.set_error()
    
    ju.check_abort_file()
    ju.test_abort()
    
    t2= np.float(time.time())
 
    if ju.mpirank==0:
        if PIDST>0: 
            if np.int(np.floor(st/PIDST) % 2)==0: NPIDX=PIDS1
            else: NPIDX=PIDS2
            if NPIDX!=PIDX:
                ju.logger.info('Updating PIDX to {:6.4f} from {:6.4f}'.format(NPIDX,PIDX));
                PIDX=NPIDX
                pid.update_XT(NPIDX)
                
    if PIDX>0:
        PIDP = np.float64(0)
        if ju.mpirank==0: PIDP,PIDIT,PIDDD=pid.update(X,st)
        PIDP       = ju.comm.bcast(PIDP,0)
        if np.isfinite(PIDP):
            PIDIT      = ju.comm.bcast(PIDIT,0)
            PIDDD      = ju.comm.bcast(PIDDD,0)
            p['PIDIT'] = PIDIT
            p['PIDDD'] = PIDDD
            if PIDG:
                g                       = copy.copy(PIDP)
                p['g']                  = g
                nfg['g']                = g
                problem.parameters['g'] = nfg
            else:
                U                        = copy.copy(PIDP)
                p['U']                   = U
                problem.parameters['U']  = U
                if fu is not None:
                    nfu['g'] = U*fu
                    problem.parameters['fu'] = nfu
                if fv is not None:
                    nfv['g'] = U*fv
                    problem.parameters['fv'] = nfv
                if fw is not None:
                    nfw['g'] = U*fw
                    problem.parameters['fw'] = nfw
                if fd is not None:
                    nfd['g'] = U*fd
                    problem.parameters['fd'] = nfd
    if t2>t1+30:
        ju.write_time(st)
        ju.jup()
        si = 'n:{:6d}, t:{:7.3f}, e: {:7.1e}'.format(solver.iteration,st,3600*(st-st0)/(t2-t1))
        t1=copy.copy(t2)
        st0=copy.copy(st)                                         
        smem=ju.MemStats()
        si += ju.prange('b',b)
        si += ju.prange('s',s)
        if cRe is not None: si += ', Re:{:4.0f}'.format(cRe)
        if FB  is not None: si += ', FB:{:6.4f}'.format(FB)
        if FS  is not None: si += ', FS:{:6.4f}'.format(FS)
        if FU  is not None: si += ', FU:{:6.4f}'.format(FU)
        if FV  is not None: si += ', FV:{:6.4f}'.format(FV)
        if FW  is not None: si += ', FW:{:6.4f}'.format(FW)
        if NS  is not None: si += ', NS:{:6.3f}'.format(NS)
        if X   is not None: si += ', X:{:6.3f}'.format(X)
        if PIDX>0:
            if PIDG: si += ', g: {:7.4f}'.format(g)
            else:    si += ', U: {:6.4f}'.format(U)
            ju.write_param(p)  
        if dtstats>0: ju.flush_stats()
        ju.logger.info(si)


    if dtstats>0: ju.write_stats(solver.iteration,st,t2,X,U,g,dt,cRe,NS,FB,FU,FV,FW)

if dtstats>0: ju.close_stats()
        
smem=ju.MemStats()
ju.logger.info(smem)
        
        
if PIDX>0: ju.write_param(p)

ju.test_error()
if not ju.jerror:
    ju.logger.info("Final step to save state")
    final_checkpoint = Checkpoint(ddir, checkpoint_name='final')
    final_checkpoint.set_checkpoint(solver, wall_dt=0.1, mode="overwrite", parallel=ck_parallel)
    solver.step(dt) #clean this up in the future...works for now.

ju.end_stats(timea,sim_time_start,solver,domain)  # Print statistics about the run     

if not ck_parallel:
    post.merge_process_files(str(ddir)+'/final/', cleanup=True)
    post.merge_process_files(str(ddir)+'/checkpoint/', cleanup=True)
    
if not fh_parallel:
    for task in analysis_tasks:
        post.merge_process_files(task.base_path, cleanup=True)

ju.finish(timea[0],solver)
            
ju.logger.info('Complete')

smem=ju.MemStats()
ju.logger.info(smem)

