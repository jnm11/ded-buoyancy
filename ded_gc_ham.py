#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DNS 2d and 3d Navier stokes for plumes and gravity currents
"""
DNS 2d and 3d Navier stokes for plumes and gravity currdedents
Usage: 
   ded_gc.py [options] NAME

Options:
  -h, --help            : Show this help message
      --debug=<>        : 0 non debuggging information the higher the number the more info
      --debug_flow=<>   : Debug main program flow
      --save_initial=<> : Save initial state 
      --sType=<>        : Type of simulation
      --preset=<>       : Run a simulation with preset parameters
      --setPID          : Set PID defaults
      --scheme=<>       : Integration scheme
  -R, --Re=<>           : Reynolds number  
  -F, --force=<>        : Strength of forcing
      --forceb=<>       : Multiplier to strength of forcing for bouyancy
      --forces=<>       : Multiplier to strength of forcing for sedimentation
  -C, --conservative=<> : Conservative formulation of advection terms
      --gclf=<>         : Add body force to match divergence around x=0
      --gcrf=<>         : Add body force to match divergence around x=L
      --fuspecial=<>    : Special test case of psi = sqrt(L^2+H^2) U sin(2 Pi x/L) sin(2 Pi z/H)
      --fuspecial1=<>   : parameter 1 for fuspecial
      --fuspecial2=<>   : parameter 2 for fuspecial
      --clipu=<>        : Clip velocities to less than this absolute value
      --clipB=<>        : Clip buoyancy between 0 and B or B and 0
      --clipS=<>        : Clip sediment between 0 and B or B and 0
      --rescale=<>      : rescale so that velocities have rms of 1
      --minB=<>         : Minimum B
      --maxB=<>         : Maximum B
      --minS=<>         : Minimum S
      --maxS=<>         : Maximum S
      --m1=<>           : MPI grid 1
      --m2=<>           : MPI grid 2
      --maxdt=<>        : Maximum time step
      --mindt=<>        : Minimum time step
      --Scb=<>          : Schmidt number for bouyancy
      --Scs=<>          : Schmidt number for sediment
      --Peb=<>          : Peclet number for bouyancy sets Scb
      --Pes=<>          : Peclet number for sediment sets Scs
  -H, --height=<>       : Box height       
  -W, --width=<>        : Box width        
  -L, --length=<>       : Box length       
      --radius=<>       : Inlet radius      
      --fixI=<>         : Fix fluid injection profile
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
      --oddg=<>         : g is a field with odd parity
  -T, --Time=<>         : Simulation time  
  -f, --forcing=<>      : Forcing type     
      --f7s=<>          : Adjustment factor for velocity smoothing in f7 simulation
      --xa=<>           : Start of region to fill with material
      --xb=<>           : End of region to fill with material
      --dimple=<>       : Amplitude of dimples to add
      --dimplewy=<>     : y wavenumber decay
      --dimplewz=<>     : z wavenumber decay
      --ddiv=<>         : number dynamic divergence modes
      --db=<>           : dynamic buoyancy
      --resetdb=<>      : reset dynamic buoyancy
      --r0=<>           : Inside radius for adding fluid
      --x0=<>           : Start of velocity uniform region
      --x1=<>           : End of velocity uniform region
      --x2=<>           : Start of density region
      --x3=<>           : end of density region
      --x4=<>           : 
      --x5=<>           : 
      --x6=<>           : 
      --x7=<>           : 
      --wu=<>           : force u velocity
      --wv=<>           : force v velocity
      --wwl=<>          : force w velocity left
      --wwr=<>          : force w velocity right
      --fbmult=<>       : force bouyancy field using multiplication
      --fbmin=<>        : force bouyancy field using min
      --fbmax=<>        : force bouyancy field using max
      --fsmult=<>       : force sediment field using multiplication
      --fsmin=<>        : furce sediment field using min
      --fsmax=<>        : force sediment field using max
      --xn1=<>          : start noise forcing region
      --xn2=<>          : 
      --xn3=<>          : 
      --xn4=<>          : end noise forcing region
      --Wx=<>           : x transition width
      --Wy=<>           : y transition width
      --Wz=<>           : z transition width
      --Wr=<>           : radial transition width
  -c, --dtcheck=<>      : Checkpoint time minutes       
  -q, --angle=<>        : Slope angle (degrees)    
  -g, --gravity=<>      : Gravity
      --lbc=<>          : Lower velocity boundary condition 
      --rbc=<>          : Upper velocity boundary condition 
      --blbc=<>         : Bouyancy left boundary condition 
      --brbc=<>         : Bouyancy right boundary condition 
      --slbc=<>         : Sediment left boundary condition 
      --srbc=<>         : Sediment right boundary condition 
  -t, --time=<>         : Max runtime 0 means infinite 
      --signals=<>      : Catch signals      
      --parallel=<>     : Use parallel HDF   
  -s, --start_new_files : Start new files while checkpointing 
  -r, --restart         : Restart  
      --rfn=<>          : Restart file name
  -p, --param           : Load from parameter file
      --pfn=<>          : Parameter file name              
      --reset           : reset solve sim time and integrated variables
  -B, --buoyancy=<>     : Buoyancy
  -S, --sediment=<>     : Sediment
      --SV=<>           : Sediment velocity
      --Ski=<>          : Stokes number for particle inertial effects viscous term
      --Skp=<>          : Stokes number for particle inertial effects pressure term
      --Skg=<>          : Stokes number for particle density effects
      --bV=<>           : Buoyancy velocity
      --dnoise=<>       : amplitude for direct noise forcing
  -n, --noise=<>        : amplitude    for noise forcing
      --noiseL=<>       : length scale for noise forcing
      --noised=<>       : xy, xz, or yz plane to impose noise forcing from stream function
      --noiseT=<>       : time scale for noise forcing
      --inoise=<>       : initial velocity noise
      --bnoise=<>       : initial bouyancy noise
      --snoise=<>       : initial sediment noise
      --sinitial=<>     : initial sediment value
      --bntype=<>       : type of initial bouyancy noise, gaussian or uniform
  -i, --inlet=<>        : inlet type gaussian, circle, square, 
      --AA=<>           : anti-aliasing factor for computation
      --AAS=<>          : anti-aliasing factor for saving state
      --AAJ=<>          : anti-aliasing factor for saving data
      --VMEM=<>         : set maximum VMEM in GB (=0 unlimited)  
      --avrg=<>         : average variables
      --dtb=<>          : time frequency for writing out 3d bouyancy data
      --dts=<>          : time frequency for writing out 3d sediment data
      --dtd=<>          : time frequency for writing out 3d divergence data
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
      --dtjadv=<>       : Flag for recording advective terms
      --dtjsw=<>        : Flag for recording strain rate and vorticity
      --dtjd1=<>        : Flag for recording first order derivatives
      --dtjd2=<>        : Flag for recording second order derivatives
      --dtj=<>          : time frequency for writing out fields
      --dtjdivyz=<>     : time frequency for writing out div field at top of plume simulation
      --dtjcheck=<>     : time frequency for checkpointing in direct format
      --dtjx=<>         : time frequency for writing out direct integrated x
      --dtjy=<>         : time frequency for writing out direct integrated y
      --dtjz=<>         : time frequency for writing out direct integrated z
      --dtjr=<>         : time frequency for writing out direct integrated r
      --dtjxy=<>        : time frequency for writing out direct integrated xy
      --dtjxz=<>        : time frequency for writing out direct integrated xz
      --dtjyz=<>        : time frequency for writing out direct integrated yz
      --dtjxyz=<>       : time frequency for writing out direct integrated xyz
      --dtja=<>         : time frequency for writing out averaged fields
      --dtjar=<>        : time frequency for writing out averaged integrated r
      --dtjax=<>        : time frequency for writing out averaged integrated x
      --dtjay=<>        : time frequency for writing out averaged integrated y
      --dtjaz=<>        : time frequency for writing out averaged integrated z
      --dtjaxy=<>       : time frequency for writing out averaged integrated zy
      --dtjayz=<>       : time frequency for writing out averaged integrated yz
      --dtjaxz=<>       : time frequency for writing out averaged integrated xz
      --dtjaxyz=<>      : time frequency for writing out averaged integrated xyz
      --dtjar=<>        : time frequency for writing out averaged integrated r
      --dtjb=<>         : time frequency for writing out direct density
      --dtjd=<>         : time frequency for writing out direct divergence
      --dtjs=<>         : time frequency for writing out direct sediment
      --dtjp=<>         : time frequency for writing out direct pressure
      --dtju=<>         : time frequency for writing out direct u
      --dtjv=<>         : time frequency for writing out direct v
      --dtjw=<>         : time frequency for writing out direct w
      --dtjWx=<>        : time frequency for writing out direct x vorticity
      --dtjWy=<>        : time frequency for writing out direct y vorticity
      --dtjWz=<>        : time frequency for writing out direct z vorticity
      --dtjS=<>         : time frequency for writing out direct strain
      --dtjE=<>         : time frequency for writing out direct enstrophy
      --dtjaS=<>        : time frequency for writing out averaged strain
      --dtjaE=<>        : time frequency for writing out averaged enstrophy
      --dtjsev=<>       : time frequency for writing out state_ev data
      --dtjavar=<>      : time frequency for writing out auxiliary variables
      --dtavrg=<>       : time frequency for writing out time averaged data
      --dtleft=<>       : time frequency for writing out left boundary data
      --dtright=<>      : time frequency for writing out right boundary data
      --dtpm=<>         : time frequency for writing out plume data
      --dtgc=<>         : time frequency for writing out gravity current data
      --dtslice=<>      : time frequency for outputting slice data
      --dtforce=<>      : time frequency for outputting forcing data
      --Fu=<>           : Body force in x direction
      --Fv=<>           : Body force in y direction
      --Fw=<>           : Body force in z direction
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
      --gmin=<>         : minimum g
      --gmax=<>         : maximum g
      --Umin=<>         : minimum U
      --Umax=<>         : maximum U
      --dT=<>           : Timescale for adding divergence 
      --dL=<>           : Lengthscale for adding divergence 
      --pmss=<>         : Standard deviation for moving plume source
      --pmsl=<>         : Decay constant for decorrelation of plume source
      --pmzr=<>         : Radial zero buoyancy for bouyancy in plume simulations
      --pmzt=<>         : Top zero buoyancy for bouyancy in plume simulations
      --pmzrt=<>        : Radial-top zero buoyancy for bouyancy in plume simulations
      --pmIu=<>         : Include average u velocity in divergence update
      --topd=<>         : Dynamic evolution of yz divergence top of plume
      --extx=<>         : Method for extending in x direction
      --exty=<>         : Method for extending in y direction
      --extz=<>         : Method for extending in z direction
"""

#ded_gc-1.1.py -R 1000 -W 0.2 -H 1 -N 64 --time=00:01:00 test
#ded_gc-1.1.py -r -p test
#ded_gc.py -H1 --Nx=64 -W0.1 -R100 -ck=1 001
#mpiexec -n 4 ded_gc-1.31.py --qgc gc1
#mpiexec -n 4 ded_gc-1.31.py --qplume1 pm1
#mpiexec -n 4 ded_gc-1.1.py --time=00:05:00 -N 128 -R 1000 -W 0.2 -a test
#mpiexec -n 4 ded_gc-1.1.py -r -p --time=00:01:00 test
#mpiexec -n 16 ded_gc-1.1.py -r -p --time=00:05:00 041
#mpiexec -n 4 python3 -m mpi4py $HOME/python/Dedalus/ded_gc-1.1.py --time=01:00:00 -N 128 -R 1000 -W 0.2 -a test
#mpiexec -n 4 python3 -m mpi4py $HOME/python/Dedalus/ded_gc-1.1.py -r -p --time=00:05:00 test
import time
tinitial=time.time()
import numpy as np
timea = [np.float(time.time())]
ver=np.float(1.32)

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
import jrec as ja
import jnoise as jn
from scipy.interpolate import griddata
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools  import post
from dedalus.core   import operators
from pathlib import Path
import os
import gc


if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

import logging

import warnings

wti0=np.float(time.time())      # Initial wall time

ju.logger=logging.getLogger(__name__)

debug_gc=np.bool(os.getenv('DED_GC',False))
if debug_gc:
    ju.logger.info('Debugging garbage collection');
    gc.set_debug(gc.DEBUG_LEAK)

debug_mem=np.bool(os.getenv('DED_MEM',False))
if debug_mem:
    ju.logger.info('Debugging memory')
    import psutil

    
debug_pympler=np.bool(os.getenv('DED_PYMPLER',False))
if debug_pympler:
    ju.logger.info('Debugging memory')
    from pympler import summary
    from pympler import muppy
    from pympler import tracker


debug_tracemalloc=np.bool(os.getenv('DED_TRACEMALLOC',False))
if debug_tracemalloc:
    import tracemalloc
    tracemalloc.start()
    

np.seterr(divide='ignore', invalid='ignore') # Ignore /0 and 0/0
#np.seterr(divide='print', invalid='print') # Ignore /0 and 0/0
#np.seterr(divide='raise',invalid='raise')

if args['--start_new_files']:                      start_new_files = True
else:                                              start_new_files = False
if args['--restart'] or args['--rfn'] is not None: restart = True
else:                                              restart = False
if args['--reset']:                                reset = True
else:                                              reset = False


param,name = dp.name_type(args,restart)
sType=param['sType']

bdir,ddir=ju.set_base_dir(sType,name)

ju.check_version(ver)

debug_flow   = args.pop('--debug_flow','False')
save_initial = args.pop('--save_initial','False')
resetdb      = args.pop('--resetdb','False')
pfn=args.pop('--pfn',  None)
rfn=args.pop('--rfn',  None)


if args['--param'] or pfn:
    pfn=ju.find_param_file(bdir,name,pfn)
    paramf=ju.read_param(pfn)
    param=dp.merge_param(param,paramf)

if args['--preset']: param=dp.preset(param,args['--preset'])
if args['--setPID']: param=dp.set_PID(param)



param=dp.read_args(param,args)
param=dp.validate(param)

# Parameters for loading initial state don't save
rparam=dict()
rparam['extx'] = param.pop('extx','zero') 
rparam['exty'] = param.pop('exty','zero') 
rparam['extz'] = param.pop('extz','zero') 
rparam['Tx']   = dp.get_param(param,'Tx')
rparam['Ty']   = dp.get_param(param,'Ty')
rparam['Tz']   = dp.get_param(param,'Tz')

debug    = param.pop('debug',  None)
if debug is not None: ju.logger.info('Debug level is {}'.format(debug))

ju.Barrier()
ju.logger.info('Parameters validated')



ju.write_status("Running")
ju.write_files()
ju.node_info()

iparam=copy.copy(param)

scheme   = dp.get_param(param,'scheme')
forcing  = dp.get_param(param,'forcing')
F        = dp.get_param(param,'force')
maxdt    = dp.get_param(param,'maxdt')
mindt    = dp.get_param(param,'mindt')
clipu    = dp.get_param(param,'clipu')
clipB    = dp.get_bool_param(param,'clipB')
clipS    = dp.get_bool_param(param,'clipS')
B        = dp.get_param(param,'B')
S        = dp.get_param(param,'S')
dtjdivyz = dp.get_param(param,'dtjdivyz')
f7s      = dp.get_param(param,'f7s')

fixI     = dp.get_bool_param(param,'fixI')

U1 = dp.get_param(param,'U1')
hu = dp.get_param(param,'hu')
hb = dp.get_param(param,'hb')

Tx       = dp.get_param(param,'Tx')
Ty       = dp.get_param(param,'Ty')
Tz       = dp.get_param(param,'Tz')
Nx       = dp.get_param(param,'Nx')
Ny       = dp.get_param(param,'Ny')
Nz       = dp.get_param(param,'Nz')

dim=np.int(Nx>1) + np.int(Ny>1) +np.int(Nz>1)  

Nxx       = np.int(np.round(Nx*param['AA']))
Nyy       = np.int(np.round(Ny*param['AA']))
Nzz       = np.int(np.round(Nz*param['AA']))
 
H        = param['H']
W        = param['W']
L        = param['L']
alpha    = dp.get_param(param,'alpha')
Wx  = dp.get_param(param,'Wx')
Wy  = dp.get_param(param,'Wy')
Wz  = dp.get_param(param,'Wz')
Wr  = dp.get_param(param,'Wr')
r0  = dp.get_param(param,'r0')
AA  = dp.get_param(param,'AA')
AAJ = dp.get_param(param,'AAJ')
x0  = dp.get_param(param,'x0')
x1  = dp.get_param(param,'x1')
x2  = dp.get_param(param,'x2')
x3  = dp.get_param(param,'x3')
x4  = dp.get_param(param,'x4')
x5  = dp.get_param(param,'x5')
x6  = dp.get_param(param,'x6')
x7  = dp.get_param(param,'x7')
fuspecial  = dp.get_param(param,'fuspecial')
fuspecial1 = dp.get_param(param,'fuspecial1')
fuspecial2 = dp.get_param(param,'fuspecial2')

noiseL  = dp.get_param(param,'noiseL')
noised  = dp.get_param(param,'noised')
noiseT  = dp.get_param(param,'noiseT')
noise   = dp.get_param(param,'noise' )
dnoise  = dp.get_param(param,'dnoise')
xn1     = dp.get_param(param,'xn1')
xn2     = dp.get_param(param,'xn2')
xn3     = dp.get_param(param,'xn3')
xn4     = dp.get_param(param,'xn4')

q       = dp.get_param(param,'q')
g       = dp.get_param(param,'g')
Re      = dp.get_param(param,'Re')
U       = dp.get_param(param,'U')
T       = dp.get_param(param,'T')
BV      = dp.get_param(param,'BV')    
SV      = dp.get_param(param,'SV')    
fbmult  = dp.get_bool_param(param,'fbmult')    
fsmult  = dp.get_bool_param(param,'fsmult')    
fbmin   = dp.get_bool_param(param,'fbmin')    
fsmin   = dp.get_bool_param(param,'fsmin')    
fbmax   = dp.get_bool_param(param,'fbmax')    
fsmax   = dp.get_bool_param(param,'fsmax')    
avrg    = dp.get_bool_param(param,'avrg')
pmzt    = dp.get_bool_param(param,'pmzt')
pmzr    = dp.get_bool_param(param,'pmzr')
pmzrt   = dp.get_bool_param(param,'pmzrt')
pmIu    = dp.get_bool_param(param,'pmIu')
pwu     = dp.get_bool_param(param,'wu')
pwv     = dp.get_bool_param(param,'wv')
pwwl    = dp.get_bool_param(param,'wwl')
pwwr    = dp.get_bool_param(param,'wwr')


if F is not None:
    if maxdt*F>2:
        maxdt=2/F
        ju.logger.info('Limiting maxdt due to forcing to {:8.3e}'.format(maxdt))

# Don't save parameters that setup initial state
inoise   = param.pop('inoise', None)
bnoise   = param.pop('bnoise', None)
snoise   = param.pop('snoise', None)
sinitial = param.pop('sinitial', None)
bntype   = param.pop('bntype', None)
xa       = param.pop('xa',     None)
xb       = param.pop('xb',     None)
dimple   = param.pop('dimple', None)
dimplewy = param.pop('dimplewy', None)
dimplewz = param.pop('dimplewz', None)


VMEM      = dp.get_param(param,'VMEM')
sType     = dp.get_param(param,'sType')
signals   = dp.get_param(param,'signals')
parallel  = dp.get_param(param,'parallel')
dtcheck   = dp.get_param(param,'dtcheck')
dtjcheck  = dp.get_param(param,'dtjcheck')
dtstats   = dp.get_param(param,'dtstats')

PIDX  = dp.get_param(param,'PIDX')
PIDT  = dp.get_param(param,'PIDT')
PIDG  = dp.get_bool_param(param,'PIDG')
PIDST = dp.get_param(param,'PIDST')
Wr    = dp.get_param(param,'Wr')


PIDD  = dp.get_param(param,'PIDD')
PIDI  = dp.get_param(param,'PIDI')
PIDP  = dp.get_param(param,'PIDP')
PIDS1 = dp.get_param(param,'PIDS1')
PIDS2 = dp.get_param(param,'PIDS2')

db    = dp.get_bool_param(param,'db')
ddiv  = dp.get_bool_param(param,'ddiv')
oddg  = dp.get_bool_param(param,'oddg')
topd  = dp.get_bool_param(param,'topd')

fpid  = dp.get_param(param,'fpid')
fpidu = dp.get_bool_param(param,'fpidu')
fpidv = dp.get_bool_param(param,'fpidv')
fpidw = dp.get_bool_param(param,'fpidw')
fpidb = dp.get_bool_param(param,'fpidb')
fpids = dp.get_bool_param(param,'fpids')
conservative=dp.get_bool_param(param,'conservative')
rescale=dp.get_bool_param(param,'rescale')
gclf=dp.get_bool_param(param,'gclf')
gcrf=dp.get_bool_param(param,'gcrf')

state_ev=dict()  # Extra state variables
avar=dict()      # Auxiliary variables to save if requested with dtjvar

if  f7s is None: f7s=1
        
if PIDX is not None:
    state_ev['PIDIT']=0.0
    state_ev['PIDDD']=0.0
    if PIDG: state_ev['g']=dp.get_param(param,'g')
    else:    state_ev['U']=dp.get_param(param,'U')

    if PIDG:
        if 'Umax' in param: param.pop('Umax')
        if 'Umin' in param: param.pop('Umin')
    else:
        if 'gmax' in param: param.pop('gmax')
        if 'gmin' in param: param.pop('gmin')
else:
    if 'PIDIT' in param: param.pop('PIDIT')
    if 'PIDDD' in param: param.pop('PIDDD')
    if 'Umax' in param: param.pop('Umax')
    if 'Umin' in param: param.pop('Umin')
    if 'gmax' in param: param.pop('gmax')
    if 'gmin' in param: param.pop('gmin')

tops=None
rmstopd=None
bmult=None
smult=None
bmultr=None
smultr=None
bmultl=None
smultl=None
wmult=None
bmin=None
bmax=None
smin=None
smax=None
pmif=None

dt = None
Idivzs=None
Igdw=None
u1=None
u2=None
fuzz=None
    
ju.make_lock()
if signals: ju.init_signals()



param=dp.remove_None(param)
ju.logger.info("Writing parameters")
ju.write_param(param)  

dx=L/Nx;
dy=W/Ny;
dz=H/Nz;
dV=dx*dy*dz
dr=np.sqrt(dy**2+dz**2)

# Should add a test to see if h5py supports parallel
max_writes=5
#https://stackoverflow.com/questions/47072859/how-to-append-data-to-one-specific-dataset-in-a-hdf5-file-with-h5py


     
    # 3D Boussinesq hydrodynamics

#print("Tz: {}, Nx: {}, Ny: {}, Nz: {}".format(Tz,Nx,Ny,Nz))
#print("L: {}, W: {}, H: {}".format(L,W,H))

if sType=='bg': domain,problem = ju.burgers_init(param)
else:           domain,problem = ju.boussinseq_init(param)

#gyz=ju.gridyz(param,param['AA'],domain.dist.grid_layout.ext_mesh[1:])

ju.logger.info('Calculating Integration weights')
Iw, Iwg  = ju.int_weights( 1,domain)
Iww,Iwwg = ju.int_weights(AA,domain)

xx = domain.grids(AA)[ 0]
zz = domain.grids(AA)[-1]
x  = domain.grids( 1)[ 0]  # These will contain only the local grid coordinates
z  = domain.grids( 1)[-1]
if Ny>1:
    dim=3
    yy = domain.grids(AA)[ 1]
    y  = domain.grids( 1)[ 1]
    gIy  = Iw[1]
    gIyy = Iww[1]
    yslices  = domain.dist.grid_layout.slices(scales=1)[1]
    yyslices = domain.dist.grid_layout.slices(scales=AA)[1]
    Nyl      = domain.dist.grid_layout.local_shape(scales=1)[1]
    Nyyl     = domain.dist.grid_layout.local_shape(scales=AA)[1]
    gyyslices =ju.comm.gather(yyslices,root=0)
else:
    dim=2
    yy=0
    y=0
    gIy=1
    gIyy=1
    Nyl=1
    Nyyl=1

gIx  = Iw[0]
gIxx = Iww[0]
gIz  = Iw[-1]
gIzz = Iww[-1]

xslices  = domain.dist.grid_layout.slices(scales=1)[0]
xxslices = domain.dist.grid_layout.slices(scales=AA)[0]
zslices  = domain.dist.grid_layout.slices(scales=1)[-1]
zzslices = domain.dist.grid_layout.slices(scales=AA)[-1]
gzzslices = ju.comm.gather(zzslices,root=0)
Nxl      = domain.dist.grid_layout.local_shape(scales=1)[0]
Nzl      = domain.dist.grid_layout.local_shape(scales=1)[-1]
Nxxl     = domain.dist.grid_layout.local_shape(scales=AA)[0]
Nzzl     = domain.dist.grid_layout.local_shape(scales=AA)[-1]

zcolor   =  zslices.start==0 and  xslices.start==0 and  Nyl>0
zzcolor  = zzslices.start==0 and xxslices.start==0 and Nyyl>0
zcomm    = ju.comm.Split(np.int( zcolor), ju.mpirank)
zzcomm   = ju.comm.Split(np.int(zzcolor), ju.mpirank)
 
if Ny>1:
    ycolor   =  yslices.start==0 and  xslices.start==0 and Nzl>0
    yycolor  = yyslices.start==0 and xxslices.start==0 and Nzzl>0
    ycomm    = ju.comm.Split(np.int( ycolor), ju.mpirank)
    yycomm   = ju.comm.Split(np.int(yycolor), ju.mpirank)

else:
    ycolor=True
    yycolor=True
    ycomm=ju.comm
    yycomm=ju.comm

if Ny>1:
    szg   = (Nx,Ny,Nz)
    szzg  = (Nxx,Nyy,Nzz)
    szl   = (Nxl,Nyl,Nzl)
    szzl  = (Nxxl,Nyyl,Nzzl)
    szx   = (Nx,1,1)
    szy   = (1,Ny,1)
    szz   = (1,1,Nz)
    szxx  = (Nxx,1,1)
    szyy  = (1,Nyy,1)
    szzz  = (1,1,Nzz)
    szxl  = (Nxl,1,1)
    szyl  = (1,Nyl,1)
    szzl  = (1,1,Nzl)
    szxxl = (Nxxl,1,1)
    szyyl = (1,Nyyl,1)
    szzzl = (1,1,Nzzl)
 
else:
    szg   = (Nx,Nz)
    szzg  = (Nxx,Nzz)
    szl   = (Nxl,Nzl)
    szzl  = (Nxxl,Nzzl)
    szx  = (Nx,1)
    szz  = (1,Nz)
    szxx = (Nxx,1)
    szzz = (1,Nzz)
    szxl  = (Nxl,1)
    szzl  = (1,Nzl)
    szxxl = (Nxxl,1)
    szzzl = (1,Nzzl)

gz  = np.zeros( szz,dtype=np.float64)
gzz = np.zeros(szzz,dtype=np.float64)
    
if ycolor:
    zcounts  = np.array(ycomm.gather(Nzl,  root=0))
    ycomm.Gatherv(sendbuf=z, recvbuf=(gz, zcounts), root=0)
else: zcounts=np.array([])
if yycolor:
    zzcounts = np.array(yycomm.gather(Nzzl, root=0))
    yycomm.Gatherv(sendbuf=zz, recvbuf=(gzz, zzcounts), root=0)
    zzcounts = yycomm.bcast(zzcounts, root=0)
else: zzcounts=np.array([])

topdfield=True

if topd:

    #https://materials.jeremybejarano.com/MPIwithPython/collectiveCom.html#scatterv-and-gatherv
    if topdfield:
        state_ev['topd'] = gyz.domain.new_field()
        avar['topu']     = gyz.domain.new_field()
        avar['topv']     = gyz.domain.new_field()
        avar['topw']     = gyz.domain.new_field()
        avar['topd']     = gyz.domain.new_field()
        #ju.logger.info('topd.shape {}'.format(avar['topd']['g'].shape))
        #quit(1)
    else:
        state_ev['topd']=np.ones((Nyy,Nzz))
        topu=np.zeros((Nyy,Nzz),dtype=np.float64)
        topv=np.zeros((Nyy,Nzz),dtype=np.float64)
        topw=np.zeros((Nyy,Nzz),dtype=np.float64)
        topd=np.zeros((Nyy,Nzz),dtype=np.float64)
        avar['topu']=topu
        avar['topv']=topv
        avar['topw']=topw
        avar['topd']=topd
    #yyzzstart  = np.array(ju.comm.gather(list(domain.dist.grid_layout.start(scales=AA)[1:]),root=0))
    #yyzzcounts = np.array(ju.comm.gather(list(domain.dist.grid_layout.local_shape(scales=AA)[1:]),root=0))
    #    if ju.mpirank==0:
    #ju.logger.info('yyzzstart[:,0]: {}'.format(yyzzstart[:,0]))
    #ju.logger.info('yyzzstart[:,1]: {}'.format(yyzzstart[:,1]))
    #ju.logger.info('yyzzcounts[:,0]: {}'.format(yyzzcounts[:,0]))
    #ju.logger.info('yyzzcounts[:,1]: {}'.format(yyzzcounts[:,1]))
    #yyzzstart  = tuple(yyzzstart[:,0]+yyzzstart[:,1]*Nyy)
    #yyzzcounts = tuple(yyzzcounts[:,0]*yyzzcounts[:,1])
    #yyzzstart  = tuple(np.cumsum(yyzzcounts)-yyzzcounts[0])
    #yyzzstart  = tuple(tuple(yyzzstart[:,0]),tuple(yyzzstart[:,1]))
    #yyzzcounts = tuple(tuple(yyzzcounts[:,0]),tuple(yyzzcounts[:,1]))
    #ju.logger.info('yyzzstart:      {}'.format(yyzzstart))
    #ju.logger.info('yyzzcounts:     {}'.format(yyzzcounts))
        
    
    #ju.logger.info('rr.shape {} Nyyl {} Nzzl {} Nyyl*Nzzl {}'.format(rr.shape,Nyyl,Nzzl,Nyyl*Nzzl))
    #ju.comm.Gatherv(sendbuf=rr.reshape(Nyyl,Nzzl),  recvbuf=(grr,yyzzcounts,yyzzstart,ju.MPI.DOUBLE),root=0)
    #ju.logger.info('rr.shape {} Nyyl {} Nzzl {} Nyyl*Nzzl {}'.format(rr.shape,Nyyl,Nzzl,Nyyl*Nzzl))
    #yyzzstart  = ju.comm.bcast(yyzzstart,0)
    #yyzzcounts = ju.comm.bcast(yyzzcounts,0)
    
    #ju.logger.info('yyzzcounts: {} {}'.format(yyzzcounts.shape,yyzzcounts))
       
if not np.isscalar(Iw[-1]):
    if ycolor:
        gIz  = np.zeros((Nz, ),dtype=np.float64)
        ycomm.Gatherv(sendbuf=Iw[-1],  recvbuf=(gIz,   zcounts), root=0)
    if yycolor:
        gIzz = np.zeros((Nzz,),dtype=np.float64)
        yycomm.Gatherv(sendbuf=Iww[-1], recvbuf=(gIzz, zzcounts), root=0)
gIz  = ju.comm.bcast(gIz,0)
gIzz = ju.comm.bcast(gIzz,0)
gz   = ju.comm.bcast(gz,0)
gzz  = ju.comm.bcast(gzz,0)

if ju.mpirank==0 and False:
    ju.logger.info('gIz {}'.format(gIz))
    ju.logger.info('gIzz {}'.format(gIzz))
    ju.logger.info('gz {}'.format(gz))
    ju.logger.info('gzz {}'.format(gzz))
    
if Ny>1:
    gy  = np.zeros( szy,dtype=np.float64)
    gyy = np.zeros(szyy,dtype=np.float64)
    if zcolor:
        ycounts = np.array( zcomm.gather(Nyl,  root=0))
        zcomm.Gatherv(sendbuf=y, recvbuf=(gy, ycounts), root=0)
    if zzcolor:
        yycounts = np.array(zzcomm.gather(Nyyl, root=0))
        zzcomm.Gatherv(sendbuf=yy, recvbuf=(gyy, yycounts), root=0)
    if not np.isscalar(Iw[1]):
        if zcolor:
            gIy = np.zeros((Ny,),dtype=np.float64)
            zcomm.Gatherv(sendbuf= Iw[1], recvbuf=( gIy,  ycounts), root=0)
        if zzcolor:
            gIyy = np.zeros((Nyy,),dtype=np.float64)
            zzcomm.Gatherv(sendbuf=Iww[1], recvbuf=(gIyy, yycounts), root=0)
    gIy =ju.comm.bcast(gIy,0)
    gIyy=ju.comm.bcast(gIyy,0)
    gy =ju.comm.bcast(gy,0)
    gyy=ju.comm.bcast(gyy,0)



nfx=10
nfy=10
nfz=10
nfr=10

maxkx=0.0
maxky=0.0
maxkz=0.0
if Tx != 'Cheb' and Nx>1: maxkx = max(domain.bases[0].wavenumbers)   # maximum  x wavenumber
if Ty != 'Cheb' and Ny>1: maxky = max(domain.bases[1].wavenumbers)   # maximum  y wavenumber
if Tz != 'Cheb' and Nz>1: maxkz = max(domain.bases[-1].wavenumbers)  # maximum  z wavenumber
maxkr=np.sqrt(maxky**2+maxkz**2)
xmaxkx  = ju.ffta_x(maxkx,nfx)  # Find the width of a region using nfx and maxkx
xmaxkx2 = ju.ffta_x(maxkx/2,nfx)  # Find the width of a region using nfx and maxkx

# Does the divergence change dynamically
dfd = ddiv or (sType=='pm' and ('pmss' in param or  forcing==7))
if dfd: ju.logger.info('Dynamic divergence')

JA=ja.jrec(domain,ddir,param,reset,dfd)
ju.logger.info('ja.jrec(domain,ddir,param,reset,Iw) complete')

   

 
dxnoise = None
dynoise = None
dznoise = None

if dnoise is not None:
    ju.logger.info('noised: {}'.format(noised))
    ju.logger.info('noiseL: {}'.format(noiseL))
    ju.logger.info('noiseT: {}'.format(noiseT))
    ju.logger.info('dnoise: {}'.format(noiseT))
    if 'x' in noised:
        dxnoise = jn.jnoise('xyz',[noiseL,noiseL,noiseL],1/noiseT,dnoise,domain)
        problem.parameters['noiseu'] = dxnoise.f
        state_ev['dxnoise']=dxnoise.f
    if 'y' in noised:
        dynoise = jn.jnoise('xyz',[noiseL,noiseL,noiseL],1/noiseT,dnoise,domain)
        problem.parameters['noisev'] = dynoise.f
        state_ev['dynoise']=dynoise.f
    if 'z' in noised:
        dznoise = jn.jnoise('xyz',[noiseL,noiseL,noiseL],1/noiseT,dnoise,domain)
        problem.parameters['noisew'] = dznoise.f
        state_ev['dznoise']=dznoise.f
    
    #djnoise.print()
    #slices    = domain.dist.coeff_layout.slices(scales=1)
    #kx        = domain.bases[0].wavenumbers[slices[0],np.newaxis,np.newaxis]
    #ky        = domain.bases[1].wavenumbers[np.newaxis,slices[1],np.newaxis]
    #kz        = domain.bases[2].wavenumbers[np.newaxis,np.newaxis,slices[2]]
    #dnoisewx  = np.exp(-(noiseL*kx)**2)
    #dnoisewy  = np.exp(-(noiseL*ky)**2)
    #dnoisewz  = np.exp(-(noiseL*kz)**2)
    #Idnoisewx = np.sqrt(ju.jsum(dnoisewx**2))/2#/domain.dist.coeff_layout.global_shape(scales=1)[0]
    #Idnoisewy = np.sqrt(ju.jsum(dnoisewy**2))/2#/domain.dist.coeff_layout.global_shape(scales=1)[1]
    #Idnoisewz = np.sqrt(ju.jsum(dnoisewz**2))/2#/domain.dist.coeff_layout.global_shape(scales=1)[2]
    #dnoisewx  = dnoisewx/Idnoisewx
    #dnoisewy  = dnoisewy/Idnoisewy
    #dnoisewz  = dnoisewz/Idnoisewz
    #dnoisef1   = domain.new_field()
    #ju.logger.info('dnoisewx.shape {}'.format(dnoisewx.shape))
    #ju.logger.info('dnoisewy.shape {}'.format(dnoisewy.shape))
    #ju.logger.info('dnoisewz.shape {}'.format(dnoisewz.shape))

elif (noise is not None) and ( (sType=='pm' and not forcing==1) or (sType=='gc') ):
    if noised=='xy':
        problem.substitutions['noiseu'] = ' wnoise*dy(noisef2)'
        problem.substitutions['noisev'] = '-wnoise*dx(noisef2)'
    elif noised=='xz':
        problem.substitutions['noiseu'] = ' wnoise*dz(noisef2)'
        problem.substitutions['noisew'] = '-wnoise*dx(noisef2)'
    else:         
        problem.substitutions['noisev'] = ' wnoise*dz(noisef2)'
        problem.substitutions['noisew'] = '-wnoise*dy(noisef2)'
            
if noise is not None:
    problem.parameters['noise']  = noise
    problem.parameters['noiseT'] = noiseT
    problem.parameters['noiseL'] = noiseL
    
problem.parameters['kappau']    = 1/Re
problem.parameters['fnu']   = 0/Re/1000
problem.parameters['L']     = L
problem.parameters['W']     = W
problem.parameters['H']     = H
if U is not None: problem.parameters['U']     = U
problem.parameters['Volume'] = W*H*L

#problem.meta[:][‘z’][‘dirichlet’] = True
#if PIDX==0 or PIDG: problem.meta['U']['x','y']['constant']


# Get pointers to the grid coordinates
#print("domain.grid(0) {} ".format(dir(domain.grid(0))))
#print("domain.dist.grid_layout.local_shape {} ".format(dir(domain.dist.grid_layout)))
#print("domain.dist.grid_layout.local_shape {} ".format(domain.dist.grid_layout.local_shape(scales=1)))
#domain.dist.grid_layout.local_shape(scales=AA)
#print("domain.dist.grid_layout.local_shape {}".format(domain.dist.grid_layout.local_shape(scales=AA)))
#print('xx.shape {}'.format(xx.shape));
#domain.dist.grid_layout(scales=1)
#print("domain.dist.grid_layout".format(domain.dist.grid_layout))




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

if Ty=='SinCos': Py=1
else:            Py=0
if Tz=='SinCos': Pz=1
else:            Pz=0



if Ny>1: dA=dz2*L/Nx*W/Ny
else:    dA=dz2*L/Nx

Edy=dy/W
Edz=dz2/H

def ssnoise(f,wx,wy,wz,sc):
    lshape = domain.dist.coeff_layout.local_shape(scales=sc)
    f['c'] = wx*wy*wz*rand.standard_normal(lshape)
    f.set_scales(sc)
    return(f['g'])

def anoise(sc):
    lshape = domain.dist.grid_layout.local_shape(scales=sc)
    x = rand.standard_normal(lshape)
    return(x)

def unoise(low,high,sc):
    lshape = domain.dist.grid_layout.local_shape(scales=sc)
    x = rand.uniform(low,high,lshape)
    #if x.size>0: print('x.min {:} x.max {:}'.format(x.min(),x.max()))
    return(x)


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
        #slices = domain.dist.grid_layout.slices(scales=AA)
        x = np.divide(anoise(AA),np.sqrt(dA))#[slices[1]]))
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
r    = None
wr   = None
rr   = None
wrr  = None
wux1 = None
wux2 = None
wux3 = None
wux4 = None
wux5 = None
wux6 = None
wux7 = None
Bwx  = None
Bwxx = None
Iux1 = None
Iux2 = None
Iux3 = None
Iux4 = None
Iux5 = None
nwdr  = None
inwdr = None

fb = None
fs = None
fu = None
fv = None
fw = None
fd = None


U2 = None
wnoise = None
nwb = None
nwd = None
nws = None
nwu = None
nwv = None
nww = None

nfg = None
nfb = None
nfs = None
nfu = None
nfv = None
nfw = None
nfd = None
rfu = None
rfv = None
rfw = None
nsvx=None

ifb = None
ifs = None
ifu = None
ifv = None
ifw = None

nyz = None
EN=None

c1=None
c2=None
u1=None
u2=None

BM=1
SM=1

dxx=L/Nxx
dyy=W/Nyy
dzz=H/Nzz
dAA=dyy*dzz
dVV=dxx*dAA
 
mun=None
vf=None

if forcing is not None and (B is not None or S is not None):
    if not fbmult:
        if not fbmax:              nfb = domain.new_field() # Target b
        if not fbmin or not fbmax:
            nwb = domain.new_field() # Target b
            ju.set_parity(nwb,Tx,Ty,Tz,1,1,1)
                
    
if sType=='pm':
    if 'pmss' in param: # moving source for the plume
        state_ev['SY']=param['pmss']*np.random.normal()
        state_ev['SZ']=param['pmss']*np.random.normal()
    R     = param['radius']
    RL=2*R # length scale for calculating Reynolds number
    inlet = param['inlet']
    IA = math.pi*R**2 # Target area
    R0=R
    InletA=1
    if ju.mpirank==0:  (R,InletA)=ju.solve_inlet(gyy,gzz,inlet,R,maxkr,nfr)
    R=ju.comm.bcast(R,0)
    InletA=ju.comm.bcast(InletA,0)
    ju.logger.info('Radius changed from {:6.4f} to {:6.4f}, amplitude {:6.4f}'.format(R0,R,InletA))
    ( wr,wrmax, r)=ju.pinlet( y, z,inlet,R,maxkr,nfr)
    (wrr,wrmax,rr)=ju.pinlet(yy,zz,inlet,R,maxkr,nfr)
    ju.check_error()
    wr =InletA*wr
    wrr=InletA*wrr
    BM=InletA*wrmax
    SM=InletA*wrmax
    #Adjust area
    if Ty=='SinCos':
        IA=IA/2
        HW=W
        FW=2*W
    else:
        HW=W/2
        FW=W

    if Tz=='SinCos':
        IA=IA/2
        HH=H
        FH=2*H
    else:
        HH=H/2
        FH=H

    IQ   = ju.jsum(wr     )*dy*dz/IA
    IQQ  = ju.jsum(wr*wr  )*dy*dz/IA
    IIQ  = ju.jsum(wrr    )*dyy*dzz/IA
    IIQQ = ju.jsum(wrr*wrr)*dyy*dzz/IA
    IAU  = IA*U
    IAUU = IA*U**2
    if noise is not None:
        wnoise = domain.new_field()  # weighting for noise
        wnoise['g']=ju.stepwf(x,nfx,xn1,xn2,xn3,xn4,L)*wr/wrmax

    if forcing is not None and (B is not None or S is not None):
        if forcing !=7:
            Bwxx=ju.stepwf( xx, nfx, x5, x6,    x7,    x7+(x6-x5),L)
            Bwx =ju.stepwf(  x, nfx, x5, x6,    x7,    x7+(x6-x5),L)
            nfb['g'] = wr*Bwx
            nwb['g'] = ju.stepwf( x, nfx, x4, x4+Wx, x7-Wx, x7,L)
            
        if xb is not None:
            if B is not None:
                domain.new_field()
                ifb['g'] = wr
            if S is not None: ifs = ifb

        if S is not None and B is not None:
            problem.substitutions['ws']   = 'wb'
            problem.substitutions['fs']   = 'fb'
            Swxx = Bwxx
            Swx  = Bwx


    if forcing==7:
        # u velocity field is odd
        # g          is constant
        # b and S are odd
        # Density is given inlet profile between x6 and x7
        # Density is set to zero between between x4 and x5
        # Fluid is inserted between x2 and x3 and horizontal velocities forced to zero 
        # Fluid is removed  between x0 and x1
        nwv  = domain.new_field() # Weight for v and w
        nfd  = domain.new_field() # divergence field
        M=7
        if Tx != 'SinCos':
            ju.logger.info('Tx must be SinCos for forcing 7')
            ju.set_error('Tx must be SinCos for forcing 7')
            ju.test_error()
        #wux1 = ju.ffta_heaviside(  xx,maxkx/2,nfx)[1]             # Bottom fluid injection region
        #wux2 = ju.ffta_heaviside(L-xx,maxkx/2,nfx)[1]             # Top fluid removal
        #wux3 = ju.ffta_heaviside(L-xx,maxkx/4,nfx)[1]             # Top fluid test region                                                 
        #math.pi/2=x*k/M/2
        #M*math.pi  =x*maxkx/2
        wux1 = ju.ffta_df(  xx, maxkx/2, M)[0]  # Bottom fluid source region
        wux2 = ju.ffta_df(L-xx, maxkx/2, M)[0]  # Top    fluid sink   region
        wux3 = ju.ffta_df(L-xx, maxkx/4, M)[0]  # Top    fluid test   region                                                 
        wux4 = ju.ffta_heaviside(xx-xmaxkx,maxkx/2,nfx)[0]-ju.ffta_heaviside(xx-L+xmaxkx,maxkx/2,nfx)[0] # Side fluid injection region
        wux6 = ju.ffta_heaviside(xx,       maxkx/2,nfx)[0]-ju.ffta_heaviside(xx-L,       maxkx/2,nfx)[0] # Side fluid injection region
        wux5 = ju.ffta_df(  xx, maxkx/4, M)[0] + ju.ffta_df(  L-xx, maxkx/2, M)[0]  # Top and bottom region in which to force horizontal velocities to zero                          
        Iwux1=wux1.sum()*dxx
        Iwux2=wux2.sum()*dxx
        Iwux3=wux3.sum()*dxx
        wux1=wux1/Iwux1
        wux2=wux2/Iwux2
        wux3=wux3/Iwux3
        wux5=wux5/wux5.max()

        if fixI: # Fix injection profile
            X0    = 2*M*math.pi/maxkx # Top of fluid source region
            #pmif=np.power(np.maximum(xx-X0,0),2/3)*ju.ffta_heaviside(L-X0-xx,maxkx,nfx)[0]
            pmif=np.power(np.maximum(xx,X0),-1/5)*ju.ffta_heaviside(L-X0-xx,maxkx,nfx)[0]*ju.ffta_heaviside(xx-X0,maxkx,nfx)[0]
            pmif=pmif/np.sqrt((pmif**2).sum()*dxx)
        
        nwv.set_scales(AA)
        nwv['g'] = wux5
        problem.substitutions['ww']   = 'wv'
        

        nfd.set_scales(AA)
        nfd['g'] = copy.copy(wrr*wux1-wux2*IIQ*IA/(H*W))
        #ju.logger.info('Infd {:}'.format(dVV*ju.jsum(nfd['g']/IA)))
        #ju.logger.info('wux1 {:}'.format(dxx*wux1.sum()))
        #ju.logger.info('wux2 {:}'.format(dxx*wux2.sum()))
        #ju.logger.info('wux3 {:}'.format(dxx*wux3.sum()))

        
        # Create function for injecting fluid that integrates to 1 across an area
        RR=min(HW,HH)-param['r0']
        theta=np.arctan2(yy,zz)
        #ju.logger.info('RR: {} Wr: {} HW: {} HH: {} q: {}'.format(RR,Wr,HW,HH,q))
        state_ev['divx']=np.zeros(szxx,dtype=np.float64)
        onewdr  = 1-ju.smoothwf(rr,nfx,RR-xmaxkx/2,RR+xmaxkx)
        nwdr    = np.nan_to_num(ju.smoothwf(rr,nfx,RR-xmaxkx/2,RR+xmaxkx)/(np.minimum(HW/np.abs(np.sin(theta)),HH/np.abs(np.cos(theta)))-RR));
        Inwdr = ju.jsum(nwdr)*dAA # RR*pi+2*HH*asinh(HW/HH)+2*HW*asinh(HH/HW);
        nwdr  = nwdr/Inwdr
        del(theta)
        Anwdr = 1/(ju.jsum(nwdr**2)*dAA)
        #ju.logger.info('Anwdr: {:6.3f}'.format(Anwdr/IA))


        if nfb is None: nfb = domain.new_field()
        if oddg: nfb.meta['x']['parity'] =  1
        else:    nfb.meta['x']['parity'] = -1
        nfb.set_scales(AA)

        if oddg: Bwxx = ju.ffta_heaviside(L/2-xx,maxkx/2,nfx)[0]
        else:    Bwxx = 2*ju.ffta_heaviside(xx,maxkx/2,nfx)[0]-ju.ffta_heaviside(xx-L/2,maxkx/2,nfx)[0]-1
        nfb['g'] = wrr*Bwxx # Buoyancy field

        #nwb['g'] = ju.ffta_heaviside(xmaxkx-xx,maxkx/2,nfx)[0]
        #ju.ffta_heaviside(xmaxkx-xx,maxkx/2,nfx)[0]
        wux7=ju.ffta_df(xx,maxkx/4  ,7)[0]           # Bottom buoyancy forcing region
        wux7=wux7/wux7.max()

        if fbmult:
            if B is not None: problem.substitutions['wb'] = 'wv'
            if S is not None: problem.substitutions['ws'] = 'wv'
            M=7
            wux8=ju.ffta_df(L-xx,maxkx  ,M)[0]   # Top buoyancy forcing region
            wux8=wux8/wux8.max()
            MM=np.int64(np.ceil(M*math.pi/maxkx/dxx))
            ju.logger.info('MM: {}, xx[Nxx-MM]: {}'.format(MM, xx[Nxx-MM]))
            if B is not None:
                bmultr  = np.minimum(1,np.maximum(0,wux8[slice(Nxx-MM, Nxx)]))
                bmultrs = slice(Nxx-MM, Nxx)

            if S is not None:
                smultr  = np.minimum(1,np.maximum(0,wux8[slice(Nxx-MM, Nxx)]))
                smultrs = slice(Nxx-MM, Nxx)
        else:
            nwb.set_scales(AA)
            nwb['g'] = wux7
            #wux8=ju.ffta_heaviside(L-xx,maxkx/2,nfx)[1]  # Top    buoyancy forcing region
            wux8=ju.ffta_df(L-xx,maxkx  ,M)[0]   # Top buoyancy forcing region
            wux8=wux8/wux8.max()
            if pmzt:  nwb['g'] += wux8
            if pmzrt: nwb['g'] += wux8*ju.smoothwf(rr,nfx,RR-xmaxkx/2,RR+xmaxkx)
            if pmzr:  nwb['g'] += ju.smoothwf(rr,nfx,RR-xmaxkx/2,RR+xmaxkx)
            nwb['g'] = ju.clip01(nwb['g'])
       
    if forcing==6: # Dynamically adjust the divergence so as to make the velocity field correct

        nfu  = domain.new_field() # Target u
        nwu  = domain.new_field() # Weight for u
        nfd  = domain.new_field() # divergence field
        rfu  = domain.new_field() # Forcing field
        #dudx = domain.new_field() 
 
        #      0
        # [x0,x1] velocity transition region [x2,x3] detransition to ensure periodicity
        # [0 x2 x3 x4] forcing region
        
        # Velocity is forced between x0 and x3
        # Gaussian profile is forced between x1-x2 and x3
        # Divergence can vary between x0 and x1
        # Divergence is added to match forced Gaussian field
        
        nwu.set_scales(AA)
        nfu.set_scales(AA)
        nfd.set_scales(AA)
        # x0 x1 region where velocity is sampled
        # x1 x2 region where constant velocity profile is imposed
        # x2 x3 region where gaussian velocity profile is imposed
        wux1     = ju.stepwf(xx,nfx,   Wx,   2*Wx,x2+2*Wx,x2+3*Wx, L)   # region where divergence is non zero for dx(fu)
        wux2     = ju.stepwf(xx,nfx,   x0, x0+Wx, x3-Wx,    x3, L)          # Region where velocity is imposed  x1-x2 constant, x2-x3 Gaussian
        wux3     = ju.stepwf(xx,nfx,   x2, x2+Wx, x2+Wx,x2+2*Wx,L)          # region to test where velocity should be Gaussian
        wux4     = ju.stepwf(xx,nfx,   x1,    x2,    x3+4*Wx,x3+x2-x1+4*Wx, L)    # Switch between Gaussian velocity and constant velocity
        wux5     = ju.stepwf(xx,nfx,x3-Wx,    x3,  L-Wx,     L, L)          # region for injecting fluid at edges
        wnx      = ju.stepwf(x,nfx,xn1,xn2,xn3,xn4,L)
        wnxx     = ju.stepwf(xx,nfx,xn1,xn2,xn3,xn4,L)
        nwu['g'] = wux2

        rwux5=np.reshape(wux5,(Nxx,))
        
        Iux1     =  wux1/(dxx*H*W*np.sum(wux1,axis=0))       # Make integral dV=1
        Iux2     =  wux2/(dxx*H*W*np.sum(wux2,axis=0))       # Make integral dV=1
        Iux3     =  wux3/(dxx*H*W*np.sum(wux3,axis=0))       # Make integral dV=1
        Iux4     =  wux4/(dxx*H*W*np.sum(wux4,axis=0))       # Make integral dV=1
        Iux5     =  wux5/(dxx*H*W*np.sum(wux5,axis=0))       # Make integral dV=1
 
        state_ev['qvx']   = np.zeros((Nxx,),dtype=np.float)
        state_ev['Fx']    = 0
        problem.parameters['Fx']=state_ev['Fx']

        problem.parameters['rfu'] = rfu
        # Half height and half width
        
        # Choose a weight function so that the r integral is independent of theta
        # and integrates to 1
        RR=min(HW,HH)-param['r0']
        q=np.arctan2(yy,zz)
        nwdr = np.nan_to_num(ju.smoothwf(rr,nfx,RR-Wr/2,RR+Wr/2)/(min(HW/abs(sin(q)),HH/abs(cos(q)))-RR));
        Inwdr = ju.jsum(nwdr)*dAA # RR*pi+2*HH*asinh(HW/HH)+2*HW*asinh(HH/HW);
        nwdr = nwdr/Inwdr
        del(q)
        
        #nwdr  = np.nan_to_num(ju.smoothwf(rr,nfx,R1,R2)/(rr/np.maximum(np.abs(zz)/H,np.abs(yy)/W)-R1-R2))
        #nwdr = 1-(1 - ju.smoothwf(np.abs(yy),nfx,W/2-2*Wr,W/2-Wr))*(1-ju.smoothwf(np.abs(zz),nfx,H/2-2*Wz,H/2-Wz))
        #nwdr=ju.remove_q_dep(yy,zz,nwdr,6)
        #nwdr = nwdr/np.sqrt(1+np.abs(2*yy*zz/np.maximum(1e-10,rr**2)))

        IAA=ju.jsum(wrr)*dAA
        inwdr=jtheta.periodic_gaussianu(1/(param['radius']+0.1*(L-x3)),yy,zz,W,H)
        #inwdr=1-ju.smoothwf(rr,nfx,R1,R2)
        inwdr = H*W*inwdr/(dAA*ju.jsum(inwdr))

        
        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'

        # Not necessary but ensures that the force data files are correct
        
        mu=ju.jsum(wrr)/ju.jsum(inwdr)
        nfu['g'] = wrr*wux4+mu*(1-wux4)*inwdr      # Impose constant region and then Gaussian region
        nfu.differentiate('x',out=nfd)
        nfd['g'] = nfd['g'] * wux1

        inwdr2=1-ju.smoothwf(r,nfx,R1,R2)
        inwdr2 = H*W*inwdr2/(dA*ju.jsum(inwdr2))
        #ifu = mu*inwdr2+(wr-mu*inwdr2)*ju.stepwf(x,nfx,x1,    x2,    L-Wx,L, L)
        #ifv=0
        #ifw=0
        
    elif forcing==5:
        # x0 = 0.5 [0.0 0.5] negative div(u)
        # x1 = 1.0 [0.5 1.0] zero velocity
        # x2 = 1.5 [1.0 1.5] positive div(u)
        # x3 = 2.0 [1.5 2.0] enforced u,v,w        
        
        nfu = domain.new_field() # Target u
        nfd = domain.new_field() # Divergence
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
        
        nwu['g'] = ju.cubici( x, x0, x1) - ju.cubici( x, x2, x3)   
        nfu['g'] = fx1*wr 

        nfd['g'] = fd1     + U2*fd2
        TQ1 = dV*ju.jsum(fd1)/IA
        TQ2 = dV*U2*ju.jsum(fd2)/IA
        TQ3 = dV*ju.jsum(nfd['g'])/IA
        
        if Ny>1: problem.substitutions['wv']   = 'wu'
        problem.substitutions['ww']   = 'wu'
    elif forcing==1:
        TB   = IA
        Q    = IA
        ifw  = IA/(L*W)
    elif forcing==2:  # stream function forcing between two gaussians
        # periodic-gaussian.mw
        # x0 = 0.5 [0.0 0.5] gaussian 1
        # x1 = 1.0 [0.5 1.0] switch region
        # x2 = 1.5 [1.0 1.5] switch region
        # x3 = 2.0 [1.5 2.0] gaussian 2         

        nfu = domain.new_field() # Target u
        if Ny>1: nfv = domain.new_field() # Target v
        nfw = domain.new_field() # Target w

        #      0
        # x0 = 0.5
        # x1 = 1.0
        # x2 = 1.5
        # x3 = 2.0
        fx1,Ffx1,dfx1=ju.smoothw(x,nfx,x0,x2)
        fx2,Ffx2,dfx2=ju.smoothw(x,nfx,x3,L)
        fx=fx1-fx2
        bb=0.5
        dfx=(dfx1-dfx2)/R*(1+bb)/(1+bb*fx)**2
        #ib=(1+bb)*fx/(1+bb*fx)/R
        ib=fx
 
 
        U1=1;
        U2=math.pi*R**2/(H*W)
        
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


        nfu = domain.new_field() # Target u
        nfv = domain.new_field() # Target v
        nfw = domain.new_field() # Target w
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
        

        nn=5
        if forcing==3:   nfu['g'],nfv['g'],nfw['g'] = ju.weighted_gaussian(c,dc,y,z,h,1,W,H,nn)
        else:            nfu['g'],nfv['g'],nfw['g'] = jtheta.periodic_gaussian(c,dc,y,z,1,W,H,R)
        
        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'

        if forcing==3: fu,fv,fw = ju.weighted_gaussian(c0,0,y,z,h,1,W,H,nn) # Forcing function at radius
        else:          fu,fv,fw = jtheta.periodic_gaussian(c0,0,y,z,1,W,H,R)
 
        #ifu=nfu['g']
        #ifv=nfv['g']
        #ifw=nfw['g']


        del(fu,fv,fw,fx,dfx,c,c0,dc,nn)


if sType=="bg": 
    if xa is not None and xb is not None: ifu = ju.stepwf(x,nfx,xa-Wx/2,xa+Wx/2,xb-Wx/2,xb+Wx/2,L)

if sType=='gc':
    if U1 is None: U1=0
    hu=min(H,max(0,hu))
    RL=hu
    if forcing is not None:
        # This is what the forcing should achieve at the inlet
        problem.parameters['EFB']=0
        problem.parameters['EFS']=0
        problem.parameters['EFU']=-param['U']
        problem.parameters['EFV']=0
        problem.parameters['EFW']=0
        problem.parameters['EFN']=Nxx-1

        U0=-1
        #ju.logger.info('H {} hu {} U {} U0 {} U1 {}'.format(H,hu,U,U0,U1))
        if   hu <= 0: U2 = U0
        elif hu >= H: U2 = 0.0
        else:         U2 = (H*U0-hu*U1)/(H-hu);   # h*U1+(H-h)*U2=H*U0
    else:
        U0 = None
        U1 = None
        U2 = None
    if noise is not None:
        if xn4>xn1:
            wnoise = domain.new_field()  # weighting for streamfunction
            if xn1==xn2 and xn3==xn4: 
                wnoise['g'] = ju.cubich(x,xn1,xn4)
                ju.logger.info('wnoise: cubich {:5.2f}-{:5.2f}'.format(xn1,xn4))
            else:
                wnoise['g'] = ju.stepwf(nwx,xn1,xn2,xn3,xn4,L)
                del(TTWN,TTTWN)
                ju.logger.info('wnoise: stepw {:5.2f} {:5.2f} {:5.2f} {:5.2f}'.format(xn1,xn2,xn3,xn4))
                
            if Tx=='SinCos': wnoise.meta['x']['parity'] = 1
        else:
            wnoise=1
            ju.logger.info('wnoise: noise everywhere')
 
    # All the buoyancy forcing is the same controlled by x4,x5,x6,x7 and x0 and Wx
    # weight region is x0,x0+Wx to x6,x7
    # b=B    region is x4,x5 past x7 and 0 to hb+Wz/2,hb-Wz/2
    if hb is not None:
        if hb>=H: hbz=np.ones_like(z)
        else:     hbz=ju.smoothwf(z,nfz,hb+Wz/2,hb-Wz/2)
        if hb>=H: hbzz=np.ones_like(zz)
        else:     hbzz=ju.smoothwf(zz,nfz,hb+Wz/2,hb-Wz/2)
    if xb is not None:
        ifb=domain.new_field()
        if Tx == "SinCos": ifb.meta['x']['parity'] = 1
        if Ty == "SinCos": ifb.meta['y']['parity'] = 1
        if Tz == "SinCos": ifb.meta['z']['parity'] = 1
        ifb.set_scales(AA)
        if xa is not None: ifb['g'] = hbzz*(ju.ffta_heaviside(xx-xa,maxkx,nfx)[0]-ju.ffta_heaviside(xx-xb,maxkx,nfx)[0])
        else:              ifb['g'] = hbzz*ju.ffta_heaviside(xb-xx,maxkx,nfx)[0]
    if Tz=='Cheb':
        if param['lbc']=='slip':
            if   hu<=0: hu1=np.zeros_like(z)
            elif hu>=H: hu1=np.ones_like(z)
            else:       hu1=ju.smoothwf(z,nfz,hu+Wz/2,hu-Wz/2)
        else:
            if   hu<=0: hu1=np.zeros_like(z)
            elif hu>=H: hu1=ju.smoothwf(z,nfz,0-Wz/2,0+Wz/2)
            else:       hu1=ju.smoothwf(z,nfz,-Wz/2,Wz/2)-ju.smoothwf(z,nfz,hu-Wz/2,hu+Wz/2)
        hu2=1-hu1
    if forcing is not None:
        if forcing!=7:
            W76=x7-x6
            nwb['g'] =       ju.stepwf(x,nfx,x0,x0+Wx,x6,x7,      L) # Weight field
            nfb['g'] = hbz*ju.stepwf(x,nfx,x4,x5,   x7+W76,x7+2*W76,L) # Buoyancy field
    else:
        if dimple is not None:
            #ju.write_hdf5('kx.hdf5','kx',domain.bases[0].wavenumbers) #kx = (0:Nx-1)*pi/L
            #ju.write_hdf5('ky.hdf5','ky',domain.bases[1].wavenumbers) #kx = (0:Nx-1)*pi/L
            #ju.write_hdf5('kz.hdf5','kz',domain.bases[1].elements)
            dmslice = domain.dist.coeff_layout.slices(scales=1)
            dmshape = domain.dist.coeff_layout.local_shape(scales=1)
            nyz=domain.new_field()
            if Tx == "SinCos": nyz.meta['x']['parity']=1
            if Ty == "SinCos": nyz.meta['y']['parity']=1
            dmky     = domain.bases[1].wavenumbers[np.newaxis,dmslice[1],np.newaxis]
            dmkz     = domain.bases[2].elements[   np.newaxis,np.newaxis,dmslice[2]]*math.pi/H
            #ju.logger.info('dmky [{:8.3f} {:8.3f}]'.format(ju.jmin(dmky),ju.jmax(dmky)))
            #ju.logger.info('dmkz [{:8.3f} {:8.3f}]'.format(ju.jmin(dmkz),ju.jmax(dmkz)))
            nyz['c'] = dmky*dmkz*np.exp(-(dimplewy*dmky+dimplewz*dmkz))*rand.standard_normal(dmshape)
            nyz['c'][1:-1,:,:]=0
            nyz.set_scales(1)
            dmshape = domain.dist.grid_layout.local_shape(scales=1)
            dmshape[0]=1
            nyz=nyz['g'][0,:,:]
            nyz=np.reshape(dimple*nyz/np.sqrt(ju.jsum(nyz**2)/(Ny*Nz)),dmshape)
            ifb['g']=hbz*ju.smoothwf(xb+nyz-x,nfx,-Wx/2,Wx/2)
       

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
        
        if Ny>1: nwv['g']  = w4
        nww['g']  = (w1 + w2)
        nwu['g']  = (w1 + w2)
        nfu['g']  = U0*(1-W2) + W2*(U1*hu1+U2*hu2)
    elif forcing==2:
        w1=ju.indicatorp(x,x0,   x1,   Dx,L)   # This region force uniform velocity profile
        W1=ju.indicatorp(x,x0-Wx,x1+Wx,Dx,L) 
        nwu = domain.new_field()  # Weight for u
        fu        = W1
        
        if Ny>1: nwv['g']  = w1
        
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
        xwu = ju.stepwf(x,nfx,x0,x0+Wx,x2,x3      ,L)  # velocity weight
        xfu = ju.stepwf(x,nfx,x1,   x2,x3,x3+x2-x1,L)  # velocity value
        
        
        nfu['g'] = (xfu*(U1*hu1+U2*hu2)-(1-xfu))
        nfu.differentiate('x', out=nfd)
        nfd['g']=np.where((xx-x0)*(x3-xx)>0,nfd['g'],0)
        nwu['g'] = xwu


        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'
        if xb is not None:
            if xb>x2: xfu=ju.stepwf(x,nfx,x1,x2,xb-Wx,xb+Wx,L)
        #ifu=(xfu*(U1*hu1+U2*hu2)-(1-xfu))
    
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
    
        nwu['g'] =  ju.stepwf(x,nfx,x0,x0+Wx,x1,x2,L)     # Velocity   weight
 
        #hu*U1+(H-hu)*U2=-H*U
        U2=-(H+hu*U1)/(H-hu)
        dw=ju.dcubici(x,x4,x3)  
        nfd['g'] = -dw*(1+U1*hu1+U2*hu2)
  
        if PIDX is not None and not PIDG:
            nfu = domain.new_field()  # Velocity function
            nfu['g'] = -1
        else:
            problem.substitutions['fu'] = '-U'
 


        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'
        
 
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
          
        xfu =ju.stepwf(x,nfx,x0,x1,x2,x3,L) # Velocity forcing weight
   
        nwu['g'] = xfu


        if Ny>1: problem.substitutions['wv'] = 'wu'

        problem.substitutions['ww'] = 'wu'

        if PIDX is not None and not PIDG:
            nfu = domain.new_field()  # Velocity function
            nfu['g']=-1
        else: problem.substitutions['fu'] = '-U'
        
        if xb is not None:
            if xb>x1: xfu=ju.stepwf(x,nfx,x0,x1,xb-Wx,xb+Wx,L)
        #ifu=xfu*(U1*hu1+U2*hu2)+(1-xfu)*U0
         
        
    elif forcing==6:
        # Streamfunction forcing with velocity controlled by x0,x1,x2 and x3
        # buoyancy forcing controlled by x4, x5, x6 and x7

        nwu  = domain.new_field()  # Velocity weight 
        psi  = domain.new_field()  # stream function
        
        if Ty=="SinCos": psi.meta['y']['parity'] = 1
        fx,Fx,dfx   = ju.stepw(x,nfx,x1,x3,   x3,x3+x3-x1,L)
        xwu,Iwu,Dwu = ju.stepw(x,nfx,x0,x0+Wx,x2,x3,      L)
        nwu['g'] = xwu   # Velocity weight
    
        fx1,Fx1,dfx1 = ju.smoothw(H,nfz,hu-Wz/2,hu+Wz/2)
        fx0,Fx0,dfx0 = ju.smoothw(0,nfz,hu-Wz/2,hu+Wz/2)
        Fx1-=Fx0
        FF = (U1-U0)*H-U1*Fx1
        F2 = Fx1

        
        if param['lbc']=='noslip':
            fx1,Fx1,dfx1 = ju.smoothw(H,nfz,-Wz/2,+Wz/2)
            fx0,Fx0,dfx0 = ju.smoothw(0,nfz,-Wz/2,+Wz/2)
            Fx1-=Fx0
            FF += 2*(U0-U1)*(H-Fx1)
        if param['rbc']=='noslip':
            fx1,Fx1,dfx1 = ju.smoothw(H,nfz,H-Wz/2,H+Wz/2)
            fx0,Fx0,dfx0 = ju.smoothw(0,nfz,H-Wz/2,H+Wz/2)
            Fx1-=Fx0
            FF +=  2*Fx1*U0
            F2 += -2*Fx1
        #U0*H + fx*FF + fx*F2*U2 = U0*H -> FF + F2*U2=0
        U2=copy.copy(-FF/F2)
        
        fx1,Fx1,dfx1 = ju.smoothw(z,nfz,hu-Wz/2,hu+Wz/2)
        fx0,Fx0,dfx0 = ju.smoothw(0,nfz,hu-Wz/2,hu+Wz/2)
        Fx1-=Fx0
        FF = (U1-U0)*z+Fx1*(U2-U1)
        if param['lbc']=='noslip':
            fx1,Fx1,dfx1 = ju.smoothw(z,nfz,-Wz/2,+Wz/2)
            fx0,Fx0,dfx0 = ju.smoothw(0,nfz,-Wz/2,+Wz/2)
            Fx1-=Fx0
            FF += 2*(U0-U1)*(z-Fx1)
        if param['rbc']=='noslip':
            fx1,Fx1,dfx1 = ju.smoothw(z,nfz,H-Wz/2,H+Wz/2)
            fx0,Fx0,dfx0 = ju.smoothw(0,nfz,H-Wz/2,H+Wz/2)
            Fx1-=Fx0
            FF += 2*(U0-U2)*Fx1
        psi['g'] = U0*z + fx*FF

        nfu = psi.differentiate('z')  
        nfw = psi.differentiate('x')
        nfw['g'] = -nfw['g']

        

        problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'

        
        if xb is not None:
            if xb>x3: xfu=ju.stepw(x,nfx,x1,x3,xb-Wx/2,xb+Wx/2,L)
        #ifu=xwu*(U1*hu1+U2*hu2)-(1-xwu)
        
        del(fx1,Fx1,dfx1,fx,Fx,Iwu,Dwu,psi)
    elif forcing==7: # Forces only by adding and removing fluid Tx should be sincos
        nfd = domain.new_field()  # div(u)
        nfd.set_scales(AA)
        nfd.meta['x']['parity']=1

        
        #[u1,wux4]=ju.ffta_heaviside(xx,maxkx/f7s,Anfx)[0:2]
        #[wux4,dwux4,I4,u1]=ju.ffta_heaviside(f7s*xmaxkx-xx,maxkx/f7s,nfx)
        #ju.logger.info('wux4 {}'.format(wux4))
        #ju.logger.info('f7s*xmaxkx {:7.4f}'.format(f7s*xmaxkx))
        #ju.logger.info('I4: {:7.4f}, Iwux4a: {:7.4f}, Iwux4: {}'.format(I4,Iwux4a,Iwux4[-1]))
        [wux4,uu]=ju.ffta_df(xx,maxkx/2,7)
        Iwux4=dxx*wux4.sum()
        wux4=wux4/Iwux4
        uu=uu/Iwux4
        rful=uu*wux4
        wux5=np.flip(wux4)
        rfur=np.flip(-rful)
        wux4=np.reshape(wux4,szxx)
        wux5=np.reshape(wux5,szxx)
        if   hu==0: fuzz=np.ones_like(zz)
        elif hu==H: fuzz=np.zeros_like(zz)
        else:       fuzz=ju.smoothwf(zz,nfz,hu-Wz,hu+Wz)/(1-hu/H)
        if not gclf: rful=None
        if not gcrf: rfur=None

        Ndf=7
        math.pi*Ndf/hb
        #kwux3=maxkx/2
        M = 7
        MM=np.int64(np.floor(hb/dxx))


        kwux3=2*M*math.pi/hb
        wux3=ju.ffta_df(xx,kwux3,M)[0]      # Fix density profile
        wux3=wux3/wux3.max()
        ju.logger.info('MM: {}, xx: [{} {}], w [{},{}]'.format(MM, xx[MM], xx[Nxx-MM],wux3[MM-1],wux3[MM]))
        if db:
            wux2=ju.ffta_df(xx-hb/2,kwux3,M)[0]      # Sample density profile
            wux2=wux2/(wux2.sum()*dxx)
            gbp            = np.zeros((Nzz,), dtype=np.float64)
            state_ev['bp'] = np.zeros((Nzz,), dtype=np.float64)
            #print('hbzz: {}'.format(hbzz))
            if yycolor:
                yycomm.Gatherv(sendbuf=hbzz, recvbuf=(state_ev['bp'], zzcounts), root=0)
                #print('rank: {} yycomm.rank {} state_ev[bp] {} zzcounts {}'.format(ju.mpirank,yycomm.Get_rank(),state_ev['bp'], zzcounts))
            sbp=state_ev['bp']
            state_ev['bp'] = ju.comm.bcast(state_ev['bp'],0)
            hbzz=((state_ev['bp'][zzslices]).reshape(szzzl))
        if fbmult:
            
            bmultls = slice(0, MM)
            bmultl  = wux3[bmultls]

            wux8=ju.ffta_df(L-xx,maxkx  ,M)[0]   # Top buoyancy forcing region
            wux8=wux8/wux8.max()
            MM=np.int64(np.ceil(M*math.pi/maxkx/dxx))
            bmultrs = slice(Nxx-MM, Nxx)
            bmultr  = np.minimum(1,np.maximum(0,wux8[bmultrs]))
            

            if PIDX is not None:
                if PIDX >= xx[Nxx-MM]:
                    ju.logger.info('PIDX {:6.3f} is greater than forcing point {}'.format(PIDX, xx[Nxx-MM]))
                    ju.set_error('PIDX {:6.3f} is greater than foring point {}'.format(PIDX, xx[Nxx-MM]))
            ju.check_error()
                
        else:
            if not fbmax:
                nwb.set_scales(AA)
                nfb.set_scales(AA)
            if not fbmin: nwb.set_scales(AA)

            if fbmax: bmaxf = wux3*hbzz
            else:
                nwb['g'] += wux3
                if fbmin: nfb['g'] = hbzz * ju.ffta_heaviside(L/2-xx,maxkx/2,nfx)[0]
                else:     nfb['g'] = hbzz * wux3
                if fbmin:  bmin  = ju.ffta_heaviside(L-0.5-xx,maxkx/2,nfx)[0]
                else:  nwb['g'] += np.flip(wux3)

        if not ddiv:
            wux1=wux4.flatten()
            nfd['g']  = wux5 - wux4*fuzz
        else:
            xxf=xx.flatten()
            wux1 = wux4.flatten()
            wux1 = wux1/(wux1.sum()*dxx)
            Iwux1=np.sqrt((wux1**2).sum()*dxx)

            if pwu or pwv or pwwl or pwwr:
                nww = domain.new_field()  # div(u)
                ju.set_parity(nww,Tx,Ty,Tz,1,1,1)
                nww.set_scales(AA)
                nww['g']=0
                if pwwl: nww['g'] += wux1.reshape(szxxl)/wux1.max()
                if pwwr: nww['g'] += np.flip(wux1).reshape(szxxl)/wux1.max()
                nww['g'] = ju.clip01(nww['g'])
                if pwwl or pwwr:
                    if pwu: problem.substitutions['wu'] = 'ww'
                    if pwv: problem.substitutions['wv'] = 'ww'
                elif pwu:
                    nwu=nww
                    del(nww)
                    if pwv: problem.substitutions['wv'] = 'wu'
                
            
            if gclf or gcrf:
                rfu=domain.new_field() # Forcing field
                ju.set_parity(rfu,Tx,Ty,Tz,-1,1,1)
                rfu.set_scales(AA)
            if gclf:
                if ddiv:
                    rfu['g']=wux1.reshape(szxx)
                    if dim==2: rful = (rfu.integrate('x')['g']*rfu['g'])[:,0]
                    else:      rful = (rfu.integrate('x')['g']*rfu['g'])[:,0,0]
                    
            if gclf or gcrf: problem.parameters['rfu']=rfu
            if ddiv:
                gdw              = np.zeros((Nzz,), dtype=np.float64)
                state_ev['divz'] = np.zeros((Nzz,), dtype=np.float64)
                #if yycolor: yycomm.Gatherv(sendbuf=fuzz, recvbuf=(state_ev['divz'][:,0], zzcounts), root=0)
            state_ev['divz']=-ju.smoothwf(gzz.flatten(),nfz,hu-Wz,hu+Wz)/(1-hu/H)
            ju.logger.info('wux5.shape {} wux1.shape {} divz.shape {} szzz {}'.format(wux5.shape,wux1.shape,state_ev['divz'].shape,szzzl))
            nfd['g']    = wux5 + wux1.reshape(szxx)*state_ev['divz'][zzslices].reshape(szzzl)
            if   gclf and gcrf: rfu['g'] = U*(rfur + rful.reshape(szxx)*state_ev['divz'][zzslices].reshape(szzzl)**2)
            elif gclf:          rfu['g'] = U*(       rful.reshape(szxx)*state_ev['divz'][zzslices].reshape(szzzl)**2)
            elif gcrf:          rfu['g'] = U*(rfur)
            divz0=state_ev['divz']
            del(wux4)
            
 
        
if fuspecial is not None:
    problem.parameters['wu'] = 1
    problem.parameters['wv'] = 1
    problem.parameters['ww'] = 1
    problem.parameters['Fu']=param['force']
    problem.parameters['Fv']=param['force']
    problem.parameters['Fw']=param['force']

    nfu=domain.new_field()
    nfw=domain.new_field()
    if fuspecial=='sinsin':
        nfu['g'] = +np.sin(fuspecial1*2*math.pi/L*x)*np.cos(2*math.pi/H*z)
        nfw['g'] = -np.cos(fuspecial2*2*math.pi/L*x)*np.sin(2*math.pi/H*z)
        ju.logger.info('U {}'.format(U))
if PIDG:
    nfg = domain.new_field()  # gravity function
    nfg['g'] = g
    nfg.meta['x']['constant']=True
    if dim==3: nfg.meta['y']['constant']=True
    nfg.meta['x']['constant']=True
    ju.set_parity(nfg,Tx,Ty,Tz,1,1,1)
    problem.parameters['g'] = nfg
elif g is not None: problem.parameters['g'] = g

if nfd is not None:
    if Tx=="SinCos": nfd.meta['x']['parity'] = 1
    if Ty=="SinCos": nfd.meta['y']['parity'] = 1
    if Tz=="SinCos": nfd.meta['z']['parity'] = 1
if nfu is not None:
    if Tx=="SinCos": nfu.meta['x']['parity'] = 1
    if Ty=="SinCos": nfu.meta['y']['parity'] = 1
    if Tz=="SinCos": nfu.meta['z']['parity'] = 1
if rfu is not None:
    if Tx=="SinCos": rfu.meta['x']['parity'] = -1
    if Ty=="SinCos": rfu.meta['y']['parity'] = 1
    if Tz=="SinCos": rfu.meta['z']['parity'] = 1

ju.check_error()


if noise is not None: problem.parameters['fnoise'] = zmn


if S is not None:
    minS=0
    maxS=SM
else:
    minS=None
    maxS=None
    clipS=False
if B is not None:
    minB=0
    maxB=BM
else:
    minB=None
    maxB=None
    clipB=False
   
if 'minB' in param: minB=param['minB']
if 'maxB' in param: maxB=param['maxB']
if 'minS' in param: minS=param['minS']
if 'maxS' in param: maxS=param['maxS']


if forcing==0 and noise is not None:
    problem.parameters['wnoise']=wnoise
    EWN=1
else:
    EWN = ju.set_wfun(problem, wnoise, 'wnoise', L*H*W, Tx, Ty, Tz, AA, JA, AAJ)
if EWN is not None: problem.parameters['EWN']=EWN
del(EWN)
 
ju.set_wfun(problem, nwd,    'wd'    , L*H*W, Tx, Ty, Tz, AA, JA, AAJ)
ju.set_wfun(problem, nwb,    'wb'    , L*H*W, Tx, Ty, Tz, AA, JA, AAJ)
ju.set_wfun(problem, nws,    'ws'    , L*H*W, Tx, Ty, Tz, AA, JA, AAJ)
ju.set_wfun(problem, nwu,    'wu'    , L*H*W, Tx, Ty, Tz, AA, JA, AAJ)
ju.set_wfun(problem, nwv,    'wv'    , L*H*W, Tx, Ty, Tz, AA, JA, AAJ)
ju.set_wfun(problem, nww,    'ww'    , L*H*W, Tx, Ty, Tz, AA, JA, AAJ)

#if fpidu and nfu is not None:
#    ju.logger.info('Setting fpidu to xwu*(-(udxx+udzz)/Re + nfu*udx+nfw*udz)')
#    udx   = nfu.differentiate('x')
#    udxx  = udx.differentiate('x')
#    udz   = nfu.differentiate('z')
#    udzz  = udz.differentiate('z')
#    rfu   = U*xwu*(-(udxx+udzz)/Re + nfu*udx+nfw*udz)
#    del(udx,udxx,udz,udzz)

if fpid is not None: problem.parameters['fpid']=fpid
#if param['fpidw'] and nfw is not None:
#    ju.logger.info('Setting fpidw to nwu*(-(wdxx+wdzz)/Re + nfu*wdx+nfw*wdz)')
#    wdx   = nfw.differentiate('x')
#    wdxx  = wdx.differentiate('x')
#    wdz   = nfw.differentiate('z')
#    wdzz  = wdz.differentiate('z')
#    rfw   = U*xww*(-(wdxx+wdzz)/Re + nfu*wdx+nfw*wdz)
#    del(wdx,wdxx,wdz,wdzz)


if nfd: Tdiv=ju.jint(nfd)
else:   Tdiv=None

BSflg = 'pmss' in param
Uflg  = PIDX is None or PIDG or sType=='pm'
   
fb  = ju.set_ffun(problem,nfb, 'fb',1,Tx,Ty,Tz,BSflg,AA, 1,JA,AAJ)
fs  = ju.set_ffun(problem,nfs, 'fs',1,Tx,Ty,Tz,BSflg,AA, 1,JA,AAJ)
fd  = ju.set_ffun(problem,nfd, 'fd',U,Tx,Ty,Tz,Uflg | BSflg, AA, 1,JA,AAJ)
fu  = ju.set_ffun(problem,nfu, 'fu',U,Tx,Ty,Tz,Uflg, AA,-1,JA,AAJ)
fv  = ju.set_ffun(problem,nfv, 'fv',U,Tx,Ty,Tz,Uflg, AA, 1,JA,AAJ)
fw  = ju.set_ffun(problem,nfw, 'fw',U,Tx,Ty,Tz,Uflg, AA, 1,JA,AAJ)
#rfu = ju.set_ffun(problem,rfu,'rfu',U,Tx,Ty,Tz,Uflg, AA,-1,JA,AAJ)
if rfu is not None: JA.savexyz('force/rfu.hdf5','rfu',rfu['g'],AA)
if rfv is not None: JA.savexyz('force/rfv.hdf5','rfv',rfv['g'],AA)
if rfw is not None: JA.savexyz('force/rfw.hdf5','rfw',rfw['g'],AA)

 
    
ju.open_stats('stats.hdf5',param,reset,problem)

if restart:
    if rfn is not None: fnrst=ju.find_start_files(bdir / Path(rfn))
    else:               fnrst=ju.find_start_files(bdir / name)
    if len(fnrst)>0:    ju.logger.info('Restart file: {}'.format(fnrst[0]))

if (not restart) or start_new_files:
    mode = "overwrite"
else:
    mode = "append"
    
if sType=='gc':
    problem.parameters['A']  = W
    fluxd='yz'
elif sType=='pm':
    problem.parameters['A']  = 1
    if   forcing==5: fluxd='yz'
    elif forcing==1: fluxd='xy'
    elif forcing==2: fluxd='yz'
    elif forcing==3: fluxd='yz'
    elif forcing==4: fluxd='yz'
    elif forcing==6: fluxd='yz'
    elif forcing==7: fluxd='yz'
    else:
        ju.logger.error('forcing should be 0 to 7 for plumes'.format(forcing))
        ju.set_error('Unknown forcing for plume')
elif sType=='bg':
    problem.parameters['A']  = 1
    if forcing==0: ju.logger.info('Burgers equation with no forcing')
else:
    ju.logger.error('sType should be gc, pm or bg'.format(sType))
    ju.set_error('Unknown simulation type')
ju.check_error()

#https://docs.python.org/3/library/resource.html
#resource.RLIMIT_AS   The maximum area (in bytes) of address space which may be taken by the process.
#resource.RLIMIT_VMEM The largest area of mapped memory which the process may occupy.
if VMEM is not None:
    memlim=np.int64(VMEM)*1024**3
    resource.setrlimit(resource.RLIMIT_AS,  (memlim,memlim))
    resource.setrlimit(resource.RLIMIT_RSS, (memlim,memlim))
    ju.logger.info('Set maximum memory: {:6.3f} GB (AS & RSS)'.format(VMEM))
 
smem=ju.MemStats()
#ju.logger.info(smem)

ju.logger.info("Running ded_gc.py: Bouyancy driven flow simulation system")
ju.run_info()
ju.logger.info("signals:            {}".format(signals))
ju.logger.info("parallel:           {}".format(parallel))

def snn(f,x):
    if x is not None: return(f.format(x))
    else: return('')

if forcing is not None:
    ju.logger.info("forcing:  {:1d}, force:    {:4.0f}, forceb:    {:4.2f}, forces:    {:4.2f}".format(forcing,F,param['forceb'],param['forces']))
    ju.logger.info('nfx: {:2d}, xmaxkx1 {:6.4f}, xmaxkx2 {:6.4f}'.format(nfx,xmaxkx,xmaxkx2))
   
    ju.logger.info(snn('x0: {:7.4f} ',x0)+snn('x1: {:7.4f} ',x1)+snn('x2: {:7.4f} ',x2)+snn('x3: {:7.4f} ',x3));
    ju.logger.info(snn('x4: {:7.4f} ',x4)+snn('x5: {:7.4f} ',x5)+snn('x6: {:7.4f} ',x6)+snn('x7: {:7.4f} ',x7));
if clipB:   ju.logger.info('clipB: BM:{:5.2f} [{:5.2f}, {:5.2f}]'.format(BM,minB,maxB))
if clipS:   ju.logger.info('clipS: SM:{:5.2f} [{:5.2f}, {:5.2f}]'.format(SM,minS,maxS))
if clipu is not None: ju.logger.info("clipu: {:4.0f},".format(clipu))
ju.logger.info("restart:            {}".format(restart))
ju.logger.info("reset:              {}".format(reset))
ju.logger.info("start_new_files:    {}".format(start_new_files))
if dtcheck is not None:  ju.logger.info("dtcheck:            {:7.4f} (min)".format(dtcheck))
if dtjcheck>0: ju.logger.info("dtjcheck:           {:7.4f} (min)".format(dtjcheck))
if avrg: ju.logger.info("avrg:               {}".format(avrg))
ju.logger.info("Wall time:          {}".format(param['time']))
ju.logger.info("Integration scheme: {}".format(scheme))
ju.logger.info("mindt:  {:8.6f}, maxdt: {:8.6f}".format(mindt, maxdt))
ju.logger.info(ju.nns( 'dx:    {:7.5f}',dx)    + ju.nns( ',  dy:    {:7.5f}',dy)    + ju.nns( ',  dz:    {:7.5f}',dz)    + ju.nns( ',  dr:    {:7.5f}',dr)) 
ju.logger.info(ju.nns( 'Wx:    {:7.5f}',Wx)    + ju.nns( ',  Wy:    {:7.5f}',Wy)    + ju.nns( ',  Wz:    {:7.5f}',Wz)    + ju.nns( ',  Wr:    {:7.5f}',Wr)) 
ju.logger.info(ju.nnsd('Wx/dx: {:7.5f}',Wx,dx) + ju.nnsd(',  Wy/dy: {:7.5f}',Wy,dy) + ju.nnsd(',  Wz/dz: {:7.5f}',Wz,dz) + ju.nnsd(',  Wr/dr: {:7.5f}',Wr,dr)) 
if Tdiv: ju.logger.info("Tdiv:  {:7.5e}".format(Tdiv))
ju.log_param(iparam,('Re','T'),'{:8.0f}')
ju.log_param(iparam,('B','S','SV','Skp','Ski','Skg','Scb','Scs','inoise','bnoise','snoise','sinitial','hu','hb','U','g','q','pmss','pmsl'),'{:8.5f}')
ju.log_param(iparam,('lbc','rbc','blbc','brbc','slbc','srbc'),'{:>7s}')
ju.log_param(iparam,('xa','xb','x0','x1','x2','x3','x4','x5','x6','x7'),'{:8.5f}')
ju.log_param(iparam,('AA','AAJ','AAS'),'{:4.2f}')
ju.log_param(iparam,('pmzt','pmzr','pmzrt','oddg'),'{}')

if fbmult: ju.logger.info("Minimum bouyancy enforced by multiplication")
if fbmin:  ju.logger.info("Minimum bouyancy enforced by min")
if fbmax:  ju.logger.info("Maximum bouyancy enforced by max")
if fsmult: ju.logger.info("Minimum bouyancy enforced by multiplication")
if fsmin:  ju.logger.info("Minimum sediment enforced by min")
if fsmax:  ju.logger.info("Maximum sediment enforced by max")
if ddiv:   ju.logger.info("Dynamic divergence")
if db:     ju.logger.info("Dynamic bouyancy")
if noise is not None:
    ju.logger.info("noise:  {:7.4f}, noiseL: {:7.5f}, noiseT: {:7.5f}, noised {}".format(noise,noiseL,noiseT,noised))
    ju.logger.info("xn1:    {:7.4f}, xn2:     {:7.4f}, xn3:    {:7.4f},  xn4: {:7.4f}".format(xn1,xn2,xn3,xn4))
if U1 is not None: ju.logger.info("U1:     {:7.4f}".format(U1))
if U2 is not None: ju.logger.info("U2:     {:7.4f}".format(U2))
if pwu:                ju.logger.info("Forcing u velocity")
if pwv:                ju.logger.info("Forcing v velocity")
if pwwl:               ju.logger.info('Forcing w velocity left')
if pwwr:               ju.logger.info('Forcing w velocity right')
if rescale:            ju.logger.info('Rescaling velocities and pressure')

if sType=='pm':
    ju.logger.info("Inlet shape:    {:}".format(inlet))
    ju.logger.info("Inlet radius:   {:7.4f}".format(R))
    ju.logger.info("Inlet velocity: {:7.4f}".format(U))
    ju.logger.info("Inlet area           IA: {:7.4f}".format(IA))
    ju.logger.info("       sqrt(IA/pi): {:7.4f}".format(np.sqrt(IA/np.pi)))
    ju.logger.info("Inlet flux:         {:7.4f} {:7.4f}".format(IQ,IIQ))
    ju.logger.info("Inlet flux squared: {:7.4f} {:7.4f}".format(IQQ,IIQQ))
    if forcing==0: ju.logger.info("TQ1:   {:7.4f} TQ2: {:7.4f} TQ3: {:7.4e}".format(TQ1,TQ2,TQ3))
if pmzt:  ju.logger.info('Buoyancy is set to zero at top');
if pmzr:  ju.logger.info('Buoyancy is set to zero at edges');
if pmzrt: ju.logger.info('Buoyancy is set to zero at top edges');

#for s in problem.substitutions: ju.logger.info('substitutions: {:7s}: {}'.format(s,problem.substitutions[s]))
#for s in problem.parameters:    ju.logger.info('parameters:    {:7s}: {}'.format(s,problem.parameters[s]))

ju.print_parity(problem)
if sType=='bg': ju.burgers_equations(domain,problem,param)
else:           ju.boussinseq_equations(domain,problem,param)


#problem.meta[:]['z']['dirichlet'] = True
if Tz=='Cheb':
    if noise is not None: # If we are using this as a stream function must be constant at the ends 
        problem.add_bc("right(noisef1) = 0")
        problem.add_bc("left(noisef1)  = 0")
    if sType=='pm':
        problem.add_bc("left(b)    = left(BF(t,x,y,h,B))*(1+left(noisef1))")
        problem.add_bc("right(bz)  = 0")
        problem.add_bc("left(v)    = 0")
        problem.add_bc("right(vz)  = 0")
        problem.add_bc("left(u)    = 0")
        problem.add_bc("right(uz)  = 0")
        problem.add_bc("left(w)    = left(BF(t,x,y,h,U))*(1+left(noisef1))")
        problem.add_bc("right(wz)  = 0", condition="(nx != 0) or (ny != 0)")
        problem.add_bc("integ_z(p) = 0", condition="(nx == 0) and (ny == 0)")
    else:
        problem.meta['w']['z']['dirichlet'] = True
        #problem.meta['u']['z']['dirichlet'] = True
        #problem.meta['uz']['z']['dirihlet'] = True
        #problem.meta['w']['z']['dirichlet'] = True
        #problem.meta['wz']['z']['dirichlet'] = True
        #problem.meta['b']['z']['dirichlet'] = True
        #problem.meta['bz']['z']['dirichlet'] = True
        if   param['lbc']=='slip':
            problem.add_bc(         "left(uz)   = 0")
            if Ny>1: problem.add_bc("left(vz)   = 0")
        elif param['lbc']=='noslip':
            if PIDX is not None and not PIDG: problem.add_bc("left(u) = left(fu)")
            elif U is not None:     problem.add_bc("left(u) = -U")
            else:                   problem.add_bc("left(u) = 0")
            if Ny>1: problem.add_bc("left(v)    = 0")
        if   param['rbc']=='slip':
            problem.add_bc(         "right(uz)  = 0")
            if Ny>1: problem.add_bc("right(vz)  = 0")
        elif param['rbc']=='noslip':
            if PIDX is not None and not PIDG: problem.add_bc("right(u) = right(fu)")
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
                
        vz=list(set(problem.variables) & set(('bz','fpidbz','fpidsz','fpiduz','fpiduz','fpidwz')))
        for v in vz:
            if v != 'sz' and v!= 'bz':
                problem.add_bc("right(" + v + ") = 0")
                problem.add_bc("left("  + v + ") = 0")
            problem.meta[v]['z']['dirichlet'] = True

        if 'sz' in problem.variables:
            if   param['srbc']=='z':  problem.add_bc("right(s)  = 0")
            elif param['srbc']=='dz': problem.add_bc("right(dzs)  = 0")
            elif param['srbc']=='dq': problem.add_bc("right(kappas*dzs+SVZ*s= 0")
            if   param['slbc']=='z':  problem.add_bc("left(s)  = 0")
            elif param['slbc']=='dz': problem.add_bc("left(dzs)  = 0")
            elif param['slbc']=='dq': problem.add_bc("left(kappas*dzs+SVZ*s= 0")

        if 'bz' in problem.variables:
            if   param['brbc']=='z':  problem.add_bc("right(b)  = 0")
            elif param['brbc']=='dz': problem.add_bc("right(dzb)  = 0")
            elif param['brbc']=='dq': problem.add_bc("right(kappab*dzb+bVZ*b= 0")
            if   param['blbc']=='z':  problem.add_bc("left(b)  = 0")
            elif param['blbc']=='dz': problem.add_bc("left(dzb)  = 0")
            elif param['blbc']=='dq': problem.add_bc("left(kappab*dzb+bVZ*s= 0")

            
            


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
ju.check_error()
ju.logger.info('Solver built for initial value problem with integration scheme {}'.format(scheme)) 
    
# Initial conditions
u=None
v=None
w=None
b=None
s=None
p=None
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
nf1=None
nf2=None

if 'u'       in problem.variables:  u     = solver.state['u']
if 'v'       in problem.variables:  v     = solver.state['v']
if 'w'       in problem.variables:  w     = solver.state['w']
if 's'       in problem.variables:  s     = solver.state['s']
if 'b'       in problem.variables:  b     = solver.state['b']
if 'p'       in problem.variables:  p     = solver.state['p']
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
if 'noisef1' in problem.variables:  nf1   = solver.state['noisef1']
if 'noisef2' in problem.variables:  nf2   = solver.state['noisef2']

#domain.dim
#Ux=domain.bases[0].Integrate(u)
#print('Ux {}'.format(Ux))
#y_basis=domain.bases[0]
#z_basis=domain.bases[0]
#x_basis.Integrate(arg0, out=out)

if restart:
    dt=0
    for a in fnrst:
        (dt,state_ev) = JA.load_state(a,solver,state_ev,rparam)
        if dt>0: break
        ju.logger.info('load_state failed with: {}'.format(a))
        dt = ju.restartf(a, solver,-1)
        if dt>0: break
        ju.logger.info('restartf failed with: {}'.format(a))
    if dt==0:
        ju.logger.info('Failed to load any initial state: {}'.format(fnrst))
        ju.set_error('Failed to load initial state')
    ju.check_error()
    if resetdb:
        state_ev['bp'] = ju.comm.bcast(sbp,0)
        ju.logger.info('Resetting density profile')
 
else:
    ju.logger.info('Initialising fields')

    # Get pointers to the simulation variables and set them to zero

    if b  is not None: b['g'] = 0
    if s  is not None: s['g'] = 0
    if u  is not None: u['g'] = 0
    if v  is not None: v['g'] = 0
    if w  is not None: w['g'] = 0
    if p  is not None: p['g'] = 0
    
    if b is not None and 'fb' in problem.parameters and 'wb' in problem.parameters:
        b.set_scales(AA)
        b['g']   = problem.parameters['fb']['g']*problem.parameters['wb']['g']

    if s is not None and  'fs' in problem.parameters and 'ws' in problem.parameters:
        s.set_scales(AA)
        s['g']   = problem.parameters['fs']['g']*problem.parameters['ws']['g']

    if bz is not None: b.differentiate('z', out=bz)
    if sz is not None: s.differentiate('z', out=sz)
    if uz is not None: u.differentiate('z', out=uz)
    if vz is not None: v.differentiate('z', out=vz)
    if wz is not None: w.differentiate('z', out=wz)

if ifb is not None:
    ifb.set_scales(AA)
    if b is not None:
        b.set_scales(AA)
        b['g'] = ifb['g']
        if bz is not None: b.differentiate('z', out=bz)
    if s is not None:
        s.set_scales('AA')
        s['g'] = ifb['g']
        if sz is not None: s.differentiate('z', out=sz)
            
if ifu is not None:
    u['g']   = U*ifu
    if uz is not None: u.differentiate('z', out=uz)
if ifv is not None:
    v['g']   = U*ifv
    if vz is not None: v.differentiate('z', out=vz)
if ifw is not None:
    w['g']   = U*ifw
    if wz is not None: w.differentiate('z', out=wz)
    
if nf2 is not None: nf2['g'] = anoise(1)



#if ju.mpirank==0: print('bz {}'.format(solver.state['bz']['g'][0,0:5]))
#quit(1)

  
if bnoise is not None:
    b.set_scales(1)
    if bntype=='gaussian':  b['g']=b['g']*(1+bnoise*anoise(1))
    if bntype=='uniform':   b['g']=b['g']*(1+unoise(-bnoise,bnoise,1))
    ju.logger.info('Adding noise to bouyancy: type: {:}, amplitude {:8.3f}'.format(bntype,bnoise))
    if bz is not None: b.differentiate('z', out=bz)
    b.set_scales(AA)

if sinitial is not None: s['g']=sinitial
if snoise is not None:
    s.set_scales(1)
    if forcing==0:
        if bntype=='gaussian':  s['g']=sinitial+snoise*anoise(1)
        if bntype=='uniform':   s['g']=sinitial+unoise(-snoise,snoise,1)
    else:
        if bntype=='gaussian':  s['g']=s['g']*(1+snoise*anoise(1))
        if bntype=='uniform':   s['g']=s['g']*(1+unoise(-snoise,snoise,1))
    if sz is not None: s.differentiate('z', out=sz)        
    ju.logger.info('Adding noise to sediment: type: {:}, amplitude {:8.3f}'.format(bntype,snoise))
    s.set_scales(AA)
    ju.logger.info('Sediment range: type: [{:5.3f} {:5.3f}]'.format(s['g'].min(),s['g'].max()))
    
if inoise is not None:
    ju.logger.info('Adding noise to velocities: with amplitude {:8.3f}'.format(inoise))
    ju.add_zm_noise(u,inoise,'u',W*H*L,rand)
    ju.add_zm_noise(v,inoise,'v',W*H*L,rand)
    ju.add_zm_noise(w,inoise,'w',W*H*L,rand)
    if uz is not None: u.differentiate('z', out=uz)
    if vz is not None: v.differentiate('z', out=vz)
    if wz is not None: w.differentiate('z', out=wz)

maxu=0
maxv=0
maxw=0
if u is not None: maxu=ju.jmax(np.abs(u['g']))
if v is not None: maxv=ju.jmax(np.abs(v['g']))
if w is not None: maxw=ju.jmax(np.abs(w['g']))

ju.logger.info('Maximum velocities: [{:6.3f},{:6.3f},{:6.3f}]'.format(maxu,maxv,maxw))

if dt is None:
    if maxu==0 and maxv==0 and maxw==0: dt =np.sqrt(mindt*maxdt)
    else:                               dt = 0.5/max(maxu/dx,maxv/dy,maxw/dz)
dt=max(mindt,min(maxdt,dt))
ju.logger.info('Initial timestep:   [{:6.4f},{:6.4f},{:6.4f}]'.format(mindt,dt,maxdt))

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
    ju.logger.info('Reset complete')

# This needs to be after the state has been loaded and possibly reset so that solver.sim_time is correct    
JA.set_solver(solver,problem,param,reset,state_ev,avar)
if save_initial: JA.save_state(dt,state_ev,fn='initial.hdf5')
if debug_flow: ju.logger.info('ded_gc: ja.jit starting') 
JY=ja.jit(JA,JA.javrg, False, False,  True, False, 0,   None,   param)
if debug_flow: ju.logger.info('ded_gc: ja.jit finished') 
if dtjdivyz is not None:
    dtjdivyzn=np.int64(np.ceil(solver.sim_time/dtjdivyz))
    ju.logger.info('jrec: nm: dtjdivyz dt: {:5.3f} {:5d}'.format(dtjdivyz,dtjdivyzn))

if debug_flow: ju.logger.info('ded_gc: saved xz') 
#if bminmax: JA.savexz('bmax.hdf5',  'bmaxf',bmaxf, AA)
#if sminmax: JA.savexz('smaxf.hdf5', 'smaxf',smaxf, AA)
JA.saved('fxz.hdf5',('bmaxf','smaxf'),'xz',locals(), AA)


if debug_flow: ju.logger.info('ded_gc: ju.write_hdf5_vars') 
ju.write_hdf5_vars('force/fx.hdf5',('pmif','rfur','rful','Bwxx','wmult','bmultl','bmultr','wux1','wux2','wux3','wux4','wux5','wux6','wux7','wux8','Iux1','Iux2','Iux3','Iux4','Iux5','sminf','bminf'),locals())

#ju.logger.info('MM {}, M {}, xx[MM] {}, hb {}'.format(MM,M,xx[MM],hb))
#ju.logger.info('wux3[0:MM] {}'.format(wux3[0:MM]))
#quit(1)

if debug_flow: ju.logger.info('ded_gc: saved yz') 
JA.saved('fyz.hdf5',('inwdr','nwdr','onewdr','rr','wrr','nyz'),'yz',locals(),AA)

   
del(ifb,ifs,ifu,ifv,ifw)

if dim==3:
    if 'oddcx' in problem.parameters: JA.saven('force/oddcx.hdf5','oddcx',problem.parameters['oddcx']['g'][:,0,0].flatten(),AA,'x')
    if 'oddcy' in problem.parameters: JA.saven('force/oddcy.hdf5','oddcy',problem.parameters['oddcy']['g'][0,:,0].flatten(),AA,'y')
    if 'oddcz' in problem.parameters: JA.saven('force/oddcz.hdf5','oddcz',problem.parameters['oddcz']['g'][0,0,:].flatten(),AA,'z')
else:
    if 'oddcx' in problem.parameters: JA.saven('force/oddcx.hdf5','oddcx',problem.parameters['oddcx']['g'][:,0].flatten(),AA,'x')
    if 'oddcz' in problem.parameters: JA.saven('force/oddcz.hdf5','oddcz',problem.parameters['oddcz']['g'][0,:].flatten(),AA,'z')

#JA.savexyz('oddcx.hdf5','oddcx',problem.parameters['oddcx']['g'],1)
 

# Integration parameters
solver.stop_sim_time  = T
vv=()
if u is not None: vv+=('u',)
if v is not None: vv+=('v',)
if w is not None: vv+=('w',)

vs=''
if u is not None: vs+=' u*u'
if v is not None: vs+='+v*v'
if w is not None: vs+='+w*w'

if maxdt>mindt:
    CFL=flow_tools.CFL(solver, initial_dt=dt, cadence=5, safety=0.8, max_change=1.5, min_change=0.5, max_dt=maxdt, min_dt=mindt)
    CFL.add_velocities(vv)


solver.stop_wall_time = ju.get_sec(param['time'])
solver.stop_iteration = np.inf

if debug_flow: ju.logger.info('ded_gc: analysis_tasks')
analysis_tasks=da.set_analysis_tasks(ddir,solver,problem,max_writes,param,mode)
ju.logger.info('Analysis tasks initialised')

#flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
#flow.add_property(vs, name='SS')
#if wq['WB'] is not None: flow.add_property("(b-fb)**2*abs(wb)", name='FB')
#if wq['WS'] is not None: flow.add_property("(s-fs)**2*abs(ws)", name='FS')
#if wq['WU'] is not None: flow.add_property("(u-fu)**2*abs(wu)", name='FU')
#if wq['WV'] is not None: flow.add_property("(v-fv)**2*abs(wv)", name='FV')
#if wq['WW'] is not None: flow.add_property("(w-fw)**2*abs(ww)", name='FW')
#if noise>0: flow.add_property("wnoise*(noiseu*noiseu + noisev*noisev + noisew*noisew)",   name='NS')

ju.check_abort()

smem=ju.MemStats()
#ju.logger.info(smem)

sim_iter_start=solver.iteration
sim_time_start=solver.sim_time
st0=sim_time_start
if PIDT is not None: nt0=np.floor(st0/PIDT)
else:                nt0=0
    

tp= np.float(time.time())


if PIDST is not None:
    if np.int(np.floor(st0/PIDST) % 2)==0: PIDX=PIDS1
    else: PIDX=PIDS2
    ju.logger.info("Stepping PID input to calibrate controller")
    ju.logger.info("PIDST:   {:7.3f},   PIDS1: {:7.4f},   PIDS2: {:7.4f}".format(PIDST, PIDS1, PIDS2))


# We need to make sure that it is fine if the front position is initially a NaN    
tX0=None

def calcX(sType,param,b,s):
    X=None
    if sType=='gc':
        if       S is not None and     B is not None: X=ju.calc_front(b+s)
        elif not S is not None and     B is not None: X=ju.calc_front(b)
        elif     S is not None and not B is not None: X=ju.calc_front(s)
        else: X=0
        if not np.isfinite(X): X=0
        #ju.logger.info("Front position X: {:7.4f}".format(X))
    return(X)

X=calcX(sType,param,b,s)
X0=X
tX0=solver.sim_time
tX1=tX0


if PIDX is not None:
    if ju.mpirank==0:
        if 'PIDIT' in param:
            state_ev['PIDIT']=param.pop('PIDIT',0.0)
            ju.logger.info('Setting PIDIT to {:7.3f}'.format(state_ev['PIDIT']))
        if 'PIDDD' in param:
            state_ev['PIDDD']=param.pop('PIDDD',0.0)
            ju.logger.info('Setting PIDDD to {:7.3f}'.format(state_ev['PIDDD']))

        if 'g' in state_ev and 'g' in param:
            state_ev['g']=param.pop('g')
            ju.logger.info('Setting g to {:7.3f}'.format(state_ev['g']))
        if 'U' in state_ev and 'U' in param:
            state_ev['U']=param.pop('U')
            ju.logger.info('Setting U to {:7.3f}'.format(state_ev['U']))
        
        if state_ev['PIDIT']==0:
            if PIDG: state_ev['PIDIT']=state_ev['g']
            else:    state_ev['PIDIT']=state_ev['U']
        if PIDG and state_ev['g']>0: g=state_ev['g']
        elif state_ev['U']>0:        U=state_ev['U']
    state_ev=ju.comm.bcast(state_ev,0)

    if PIDG: ju.logger.info("PID controller for front position is active on g: [{:7.3f},{:8.4f},{:7.3f}]".format(param['gmin'],state_ev['g'],param['gmax']))
    else:    ju.logger.info("PID controller for front position is active on U: [{:7.3f},{:8.4f},{:7.3f}]".format(param['Umin'],state_ev['U'],param['Umax']))
    ju.logger.info("PIDX: {:7.4f},  PIDIT: {:7.4f}, PIDDD: {:7.1e}".format(PIDX, state_ev['PIDIT'], state_ev['PIDDD']))
    ju.logger.info("PIDT: {:7.4f},   PIDP: {:7.4f}   PIDI: {:7.5f},  PIDD: {:7.4f}".format(PIDT, PIDP, PIDI, PIDD))
    if ju.mpirank==0:
        if PIDG: pid = jpid.pid(PIDP, PIDI, PIDD, PIDX, state_ev['PIDIT'], state_ev['PIDDD'], X,  state_ev['g'], st0, param['gmin'],param['gmax'],5,-1)
        else:    pid = jpid.pid(PIDP, PIDI, PIDD, PIDX, state_ev['PIDIT'], state_ev['PIDDD'], X,  state_ev['U'], st0, param['Umin'],param['Umax'],0.5,1)
     
UF=None
cRe=None

dforce=True

del(wnoise,nwd,nwb,nws,nwu,nwv,nww)
del(nfb,nfd,nfu,nfv,nfw)

if False:
    maxu=ju.jmax(np.abs(u['g']))
    if Ny>1: maxv=ju.jmax(np.abs(v['g']))
    else:    maxv=0
    maxw=ju.jmax(np.abs(w['g']))
    ju.logger.info('Maximum velocities: [{:6.3f},{:6.3f},{:6.3f}]'.format(maxu,maxv,maxw))

check_finite=True
fflag=False
st0=solver.sim_time
si0=solver.iteration
ju.logger.info('Starting main loop');

#ju.print_item_sizes(locals().items())
#ju.print_item_sizes(globals().items())
tsave_state=0
twrite_stats=0
tcalcX=0
tjrec=0


tstep=0
tinitial = np.float(time.time()-tinitial)
tloop=time.time()
k=0
tcheck=tloop

def display_top(snapshot, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

if fpidu is not None: fpidu.set_scales(AA) 

if ddiv:   fd = problem.parameters['fd']
if db and 'fb' in problem.parameters: fb = problem.parameters['fb'] 

if 'pmss' in param:
    if 'fu' in problem.parameters: fu = problem.parameters['fu']  
    if 'fd' in problem.parameters: fd = problem.parameters['fd']  
    if B is not None: fb = problem.parameters['fb']  
    #if S is not None: fs = problem.parameters['fb']  

 
if debug_pympler: mem_tracker = tracker.SummaryTracker()
if debug_tracemalloc:
    snapshot1 = tracemalloc.take_snapshot()
    top_stats = snapshot1.statistics('lineno')
    #print("[ Top 10 ]")
    #for stat in top_stats[:10]: print(stat)
    display_top(snapshot1)

t1= np.float(time.time())
timea.append(t1)
ju.logger.info('Initialisation time: {:8.2f}'.format(timea[1]-timea[0]))
ju.logger.info('Simulation time: {:8.3f}'.format(sim_time_start))
ju.logger.info('Wall time: {}'.format(ju.sec2dhms(time.time())))

if 'divza' in state_ev:          state_ev['divz'][:,0] = state_ev.pop('divza')
if 'divzb' in state_ev and ddiv: state_ev['divz'][:,1] = state_ev.pop('divzb')
#JA.save_state(dt,state_ev,'state3.hdf5')

if sType=='gc' and not ddiv: fd = None


if rescale:
    if dim==2: vsqr = de.operators.integrate(u*u+w*w,     'x',      'z')
    else:      vsqr = de.operators.integrate(u*u+v*v+w*w, 'x', 'y', 'z')
while solver.ok and not ju.abort:
    #    ju.MemTrack()
    #    if debug_mem: ju.MemPrint()
    gc.collect()
    if debug_gc:
        for item in gc.garbage:
            print(item)
            assert list_id == id(item)

    
    #if debug_mem: print(psutil.Process(ju.pid).memory_info())
    if debug_pympler: mem_tracker.print_diff()
    if debug_tracemalloc:
        snapshot2 = snapshot1
        snapshot1 = tracemalloc.take_snapshot()
        top_stats = snapshot1.compare_to(snapshot2, 'lineno')
        for stat in top_stats[:1]: print(stat)
        #print(top_stats[0].traceback.format())
    k=k+1
    #if k==1 or k==10: ju.MemAll()
    if k<3:  # Times for calulating efficiency
        sti=solver.sim_time # Initial simulation time
        wti=np.float(time.time())      # Initial wall time

    if rescale:
        if dim==2: sf = np.sqrt((vsqr.evaluate())['g'][0,0]/(H*L))
        else:      sf = np.sqrt((vsqr.evaluate())['g'][0,0,0]/(H*W*L))
        #ju.logger.info('rescale: sf={:7.5f}'.format(sf))
        if u  is not None: u['g']  = u['g']/sf
        if v  is not None: v['g']  = v['g']/sf
        if w  is not None: w['g']  = w['g']/sf
        if uz is not None: uz['g'] = uz['g']/sf
        if vz is not None: vz['g'] = vz['g']/sf
        if wz is not None: wz['g'] = wz['g']/sf
        if p  is not None: p['g']  = p['g']/sf**2
        
    if check_finite:
        if b is not None: fflag = fflag or not np.all(np.isfinite(b['g']))
        if s is not None: fflag = fflag or not np.all(np.isfinite(s['g']))
        if u is not None: fflag = fflag or not np.all(np.isfinite(u['g']))
        if v is not None: fflag = fflag or not np.all(np.isfinite(v['g']))
        if w is not None: fflag = fflag or not np.all(np.isfinite(w['g']))
        if fflag:
            ju.set_error('Fields not finite')
            ju.logger.info('Fields not finite')
        ju.check_error()
        
    if forcing==7 and sType=='gc':
        if db:
            b.set_scales(AA)
            if Ny > 1:
                bp=JY.int(b['g'])/W
                #print('yycolor {} b.shape {} bp.shape {}'.format(yycolor,b['g'].shape,bp.shape))
                #if ju.mpirank==0: bp = dxx*np.sum(wux2.reshape((Nxx,1))*bp,axis=(0,))
                if yycolor: bp = dxx*np.sum(wux2.reshape((Nxx,1))*bp,axis=(0,))
                #print('yycolor {} b.shape {} bp.shape {}'.format(yycolor,b['g'].shape,bp.shape))
                #quit(1)
                #bp=ju.comm.bcast(bp,0)
            else: bp = dxx*np.sum(wux2*b['g'],axis=(0,))
            #print('rank {} yyrank {} bp {} zzcounts {}'.format(ju.mpirank,yycomm.Get_rank(),bp,zzcounts))
            if yycolor:

                #print('bp {}'.format(bp))
                yycomm.Gatherv(sendbuf=bp, recvbuf=(gbp, zzcounts), root=0)
                #print('gbp {}'.format(gbp))
                
            if ju.mpirank==0:
                avar['bp']=gbp
                Igbp=(gIzz*gbp).sum()/hb
                if Igbp!=0: gbp=gbp/Igbp
                state_ev['bp']=(param['dT']*state_ev['bp']+dt*gbp)/(param['dT']+dt)
                state_ev['bp']=state_ev['bp']-state_ev['bp'].min()
                Ibp=(gIzz*state_ev['bp']).sum()/hb
                if Ibp!=0: state_ev['bp']=ju.clip01(state_ev['bp']/Ibp).reshape((Nzz,))
                if not np.all(np.isfinite(state_ev['bp'])):
                    ju.logger.info('bp not finite')
                    ju.logger.info('Ibp/hb {:6.3f}'.format((gIzz*state_ev['bp']).sum()/hb))
                    ju.srange('gbp range [{:6.3f},{:6.3f}]',gbp)
                    ju.srange('bp  range [{:6.3f},{:6.3f}]',state_ev['bp'])
                    ju.set_error('bp not finite')
            ju.check_error()
            state_ev['bp'] = ju.comm.bcast(state_ev['bp'],0)
            if fbmult:
                bps=(state_ev['bp'][zzslices]).reshape(szzzl)
                b['g'][bmultls,...] += (bps - b['g'][bmultls,...])*bmultl*bps
            else:
                if fbmax: bmaxf = ((state_ev['bp'][zzslices]).reshape(szzzl))*wux3
                else:   fb['g'] = ((state_ev['bp'][zzslices]).reshape(szzzl))
                
            if np.mod(solver.iteration,5)==0 and False:
                problem.parameters['fb'].set_scales(AA)
                problem.parameters['wb'].set_scales(AA)
                #JA.savexz('fb/fb-{:05d}.hdf5'.format(np.int(solver.iteration/5)),('fb','wb'),(problem.parameters['fb']['g'],problem.parameters['wb']['g']),AA)
                problem.parameters['fb'].set_scales(1)
                problem.parameters['wb'].set_scales(1)
                JA.savexz('fb1/fb-{:05d}.hdf5'.format(np.int(solver.iteration/5)),('fb','wb'),(problem.parameters['fb']['g'],problem.parameters['wb']['g']),1)
                ju.logger.info('fb parity {}'.format(problem.parameters['fb'].meta['x']['parity']))
                ju.logger.info('wb parity {}'.format(problem.parameters['wb'].meta['x']['parity']))
                ju.logger.info('bottom right FB {:6.1f} fb {:6.4f} wb {:6.4f} '.format(problem.parameters['Fb'],problem.parameters['fb']['g'][-1,0],problem.parameters['wb']['g'][-1,0]))
        else:
            if fbmult:
                b.set_scales(AA)
                b['g'][bmultls,...] += (hbzz - b['g'][bmultls,...])*bmultl*hbzz
                 
        if ddiv:
            w.set_scales(AA)
            WA=JY.int(w.differentiate('z')['g'])
            if yycolor:
                dw  = -dxx*np.sum(wux1.reshape((Nxxl,1))*WA.reshape((Nxxl,Nzzl)),axis=(0,))
                yycomm.Gatherv(sendbuf=dw, recvbuf=(gdw, zzcounts), root=0)
            if ju.mpirank==0:
                if Iwux1 !=0 : gdw=gdw/Iwux1 # Different normalisation
                avar['dwdz'] = -gdw/W
                state_ev['divz'] += (dt/param['dT'])*gdw
                Idivz=(gIzz*state_ev['divz']).sum()
                if Idivz!=0: state_ev['divz'] =-state_ev['divz']*H*U/Idivz
                state_ev['divz'] = state_ev['divz'].reshape((Nzz,))
                Igdw   = np.sqrt((gIzz*(gdw                   )**2).sum(axis=(0,))/(H*W))
                Idivzs = np.sqrt((gIzz*(state_ev['divz']-divz0)**2).sum(axis=(0,))/(H*W))
                Idivz   =        (gIzz* state_ev['divz']          ).sum(axis=(0,))/(H*W)
                Idivz0  =        (gIzz* divz0                     ).sum(axis=(0,))/(H*W)
                if not np.all(np.isfinite(state_ev['divz'])):
                    ju.logger.info('divz is not finite')
                    ju.logger.info('Idivzs   {}'.format(Idivzs))
                    ju.logger.info('Idivz0   {}'.format(Idivz0))
                    ju.logger.info('Idivzs0  {}'.format((gIzz*state_ev['divz'][:,0]).sum()/(H*W)))
                    ju.srange('gdw range [{:6.3f},{:6.3f}]',gdw)
                    ju.srange('divz  range [{:6.3f},{:6.3f}]',state_ev['divz'])
                    ju.set_error('bp not finite')
            ju.check_error()
            state_ev['divz'] = ju.comm.bcast(state_ev['divz'],0)
            fd['g']            = U*wux5 +            wux1.reshape(szxx)*state_ev['divz'][zzslices].reshape(szzzl) #missing U???   
            if   gclf and gcrf: rfu['g'] = U*(rfur + rful.reshape(szxx)*state_ev['divz'][zzslices].reshape(szzzl)**2)
            elif gclf:          rfu['g'] = U*(       rful.reshape(szxx)*state_ev['divz'][zzslices].reshape(szzzl)**2)
            elif gcrf:          rfu['g'] = U*(rfur)
            #ju.logger.info('fd    {:8.2e}'.format(dVV*ju.jsum(Iww[-1]*fd['g'])))
            #ju.logger.info('Iwux5 {:8.2e}'.format(dxx*wux5.sum()-1))
            #ju.logger.info('Iwux1 {:8.2e}'.format(dxx*wux1.sum()-1))

                
    if fbmax or fbmin or bmultl is not None or bmultr is not None: b.set_scales(AA)
    if fbmax:  b['g'] = np.maximum(bmaxf,b['g'])
    if fbmin:  b['g'] = np.minimum(bminf,b['g'])
 
    if fsmax or fsmin or smultl is not None or smultr is not None: s.set_scales(AA)
    if fsmax:  s['g'] = np.maximum(smaxf,s['g'])
    if fsmin:  s['g'] = np.minimum(sminf,s['g'])

    if bmultr is not None: b['g'][bmultrs,...] -= b['g'][bmultrs,...] * bmultr
    if smultr is not None: s['g'][smultrs,...] -= s['g'][smultrs,...] * smultr

    if 'pmss' in param: # moving source for the plume
        if ju.mpirank==0:
            state_ev['SY']=ju.noiseint(state_ev['SY'],dt,param['pmsl'],param['pmss'])
            state_ev['SZ']=ju.noiseint(state_ev['SZ'],dt,param['pmsl'],param['pmss'])
            RR=ju.get_radius(R,state_ev['SY'],state_ev['SZ'],Py,Pz)
        state_ev['SY'] = ju.comm.bcast(state_ev['SY'],0)
        state_ev['SZ'] = ju.comm.bcast(state_ev['SZ'],0)
        RR             = ju.comm.bcast(RR,0)
        avar['SY']=state_ev['SY']
        avar['SZ']=state_ev['SZ']
        avar['R']=RR
        wrr=InletA*ju.pinlet(yy-state_ev['SY'],zz-state_ev['SZ'],inlet,RR,maxkr,nfr)[0]
        if forcing==7:    fd['g'] = U*wrr*wux1
        if B is not None: fb['g'] = wrr*Bwxx
        #if S is not None: fs = fb

        if not ddiv and not forcing==7: fu['g'] = wrr*wux4+mu*(1-wux4)*inwdr 
        #ju.logger.info('Inlet area {:7.4f},  R={:7.4f}, Y={:7.4f}, Z={:7.4f}'.format(dyy*dzz*ju.jsum(wrr)/IA,RR,state_ev['SY'],state_ev['SZ']))
        #ju.srange('fb range [{:6.3f},{:6.3f}]',fb)
    if dxnoise is not None: dxnoise.int(dt)
    if dynoise is not None: dynoise.int(dt)
    if dznoise is not None: dznoise.int(dt)
    
    if sType=='pm' and forcing==7:
        u.set_scales(AA)
        fd.set_scales(AA)
        if ddiv:
            # find u velocity and gradient near the edges
            Iu    = ju.comm.reduce(np.sum(nwdr*u['g'],axis=(1,2)), op=MPI.SUM)
            Idudx = ju.comm.reduce(np.sum(nwdr*u.differentiate('x')['g'],axis=(1,2)), op=MPI.SUM)
            if ju.mpirank==0:
                Idudx = dAA*Idudx.reshape(szxx)
                avar['u']  = dAA*Iu.reshape((Nxx,))
                avar['Iu'] = dxx*avar['u'].sum()/L
                if pmIu: state_ev['divx']  = np.maximum(0,wux6*(state_ev['divx'].reshape(szxx)-wux4*(dt/param['dT']*Idudx+Anwdr*avar['Iu']*10*dt)))
                else: state_ev['divx']  = np.maximum(0,wux6*(state_ev['divx'].reshape(szxx)-dt/param['dT']*Idudx))
                #state_ev['divx'] = np.maximum(0,wux6*(state_ev['divx'].reshape(szxx)-dt/param['dT']*wux4*Idudx))
                
                
                if fixI: state_ev['divx'] = pmif*(state_ev['divx']*pmif).sum()*dxx
                avar['dudx']=Idudx.reshape((Nxx,))
            state_ev['divx'] = ju.comm.bcast(state_ev['divx'],0)
            fd2 = state_ev['divx']*nwdr
            Ifd2    = dVV*ju.jsum(fd2)
            state_ev['divx'] = state_ev['divx'].reshape((Nxx,))
        else:
            fd2=0
            Ifd2=0
        if topd: # Dynamic evolution of top divergence
            u.set_scales(AA)
            v.set_scales(AA)
            w.set_scales(AA)
            if topdfield:
                avar['topu']['g'] = dxx*np.sum(wux3*u['g'],axis=(0,)) # Top u
                avar['topv']['g'] = dxx*np.sum(wux3*v['g'],axis=(0,)) # Top v
                avar['topw']['g'] = dxx*np.sum(wux3*w['g'],axis=(0,)) # Top w
                avar['topd'] = avar['topv'].differentiate('y')+avar['topw'].differentiate('z')
                tops = np.sqrt(dAA*(ju.jsum(avar['topu']['g']**2+avar['topv']['g']**2)))
                rmstopd = np.sqrt(dAA*(ju.jsum(avar['topd']['g']**2)))
                fd3r = np.maximum(0, (param['dT']*state_ev['topd']['g']+avar['topd']['g']*dt)/(dt+param['dT']))
                ju.logger.info('type(fd3r) {}'.format(type(fd3r)))
                ju.logger.info('type(topd) {}'.format(type(state_ev['topd'])))
                ju.logger.info('type(topd[g]) {}'.format(type(state_ev['topd']['g'])))

                Ifd3r = dAA*ju.jsum(fd3r)
                if Ifd3r<0:
                    fd3r = np.ones((Nyyl,Nzzl))/(H*W)
                    state_ev['topd']['g']=fd3r
                else:
                    fd3r = fd3r/Ifd3r
                    state_ev['topd']['g']=fd3r/Ifd3r
            else:
                lu = dxx*np.sum(wux3*u['g'],axis=(0,)) # Top u
                lv = dxx*np.sum(wux3*v['g'],axis=(0,)) # Top v
                lw = dxx*np.sum(wux3*w['g'],axis=(0,)) # Top w
                ju.gatheryz(topu,lu,gyyslices,gzzslices)
                ju.gatheryz(topv,lv,gyyslices,gzzslices)
                ju.gatheryz(topw,lw,gyyslices,gzzslices)
                if ju.mpirank==0:
                    gtopd=ju.fft_diff(topv,dyy,Py,0)+ju.fft_diff(topw,dzz,Pz,1)
                    avar['topd']=gtopd
                    avar['topu']=topu
                    avar['topv']=topv
                    avar['topw']=topw
                    tops   = np.sqrt(dAA*((topv**2).sum()+(topw**2).sum()))
                    rmstopd = np.sqrt(dAA*((gtopd**2).sum()                ))
                    state_ev['topd'] = np.maximum(0, (param['dT']*state_ev['topd']+gtopd*dt)/(dt+param['dT']))
                state_ev['topd']=ju.comm.bcast(state_ev['topd'],0)
                fd3r=state_ev['topd'][yyslices,zzslices]
                Ifd3r = dAA*ju.jsum(fd3r)
                if Ifd3r<0: fd3r = np.ones((Nyyl,Nzzl))/(H*W)
                else:
                    fd3r = fd3r/Ifd3r
                    state_ev['topd']=state_ev['topd']/Ifd3r
        else: 
            fd3r  = dxx*np.sum(wux3*np.maximum(IAU/(H*W),u['g']),axis=(0,))
            Ifd3r = dAA*ju.jsum(fd3r)
            fd3r  = fd3r/Ifd3r
        if False:
            r3=R
            Itu   = ju.comm.reduce(np.sum(np.sum(wux3*u['g']   ,axis=(0,))), op=MPI.SUM)  
            Ituu  = ju.comm.reduce(np.sum(np.sum(wux3*u['g']**2,axis=(0,))), op=MPI.SUM)  
            if ju.mpirank==0:
                Itu   = dVV*Itu/math.Pi   # based on top-hat or \int f^2 = int f profile
                Ituu  = dVV*Ituu/math.Pi
                r3=max(R,Itu/np.sqrt(Ituu))
                avar['rtop']=r3
            r3=ju.comm.bcast(r3,0)
            fd3r = jtheta.periodic_gaussianu(np.sqrt(2*math.pi)/r3,yy,zz,FW,FH)/(W*H)
            #Ifd3ru  = dAA*ju.comm.allreduce(np.sum(fd3r),    op=MPI.SUM)
            #Ifd3ruu = dAA*ju.comm.allreduce(np.sum(fd3r**2), op=MPI.SUM)
            #ju.logger.info('r3: {:8.4f} Ifd3ru/sqrt(Ifd3ruu): {:8.4f} Ifd3ru: {:8.4f} Ifd3ruu {:8.4f} '.format(r3,Ifd3ru/np.sqrt(Ifd3ruu),Ifd3ru,Ifd3ruu))

        if dtjdivyz is not None:
            if solver.sim_time >= dtjdivyz*dtjdivyzn:
                dtjdivyzn +=1
                JA.saveyz(ddir / 'divyz/divyz-{:05d}.hdf5'.format(dtjdivyzn), 'divyz', fd3r, AA)
        fd1  = U*wrr*wux1
        Ifd1 = dVV*ju.jsum(fd1)
        fd3     = -(Ifd1+Ifd2)*wux2*fd3r.reshape((1,Nyyl,Nzzl))
        if 'topd' in state_ev:
            if topd:  state_ev['topd']['g'] = (Ifd1+Ifd2)*state_ev['topd']['g']
            else: state_ev['topd'] = (Ifd1+Ifd2)*state_ev['topd']
        fd['g'] = fd1+fd2+fd3
        Ifd3    = dVV*ju.jsum(fd3)
        Ifd     = dVV*ju.jsum(fd['g']) 
        #ju.logger.info('Ifd1: {:8.4f}, Ifd2: {:8.4f}, Ifd3: {:8.4f}, Ifd: {:9.2e}'.format(Ifd1/IAU,Ifd2/IAU,Ifd3/IAU,Ifd/IAU))
        #if np.mod(solver.iteration,10)==0: 

        #ju.logger.info('r3: {:8.4f}, IIu: {:8.4f}, Itu: {:8.4f}, Ituu: {:8.4f}, Ifd1: {:8.4f}, Ifd2: {:8.4f}, Ifd3: {:8.4f}, Ifd: {:9.2e}'.format(r3,IIu,Itu/IAU,Ituu/IAUU,Ifd1/IAU,Ifd2/IAU,Ifd3/IAU,Ifd/IAU))
              
    if ddiv and sType=='pm' and forcing==6:
        u.set_scales(AA)
        fu.set_scales(AA)
        rfu.set_scales(AA)
        fd.set_scales(AA)

        uE=dVV*ju.jsum(Iux3*(u['g']-U*wrr))     # Find where velocity should be Gaussian and constant
        Isb=g*(ju.jint(s)+ju.jint(b))  # Find total weight of bouyancy
        #ju.logger.info('Is {} Ib {} Isb {} Isu {} g {} param[g] {} dt {}'.format(ju.jint(s),ju.jint(b),Isb,ju.jint(u),g,param['g'],dt))

        state_ev['Fx']=-Isb    #state_ev['Fx']=-Isb-uE*H*W*L/dt  
        rfu['g']=state_ev['Fx']*Iux2           # Add the forcing over the entire velocity forced region
        problem.parameters['Fx']=state_ev['Fx']

        #mu  = dVV*ju.jsum(Iux0*u['g'])                                                         # Find mean incoming velocity per unit area
        #u.differentiate('x',out=dudx)
        #qv1 = dAA*np.reshape(wux5,(Nxx,))*ju.comm.allreduce(-np.sum(nwdr*dudx['g'],axis=(1,2)), op=MPI.SUM)        # Find downwelling profile
        #qv1 = dAA*np.reshape(wux5,(Nxx,))*ju.comm.allreduce(-np.sum(nwdr*u['g'],axis=(1,2)), op=MPI.SUM)        # Find downwelling profile
        qv1 = dAA*ju.comm.allreduce(np.sum(nwdr*u['g'],axis=(1,2)), op=MPI.SUM)        # Find downwelling profile
        #qv1 = dVV*np.reshape(wux5,(Nxx,))*np.sum(np.reshape(Iux5,(Nxx,))*ju.comm.allreduce(-np.sum(nwdr*u['g'],axis=(1,2)), op=MPI.SUM))
        if ju.mpirank==0:
            #qv1 = ju.poly_proj(qv2,rwux5,4)
            state_ev['qvx'] = state_ev['qvx']-dt/param['dT']*qv1
            state_ev['qvx'] = ju.poly_proj(state_ev['qvx'],rwux5,4)
            state_ev['qvx'] = np.maximum(-2*U*rwux5,np.minimum(2*U*rwux5,state_ev['qvx']))
            #ju.write_hdf5('qv1/qv1-{:05d}.hdf5'.format(solver.iteration),'qv1',qv1)
            #ju.write_hdf5('qvx/qvx-{:05d}.hdf5'.format(solver.iteration),'qvx',state_ev['qvx'])
            #state_ev['qvx']=(2*state_ev['qvx']+np.roll(state_ev['qvx'],1)+np.roll(state_ev['qvx'],-1))/4
            
        state_ev['qvx'] = ju.comm.bcast(state_ev['qvx'],0)

        EQ=np.sum(state_ev['qvx'])*dxx # Extra flux
        mu = (U*IA+EQ)/(H*W)

        if dxnoise is not None: dxnoise.f['g']=dxnoise.f['g']*wnx*wr/wrmax
        if dynoise is not None: dynoise.f['g']=dynoise.f['g']*wnx*wr/wrmax
        if dznoise is not None: dznoise.f['g']=dznoise.f['g']*wnx*wr/wrmax
        if dnoise is not None and False:
            state_ev['noise'].set_scales(1)
            ee=np.exp(-dt/noiseT)  # The length scaling is hard to get correct np.sqrt(np.sqrt(8/math.pi)*noiseL)**3

            #Wnoise=snoise(dnoisef1,dnoisewx,dnoisewy,dnoisewz,1)
            #IWnoise=np.sqrt(ju.jsum(Wnoise**2)/(Nx*Ny*Nz))
            # It will end up depending on grid spacing and possible system size
            # Fix the amplitude every step to remove this problem and give instant convergence
            state_ev['noise']['g'] = ee*state_ev['noise']['g']+np.sqrt(1-ee**2)*ssnoise(dnoisef1,dnoisewx,dnoisewy,dnoisewz,1)
            #Inoise=np.sqrt(ju.jsum(state_ev['noise']['g']**2)/(Nx*Ny*Nz))
            #state_ev['noise']['g'] =state_ev['noise']['g']/Inoise
            Snoise=dnoise*np.sqrt(2)*noiseL
            state_ev['noise'].differentiate('y',out=problem.parameters['noisew'])
            state_ev['noise'].differentiate('z',out=problem.parameters['noisev'])
            NZ=np.sqrt(ju.jsum(problem.parameters['noisew']['g']**2)/(Nxx*Nyy*Nzz))*Snoise
            NY=np.sqrt(ju.jsum(problem.parameters['noisev']['g']**2)/(Nxx*Nyy*Nzz))*Snoise
            problem.parameters['noisew']['g']= problem.parameters['noisew']['g']*wnxx*wrr/wrmax*Snoise
            problem.parameters['noisev']['g']=-problem.parameters['noisev']['g']*wnxx*wrr/wrmax*Snoise
            #JA.savexyz('noise/noise-{:05d}.hdf5'.format(solver.iteration),('noise','noisev','noisew'),(state_ev['noise']['g'],problem.parameters['noisev']['g'],problem.parameters['noisew']['g']),AA)
        #fu['g'] = U*(1+dnoise*state_ev['noise']['g'])*wrr*wux4+mu*(1-wux4)*inwdr      # Impose constant region and then Gaussian region
        fu['g'] = U*wrr*wux4+mu*(1-wux4)*inwdr      # Impose constant region and then Gaussian region
        fu.differentiate('x',out=fd)
        fd['g'] = fd['g'] * wux1  + np.reshape(state_ev['qvx'],(Nxx,1,1))*nwdr
        Iqv=dVV*ju.jsum(fd['g'])  # Now add negative divergence to make the system divergence free
        fd['g']=fd['g']-Iqv*Iux1  # This correction should be small becuase of mu imposition
        
        #ju.logger.info('uE {:9.2e}, dnoise {:8.4f}, NYZ {:8.4f}'.format(uE,dnoise,np.sqrt(NY**2+NZ**2)))
        #ju.logger.info('uE {:9.2e}, EQ {:8.4f}, s+b {:8.4f}, Fx {:8.4f}, Ifd {:9.2e}, NN {:8.4f}'.format(uE,EQ/(U*IA*H*W),Isb,state_ev['Fx'],ju.jint(fd),NN))


        #JA.savexyz('fu/fu-{:05d}.hdf5'.format(solver.iteration),'u',fu['g'],AA)
        #JA.savexyz('fd/fd-{:05d}.hdf5'.format(solver.iteration),'d',fd['g'],AA)
 
                
    if tX1 is not None: tX0=copy.copy(tX1)
    if sType=='gc' and np.isfinite(X): XO=copy.copy(X)

    if maxdt>mindt: dt = CFL.compute_dt()
    else:           dt = maxdt
    zmn.args = [dt]
        
    tt0=time.time()
    JA.update(dt)
    tjrec += np.float(time.time()-tt0)  
    
    tt0=time.time()
    solver.step(dt)
    tstep += np.float(time.time()-tt0)  
    if k==0:ju.write_time(solver.sim_time)
        
    if clipB and b is not None: b['g']=np.maximum(minB,np.minimum(maxB,b['g']))
    if clipS and s is not None: s['g']=np.maximum(minS,np.minimum(maxS,s['g']))
    
    if clipu is not None:
        if u is not None: u['g']=np.maximum(-clipu,np.minimum(clipu,u['g']))
        if v is not None: v['g']=np.maximum(-clipu,np.minimum(clipu,v['g']))
        if w is not None: w['g']=np.maximum(-clipu,np.minimum(clipu,w['g']))

    st=solver.sim_time
    si=solver.iteration
    if PIDT is not None: nt=np.floor(st/PIDT)
    else:                nt=nt0+1
    if sType=='gc':
        tt0=time.time()  
        XX=calcX(sType,param,b,s)
        tcalcX += np.float(time.time()-tt0)  
        if np.isfinite(XX): X=copy.copy(XX)
        if np.isfinite(X) and np.isfinite(X0) and tX1>tX0:
            tX1=st
            UC=(X-XO)/(tX1-tX0)
            UF = UC + np.exp(-dt)*(UF-UC) # Calculate smoothed front velocity
            cRe=RL*np.abs(U+UF)*Re
    else: cRe=RL*U*Re
         
    if cRe is not None:
        if not np.isfinite(cRe):
            ju.logger.error('Reynolds number is infinite or NaN')
            ju.set_error()
    
    ju.check_abort_file()
    ju.test_abort()
    
    t2= np.float(time.time())
 
    if ju.mpirank==0:
        if PIDST is not None: 
            if np.int(np.floor(st/PIDST) % 2)==0: NPIDX=PIDS1
            else: NPIDX=PIDS2
            if NPIDX!=PIDX:
                ju.logger.info('Updating PIDX to {:6.4f} from {:6.4f}'.format(NPIDX,PIDX));
                PIDX=NPIDX
                pid.update_XT(NPIDX)
                
    if PIDX is not None:
        if ju.mpirank==0:
            if PIDG: state_ev['g'],state_ev['PIDIT'],state_ev['PIDDD']=pid.update(X,st)
            else:    state_ev['U'],state_ev['PIDIT'],state_ev['PIDDD']=pid.update(X,st)
        state_ev['PIDIT']=ju.comm.bcast(state_ev['PIDIT'],0)
        state_ev['PIDDD']=ju.comm.bcast(state_ev['PIDDD'],0)
        if PIDG:
            state_ev['g']=ju.comm.bcast(state_ev['g'],0)
            g                       = state_ev['g']
            nfg['g']                = state_ev['g']
            problem.parameters['g'] = nfg
        else:
            state_ev['U']=ju.comm.bcast(state_ev['U'],0)
            U                        = state_ev['U']
            problem.parameters['U']  = state_ev['U']
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

    if dtstats is not None or t2>t1+30:
        tt0=time.time()  
        ss2=ju.write_stats(problem,solver,domain,dt,X,U,g,state_ev,gIxx,gIyy,gIzz,Idivzs,Igdw,tops,rmstopd)
        twrite_stats += np.float(time.time()-tt0)  
        
    if t2>t1+30:
        ju.write_time(st)
        ju.jup()
        smem=ju.MemStats()
        if dtstats is not None: ju.flush_stats()
        mdt=(st-st0)/(si-si0)
        ss1 = 'n:{:6d}, t:{:7.3f}, dt:{:7.5f}, e: {:7.2f}, {} '.format(solver.iteration,st,mdt,24*60*60*(st-sti)/(t2-wti),ju.MemShort())
        st0=copy.copy(st)                                         
        si0=copy.copy(si)                                         
        t1=copy.copy(t2)
        ju.logger.info(ss1+ss2)
    if t2>tcheck+60*dtjcheck and dtjcheck>0:
        tt0=time.time()

        JA.save_state(dt,state_ev)
        tcheck=t2
        tsave_state += np.float(time.time()-tt0)  
    
    #if np.mod(solver.iteration,100)==0: JA.savexyz('fd/fd-{:05d}.hdf5'.format(np.int(solver.iteration/100)),('fd'),(problem.parameters['fd']['g']),AA)


timea.append(np.float(time.time()))

tloop = np.float(time.time()-tloop)
tfinal=time.time()

#ju.MemTrack()
#ju.MemPrint()
#ju.MemAll()
#gc.collect()
#ju.MemAll()

 
if dtstats is not None: ju.close_stats()
        
smem=ju.MemStats()
#ju.logger.info(smem)
                
ju.test_error() 
if not ju.jerror:
    tt0=time.time()  
    ju.logger.info('Updating jrec')
    JA.update(dt)
    JA.save_state(dt,state_ev)
    tjrec += np.float(time.time()-tt0)


ju.write_time(solver.sim_time)
ju.logger.info('Saving jrec')
tt0=time.time()  
JA.save()
JA.end()
tjsave = np.float(time.time()-tt0)  

ju.end_stats(timea,sim_time_start,solver,domain)  # Print statistics about the run     

if not parallel:
    post.merge_process_files(str(ddir)+'/final/', cleanup=True)
    post.merge_process_files(str(ddir)+'/checkpoint/', cleanup=True)
    for task in analysis_tasks:
        post.merge_process_files(task.base_path, cleanup=True)

ju.finish(timea[0],solver)

smem=ju.MemStats()
#ju.logger.info(smem)


if PIDX is not None:
    if PIDG: param['g'] = state_ev['g']
    else:    param['U'] = state_ev['U']
    ju.write_param(param)

tfinal = np.float(time.time()-tfinal)
dti=np.float(time.time()-wti0)      # Total wall time
ju.logger.info('Time fractions')
ju.logger.info('  write_stats: {:5.2f}%'.format(100*twrite_stats/dti))
ju.logger.info('  calcX:       {:5.2f}%'.format(100*tcalcX/dti))
ju.logger.info('  jrec:        {:5.2f}%'.format(100*tjrec/dti))            
ju.logger.info('  jsave:       {:5.2f}%'.format(100*tjsave/dti))            
ju.logger.info('  step:        {:5.2f}%'.format(100*tstep/dti))            
ju.logger.info('  initial:     {:5.2f}%'.format(100*tinitial/dti))            
ju.logger.info('  loop:        {:5.2f}%'.format(100*tloop/dti))            
ju.logger.info('  final:       {:5.2f}%'.format(100*tfinal/dti))            
ju.logger.info('  save_state:  {:5.2f}%'.format(100*tsave_state/dti))            
ju.logger.info('Complete')

if debug_gc: gc.set_debug(0)

ju.exit(0)


