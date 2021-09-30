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
  -L, --Length=<>       : Box length       
  -W, --Width=<>        : Box width        
  -H, --Height=<>       : Box height       
  -R, --Radius=<>       : Inlet radius      
  -U, --Velocity=<>     : Flow velocity    
  -T, --Time=<>         : Simulation time  
  -f, --Forcing=<>      : Forcing type     
  -w, --WallTime=<>     : Max walltime, 0 means infinite 
  -c, --dtcheck=<>      : Checkpoint time minutes       
  -q, --Angle=<>        : Slope angle (degrees)    
  -g, --Gravity=<>      : Gravity
  -B, --Buoyancy=<>     : Buoyancy
  -S, --Sediment=<>     : Sediment
  -F, --Force=<>        : Strength of forcing
  -C, --Conservative=<> : Conservative formulation of advection terms
  -I, --Inlet=<>        : Inlet type gaussian, circle, square, 
  -N, --Noise=<>        : amplitude for noise forcing
  -Q, --Qu=<>           : Current flux 
      --save_initial=<> : Save initial state 
      --meanu=<>        : Impose a mean velocity
      --sType=<>        : Type of simulation
      --preset=<>       : Run a simulation with preset parameters
      --setPID          : Set PID defaults
      --scheme=<>       : Integration scheme
      --Re=<>           : Reynolds number  
      --forceb=<>       : Multiplier to strength of forcing for bouyancy
      --forces=<>       : Multiplier to strength of forcing for sedimentation
      --fuspecial=<>    : Special test case of psi = sqrt(L^2+H^2) U sin(2 Pi x/L) sin(2 Pi z/H)
      --fuspecial1=<>   : parameter 1 for fuspecial
      --fuspecial2=<>   : parameter 2 for fuspecial
      --clipu=<>        : Clip velocities to less than this absolute value
      --clipB=<>        : Clip buoyancy between 0 and B or B and 0
      --clipS=<>        : Clip sediment between 0 and S or S and 0
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
      --Ped=<>          : Peclet number for dynamic divergence 
      --fixI=<>         : Fix fluid injection profile
      --hu=<>           : Depth of velocity forcing region
      --hb=<>           : Depth of buoyancy forcing region
      --hs=<>           : Depth of sediment forcing region
      --wu=<>           : Width of velocity forcing region
      --wb=<>           : Width of buoyancy forcing region
      --ws=<>           : Width of sediment forcing region
      --U1=<>           : Current velocity 
      --Nx=<>           : x grid points    
      --Ny=<>           : y grid points   
      --Nz=<>           : z grid points   
      --Tx=<>           : type of x grid
      --Ty=<>           : type of y grid
      --Tz=<>           : type of z grid
      --oddg=<>         : g is a field with odd parity
      --hab=<>          : Height of region to fill with material
      --xa=<>           : Start of region to fill with material
      --xb=<>           : End of region to fill with material
      --dimple=<>       : Amplitude of dimples to add
      --dimplewy=<>     : y wavenumber decay
      --dimplewz=<>     : z wavenumber decay
      --topdivr=<>      : Calculate radial divergence at the top

#Dynamic divergence inlet condition
      --dd=<>           : Dynamic divergence evolution
      --ddT=<>          : Timescale for adding divergence 
      --ddserf=<>       : Project dynamic velocity onto error function
      --ddj=<>          : Number of derivatives to vanish
      --ddk=<>          : Over sampling ratio
      --ddq=<>          : Dynamic divergence evolution using integrated radial velocity
      --ddws=<>         : Dynamic divergence: sampling width
      --dddd=<>         : Diffuse domain dynamic divergence evolution
      --ddscl=<>        : Dynamic divergence scale
      --dbu=<>          : <False> Set velocity to zero input flux
      --ddSc=<>         : <None>  Schmidt number for dynamic divergence 
      --ddD=<>          : <False> Perform diffusion on dynamic divergence 
      --ddlogh=<>       : <False> Use logistic equation evolution to move towards hu
      --ddlogT=<>       : Time scale for logistic equation can be zero

#Injection properties for plumes
      --divxI=<>        : Scale divx so that u is zero
      --wdivxl=<>       : Scale factor for left wdivx
      --wdivxr=<>       : Scale factor for right wdivx

#Dynamic buoyancy divergence inlet condition
      --db=<>           : dynamic buoyancy
      --dbserf=<>       : dynamic buoyancy: Project onto symmetric error function
      --dbT=<>          : dynamic buoyancy: Time constant
      --dbj=<>          : dynamic buoyancy: Number of derivatives to vanish for forcing function
      --dbk=<>          : dynamic buoyancy: Over sampling ratio for forcing function
      --dbws=<>         : dynamic buoyancy: sampling width 
      --dbreset=<>      : dynamic buoyancy: reset db field
      --dbdirect=<>     : dynamic buoyancy: force logistically towards one for positive inlet velocity
      --dbz=<>          : dynamic buoyancy: adjust so that net buoyancy flux is zero
      --dbD=<>          : dynamic buoyancy: <False> Perform diffusion on dynamic divergence 
      --dblogh=<>       : dynamic buoyancy: <False> Use logistic equation evolution to move towards hb
      --dblogT=<>       : Time scale for logistic equation can be zero
      --dbin=<>         : dynamic buoyancy: <False> Use old inflow velocity condition

 
# locations of forcing regions
      --r0=<>           : Inside radius for adding fluid
      --x0=<>           : Start of velocity uniform region
      --x1=<>           : End of velocity uniform region
      --x2=<>           : Start of density region
      --x3=<>           : end of density region
      --x4=<>           : 
      --x5=<>           : 
      --x6=<>           : 
      --x7=<>           : 
      --ful=<>          : Body force u velocity left
      --fur=<>          : Body force u velocity right
      --wul=<>          : force u velocity left
      --wur=<>          : force u velocity right
      --wvl=<>          : force v velocity left
      --wvr=<>          : force v velocity right
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
      --lbc=<>          : Lower velocity boundary condition 
      --rbc=<>          : Upper velocity boundary condition 
      --blbc=<>         : Bouyancy left boundary condition 
      --brbc=<>         : Bouyancy right boundary condition 
      --slbc=<>         : Sediment left boundary condition 
      --srbc=<>         : Sediment right boundary condition 
      --signals=<>      : Catch signals      
      --parallel=<>     : Use parallel HDF   
  -s, --start_new_files : Start new files while checkpointing 
  -r, --restart         : Restart  
      --rfn=<>          : Restart file name
  -p, --param           : Load from parameter file
      --pfn=<>          : Parameter file name              
      --reset           : reset solve sim time and integrated variables
      --SV=<>           : Sediment sedimentation velocity
      --bV=<>           : Buoyancy sedimentation velocity
      --SVh=<>          : Sediment sedimentation velocity hinder factor
      --bVh=<>          : Buoyancy sedimentation velocity hinder factor
      --Ski=<>          : Stokes number for particle inertial effects viscous term
      --Skp=<>          : Stokes number for particle inertial effects pressure term
      --Skg=<>          : Stokes number for particle density effects
      --dnoise=<>       : amplitude for direct noise forcing
      --noiseL=<>       : length scale for noise forcing
      --noised=<>       : xy, xz, or yz plane to impose noise forcing from stream function
      --noiseT=<>       : time scale for noise forcing
      --inoise=<>       : initial velocity noise
      --bnoise=<>       : initial bouyancy noise
      --snoise=<>       : initial sediment noise
      --sinitial=<>     : initial sediment value
      --bntype=<>       : type of initial bouyancy noise, gaussian or uniform
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
      --cm=<>           : Calculate using centre of mass rather than front position
      --FB=<>           : Fix the buoyancy
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
      --PIDL=<>         : PID controller low pass filter time constant
      --Xmax=<>         : Terminated simulation if X exceeds this
      --Xmin=<>         : Terminated simulation if X drops below this
      --gmin=<>         : minimum g
      --gmax=<>         : maximum g
      --Umin=<>         : minimum U
      --Umax=<>         : maximum U
      --pmss=<>         : Sdb       = dp.get_bool_param(param,'db')
tandard deviation for moving plume source
      --pmsl=<>         : Decay constant for decorrelation of plume source
      --pmzr=<>         : Radial zero buoyancy for bouyancy in plume simulations
      --pmzt=<>         : Top zero buoyancy for bouyancy in plume simulations
      --pmzrt=<>        : Radial-top zero buoyancy for bouyancy in plume simulations
      --pmIu=<>         : Include average u velocity in divergence update
      --topT=<>         : Timescale for evolution of top divergence
      --topd=<>         : Set velocity sink using dynamic evolution of yz divergence
      --topdd=<>        : Set velocity sink using diffuse domain
      --topddM=<>       : Multiplier for diffuse domain sink
      --topr=<>         : Set velocity sink using momentum radius
      --topq=<>         : Set velocity sink evolving radius to reduce vr
      --topu=<>         : Set velocity sink using instanteous u velocity
      --topdrt=<>       : Clip top-radial divergence
      --intdx=<>        : Add d*int(d,dx)=u*dx(u) term to u momentum equation
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
import ffta
import dftf
tinitial=time.time()
import numpy as np
timea = [np.float(time.time())]
ver=np.float(1.32)

from docopt import docopt
import scipy.special as sp
from mpi4py import MPI
import resource
import math
import warnings
import h5dict
import jutil as ju
import serf as se
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


wti0=np.float(time.time())      # Initial wall time

DED_DEBUG_DBZ = os.getenv('DED_DEBUG_DBZ') is not None

debug_fu=np.bool(os.getenv('DED_FU',False))

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
dbreset      = args.pop('--dbreset','False')
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

ju.write_status("Running")
ju.write_files()
ju.node_info()

iparam=copy.copy(param)

scheme   = dp.get_param(param,'scheme')
forcing  = dp.get_param(param,'Forcing')
F        = dp.get_param(param,'Force')
FB       = dp.get_param(param,'FB')
maxdt    = dp.get_param(param,'maxdt')
mindt    = dp.get_param(param,'mindt')
clipu    = dp.get_param(param,'clipu')
clipB    = dp.get_bool_param(param,'clipB')
clipS    = dp.get_bool_param(param,'clipS')
dtjdivyz = dp.get_param(param,'dtjdivyz')
wdivxl   = dp.get_param(param,'wdivxl')
wdivxr   = dp.get_param(param,'wdivxr')

fixI     = dp.get_bool_param(param,'fixI')
divxI    = dp.get_bool_param(param,'divxI')

U1       = dp.get_param(param,'U1')
Qu       = dp.get_param(param,'Qu')
hu       = dp.get_param(param,'hu')
hb       = dp.get_param(param,'hb')
hs       = dp.get_param(param,'hs')
wu       = dp.get_param(param,'wu')
wb       = dp.get_param(param,'wb')
ws       = dp.get_param(param,'ws')

Tx       = dp.get_param(param,'Tx')
Ty       = dp.get_param(param,'Ty')
Tz       = dp.get_param(param,'Tz')
Nx       = dp.get_param(param,'Nx')
Ny       = dp.get_param(param,'Ny')
Nz       = dp.get_param(param,'Nz')

Nc=dict()
Nc['x']=Nx
Nc['y']=Ny
Nc['z']=Nz

dim=np.int(Nx>1) + np.int(Ny>1) +np.int(Nz>1)  

Nxx       = np.int(np.round(Nx*param['AA']))
Nyy       = np.int(np.round(Ny*param['AA']))
Nzz       = np.int(np.round(Nz*param['AA']))
 
H     = dp.get_param(param,'Height')
W     = dp.get_param(param,'Width')
L     = dp.get_param(param,'Length')
R     = dp.get_param(param,'Radius')
Inlet = dp.get_param(param,'Inlet')
q     = dp.get_param(param,'Angle')
g     = dp.get_param(param,'Gravity')
Re    = dp.get_param(param,'Re')
U     = dp.get_param(param,'Velocity')
T     = dp.get_param(param,'Time')
B     = dp.get_param(param,'Buoyancy')
S     = dp.get_param(param,'Sediment')

Wx    = dp.get_param(param,'Wx')
Wy    = dp.get_param(param,'Wy')
Wz    = dp.get_param(param,'Wz')
Wr    = dp.get_param(param,'Wr')
r0    = dp.get_param(param,'r0')
AA    = dp.get_param(param,'AA')
AAJ   = dp.get_param(param,'AAJ')
x0    = dp.get_param(param,'x0')
x1    = dp.get_param(param,'x1')
x2    = dp.get_param(param,'x2')
x3    = dp.get_param(param,'x3')
x4    = dp.get_param(param,'x4')
x5    = dp.get_param(param,'x5')
x6    = dp.get_param(param,'x6')
x7    = dp.get_param(param,'x7')
fuspecial  = dp.get_param(param,'fuspecial')
fuspecial1 = dp.get_param(param,'fuspecial1')
fuspecial2 = dp.get_param(param,'fuspecial2')

noiseL  = dp.get_param(param,'noiseL')
noised  = dp.get_param(param,'noised')
noiseT  = dp.get_param(param,'noiseT')
noise   = dp.get_param(param,'Noise' )
dnoise  = dp.get_param(param,'dnoise')
xn1     = dp.get_param(param,'xn1')
xn2     = dp.get_param(param,'xn2')
xn3     = dp.get_param(param,'xn3')
xn4     = dp.get_param(param,'xn4')

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
topdrt  = dp.get_bool_param(param,'topdrt')
pmIu    = dp.get_bool_param(param,'pmIu')
pful    = dp.get_bool_param(param,'ful')
pfur    = dp.get_bool_param(param,'fur')
pwul    = dp.get_bool_param(param,'wul')
pwur    = dp.get_bool_param(param,'wur')
pwvl    = dp.get_bool_param(param,'wvl')
pwvr    = dp.get_bool_param(param,'wvr')
pwwl    = dp.get_bool_param(param,'wwl')
pwwr    = dp.get_bool_param(param,'wwr')


if U1 is not None: U1S=U1/U
else: U1S = None


if sType=='gc' and Wz is None: Wz=4*H/Nz
if sType=='gc' and Wx is None: Wx=4*L/Nx
    
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
hab      = param.pop('hab',    None)
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


Xmax  = dp.get_param(param,'Xmax')
Xmin  = dp.get_param(param,'Xmin')

PIDX  = dp.get_param(param,'PIDX')
PIDT  = dp.get_param(param,'PIDT')
PIDL  = dp.get_param(param,'PIDL')
PIDG  = dp.get_bool_param(param,'PIDG')
PIDST = dp.get_param(param,'PIDST')
Wr    = dp.get_param(param,'Wr')


PIDD  = dp.get_param(param,'PIDD')
PIDI  = dp.get_param(param,'PIDI')
PIDP  = dp.get_param(param,'PIDP')
PIDS1 = dp.get_param(param,'PIDS1')
PIDS2 = dp.get_param(param,'PIDS2')


dd       = dp.get_bool_param(param,'dd')
ddD      = dp.get_bool_param(param,'ddD')
ddSc     = dp.get_param(param,'ddSc')
ddT      = dp.get_param(param,'ddT')
dddd     = dp.get_bool_param(param,'dddd')
ddserf   = dp.get_bool_param(param,'ddserf')
ddj      = dp.get_param(param,'ddj')
ddk      = dp.get_param(param,'ddk')
ddq      = dp.get_bool_param(param,'ddq')
ddws     = dp.get_bool_param(param,'ddws')
dddf     = dp.get_bool_param(param,'dddf')
ddscl    = dp.get_bool_param(param,'ddscl')
ddlogh   = dp.get_bool_param(param,'ddlogh')
ddlogT   = dp.get_param(param,'ddlogT',0.0)

db       = dp.get_bool_param(param,'db')
dbD      = dp.get_bool_param(param,'dbD')
dbz      = dp.get_bool_param(param,'dbz')
dbu      = dp.get_bool_param(param,'dbu')
dbws     = dp.get_bool_param(param,'dbws')
dbT      = dp.get_param(param,'dbT')
dbj      = dp.get_param(param,'dbj')
dbk      = dp.get_param(param,'dbk')
dbserf   = dp.get_bool_param(param,'dbserf')
dbdirect = dp.get_bool_param(param,'dbdirect')
dbin     = dp.get_bool_param(param,'dbin')
dblogh   = dp.get_bool_param(param,'dblogh')
dblogT   = dp.get_param(param,'ddloghT',0.0)

oddg     = dp.get_bool_param(param,'oddg')
topd     = dp.get_bool_param(param,'topd')
topdd    = dp.get_bool_param(param,'topdd')
topr     = dp.get_bool_param(param,'topr')
topq     = dp.get_bool_param(param,'topq')
topu     = dp.get_bool_param(param,'topu')
topT     = dp.get_param(param,'topT')
topdivr  = dp.get_bool_param(param,'topdivr')

topddM = dp.get_param(param,'topddM')

fpid  = dp.get_param(param,'fpid')
fpidu = dp.get_bool_param(param,'fpidu')
fpidv = dp.get_bool_param(param,'fpidv')
fpidw = dp.get_bool_param(param,'fpidw')
fpidb = dp.get_bool_param(param,'fpidb')
fpids = dp.get_bool_param(param,'fpids')
conservative=dp.get_bool_param(param,'conservative')
rescale=dp.get_bool_param(param,'rescale')

state_ev=dict()  # Extra state variables
avar=dict()      # Auxiliary variables to save if requested with dtjvar



if PIDX is not None:
    state_ev['PIDIT']=0.0
    state_ev['PIDDD']=0.0
    if PIDG: state_ev['g']=dp.get_param(param,'Gravity')
    else:    state_ev['U']=dp.get_param(param,'Velocity')

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


  
fd3=None
tops=None
rmstopd=None
bmult=None
smult=None
wbr=None
wsr=None
wbl=None
wsl=None
wmult=None
bmin=None
bmax=None
smin=None
smax=None
pmif=None

dt = None
Idds=None
Igdw=None
u1=None
u2=None
fuzz=None
    
ju.make_lock()
if signals: ju.init_signals()



param=dp.remove_None(param)
ju.logger.info("Writing parameters")
ju.write_param(param)  

if r0 is None: r0=0
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

gyz = ju.gyz
gyzJ = ju.gyzJ
gx  = ju.gx


ju.logger.info('Calculating Integration weights')
Iw, gIw  = ju.int_weights( 1,domain)
Iww,gIww = ju.int_weights(AA,domain)

Ix  = Iw[0]
Ixx = Iww[0]
Iz  = Iw[-1]
Izz = Iww[-1]

gIx  = gIw[0]
gIxx = gIww[0]
gIz  = gIw[-1]
gIzz = gIww[-1]
mgIzz = (gIzz[0:-1]+gIzz[1:])/2

xx = domain.grids(AA)[ 0]
zz = domain.grids(AA)[-1]
x  = domain.grids( 1)[ 0]  # These will contain only the local grid coordinates
z  = domain.grids( 1)[-1]
xxf=xx.flatten()
if Ny>1:
    dim=3
    yy = domain.grids(AA)[ 1]
    y  = domain.grids( 1)[ 1]
    Iy  = Iw[1]
    Iyy = Iww[1]
    gIy  = gIw[1]
    gIyy = gIww[1]
    yslices  = domain.dist.grid_layout.slices(scales=1)[1]
    yyslices = domain.dist.grid_layout.slices(scales=AA)[1]
    Nyl      = domain.dist.grid_layout.local_shape(scales=1)[1]
    Nyyl     = domain.dist.grid_layout.local_shape(scales=AA)[1]
    gyyslices =ju.comm.gather(yyslices,root=0)
else:
    dim=2
    yy=0
    y=0
    Iy=1
    Iyy=1
    gIy=1
    gIyy=1
    Nyl=1
    Nyyl=1



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

def intzz(f):
    I=0
    if yycolor: I=yycomm.reduce((Izz*f).sum(), op=MPI.SUM,root=0)
    I=ju.comm.bcast(I, root=0)
    return(I)
    
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

gz   = np.zeros( szz,dtype=np.float64)
gzz  = np.zeros(szzz,dtype=np.float64)


if ycolor:
    zcounts  = np.array(ycomm.gather(Nzl,  root=0))
    ycomm.Gatherv(sendbuf=z, recvbuf=(gz, zcounts), root=0)
else: zcounts=np.array([])
if yycolor:
    zzcounts = np.array(yycomm.gather(Nzzl, root=0))
    yycomm.Gatherv(sendbuf=zz, recvbuf=(gzz, zzcounts), root=0)
    zzcounts = yycomm.bcast(zzcounts, root=0)
else: zzcounts=np.array([])

if Ty=='SinCos': Py=1
else:            Py=0
if Tz=='SinCos': Pz=1
else:            Pz=0



      
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
fgzz = gzz.flatten()
fgz  = gz.flatten()

Dgzz = np.min(gIzz)/gIzz # Diffusion weights that conserve the quantity

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


FFTAM=17
nfx=10
nfz=nfx

FFTAx=ffta.ffta(domain.bases[0],FFTAM,nfx)
FFTAy=ffta.ffta(domain.bases[1],FFTAM,nfx)
FFTAz=ffta.ffta(domain.bases[-1],FFTAM,nfx)

# Does the divergence change dynamically
dfd = dd or (sType=='pm' and ('pmss' in param or  forcing==7))
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
        dynoise = jn.jnoise('xyz',sat[noiseL,noiseL,noiseL],1/noiseT,dnoise,domain)
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
nW  = None

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
        if not fbmax: nfb = domain.new_field() # Target b
        if not fbmin or not fbmax:
            nwb = domain.new_field() # Target b
            ju.set_parity(nwb,Tx,Ty,Tz,1,1,1)
#ju.logger.info('forcing {}'.format(forcing))
#ju.logger.info('fbmult {}'.format(fbmult))
#ju.logger.info('fbmin {}'.format(fbmin))
#ju.logger.info('fbmax {}'.format(fbmax))
#ju.logger.info('nfb {}'.format(nfb))
#ju.logger.info('nwb {}'.format(nwb))
#quit(1)
                
    
if sType=='pm':
    ftopu=None
    ftopv=None
    ftopw=None
    ftopuu=None
    if 'pmss' in param: # moving source for the plume
        state_ev['SY']=param['pmss']*np.random.normal()
        state_ev['SZ']=param['pmss']*np.random.normal()
    RL=2*R # length scale for calculating Reynolds number
    IA = math.pi*R**2 # Target area
    CA = math.pi   # Area of unit circle accounting for symmetry
    HL=L
    FL=L
    if Ty=='SinCos':
        IA=IA/2        # Inlet area 
        HW=W           # Half the width of the region
        FW=2*W         
        CA=CA/2
        HNy=2*Ny
    else:
        HW=W/2
        FW=W
        HNy=Ny
    if Tz=='SinCos':
        IA=IA/2
        HH=H
        FH=2*H
        CA=CA/2
        HNz=2*Nz
    else:
        HH=H/2
        FH=H
        HNz=Nz
    RR=min(HW,HH)-r0
    SdAA=dAA*math.pi/CA
    R0=R
    InletA=1
    if ju.mpirank==0:  (R,InletA)=ju.solve_inlet(gyy,gzz,SdAA,Inlet,R,FFTAz)
    R=ju.comm.bcast(R,0)
    InletA=ju.comm.bcast(InletA,0)
    ju.logger.info('Radius changed from {:6.4f} to {:6.4f}, amplitude {:6.4f}'.format(R0,R,InletA))
    ( wr,wrmax, r)=ju.pinlet( y, z,Inlet,R,FFTAz)
    (wrr,wrmax,rr)=ju.pinlet(yy,zz,Inlet,R,FFTAz)
    theta = np.arctan2(yy,zz)
    rrm=np.maximum(1e-10,rr)
    
    ju.check_error()
    wr  = InletA*wr
    wrr = InletA*wrr
    BM  = InletA*wrmax
    SM  = InletA*wrmax
    #Adjust area

    IQ   = ju.jsum(wr     )*dy*dz/IA
    IQQ  = ju.jsum(wr*wr  )*dy*dz/IA
    IIQ  = ju.jsum(wrr    )*dyy*dzz/IA
    IIQQ = ju.jsum(wrr*wrr)*dyy*dzz/IA
    IAU  = IA*U
    IAUU = IA*U**2

    if 'pmss' in param: avar['R']  = np.float64(0)

    if dd:
        avar['u']        = gx.new_field(Px=-1) 
        state_ev['divx'] = gx.new_field(Px=1) 
        sevdivx          = np.zeros((Nxx,), dtype=np.float64)
        avar['Iu']       = np.float64(0.0)
        Hrr=FFTAz.heaviside(rr-RR)[0]
        nwdr   = np.nan_to_num(Hrr/(np.minimum(HW/np.abs(np.sin(theta)),HH/np.abs(np.cos(theta)))-RR));
        Inwdr = ju.jsum(nwdr)*dAA 
        nwdr  = nwdr/Inwdr
        dAAnwdr=dAA*nwdr
        wdrn = 1-Hrr
        if ddq: # Flux dynamic divergence integrate vr over each horizontal region
            #nry = gyz.new_field(Py=-1,Pz= 1) # y component of r unit vector with correct symmetry
            #nrz = gyz.new_field(Py= 1,Pz=-1) # z component of r unit vector with correct symmetry
            #dftf.sawtooth(nry,gyz.domain,'y')
            #dftf.sawtooth(nrz,gyz.domain,'z')
            #nry = nry['g'][0,:,:]/rrm # We can weight our integration towards distances far from the plume
            #nrz = nrz['g'][0,:,:]/rrm
            # Weight the distances so more significance is given to far from the plume and this also ensures that things are smooth towards the origin
            # However we really would like there to be no overlap between this integration and nwdr
            dnry   = (yy*FFTAz.heaviside(RR-rr)[0]*(1-np.exp(-(rr/RR/2)**2)))[0,:,:]
            dnrz   = (zz*FFTAz.heaviside(RR-rr)[0]*(1-np.exp(-(rr/RR/2)**2)))[0,:,:]
            Adnr   = ju.jsum((dnry*yy[0,:,:]+dnrz*zz[0,:,:])/np.maximum(1e-10,rr**2))  # What sum would be get for exactly vr=q/r ?
            ju.logger.info('Dynamic divergence using vr flux integral Adnr: {}'.format(Adnr))
            dnry /= Adnr
            dnrz /= Adnr
            avar['q'] = gx.new_field(Px=1) 
            state_ev['qm']=1 # multiplier for flux
        else:
            avar['dudx']     = gx.new_field(Px=1) 
           
    if topr or topd or topu or topq: state_ev['topd'] = gyz.new_field(1,1)
 
    if ddSc is not None: # Schmidt number for diffusion of dynamic divergence
        Dc  = np.zeros((Nxx+1,),dtype=np.float64)
        Ddx = np.zeros((Nxx+1,),dtype=np.float64)
         
    if topu :
        ftopu = gyz.new_field(Py=1,Pz=1)
        avar['topu']  = ftopu
        
    if topr:
        ftopu  = gyz.new_field(Py=1,Pz=1)
        ftopuu = gyz.new_field(Py=1,Pz=1)
        state_ev['topr'] = param['Radius']
        avar['topdr']    = param['Radius']
        avar['topur']    = param['Radius']
        avar['topIu']    = np.float64(0)
        avar['topIuu']   = np.float64(0)
        avar['topId']    = np.float64(0)
        avar['topIdd']   = np.float64(0)
        avar['topuu']    = gyz.new_field(1,1)
        avar['topu']     = ftopu
        
    if topq:
        state_ev['topr'] = param['Radius']
        avar['topq']     = np.float64(0)
        r3 = np.float64(0)

    if topq or topdivr or dddd: # calculate x and y components of r unit vector making it go to zero at the edges
        oddy = gyz.new_field(Py=-1,Pz= 1)
        oddz = gyz.new_field(Py= 1,Pz=-1)
        dftf.sawtooth(oddy,gyz.domain,'y')
        dftf.sawtooth(oddz,gyz.domain,'z')
        nry = oddy['g']/rrm[0,:,:]
        nrz = oddz['g']/rrm[0,:,:]
        oddy=oddy['g'][:,0]
        oddz=oddz['g'][0,:]
        
    if topd or topq or topdivr or dddd:
        ftopv = gyz.new_field(Py=-1,Pz= 1)
        ftopw = gyz.new_field(Py= 1,Pz=-1)

    if topdivr: avar['topdivr'] = gyz.new_field(Py=1,Pz=1)

    if topd or dddd: avar['topdiv'] = gyz.new_field(Py= 1,Pz=1)
  
    if topdd:
        problem.parameters['ddfx'] = domain.new_field()
        ju.set_parity(problem.parameters['ddfx'],Tx,Ty,Tz,-1, 1, 1)
        problem.parameters['ddfx']['g'] = topddM*FFTAx.df(L -x,2) # exp(-20)=2e-9

    wdivx=1
    if wdivxl is not None: wdivx -= FFTAx.Ndf(  xxf,wdivxl)
    if wdivxr is not None: wdivx -= FFTAx.Ndf(L-xxf,wdivxr)

    if dddd: # Diffuse domain dynamic divergence
        problem.parameters['ddfy'] = domain.new_field()
        problem.parameters['ddfz'] = domain.new_field()
        problem.parameters['ddfy'].set_scales(AA)
        problem.parameters['ddfz'].set_scales(AA)
        ju.set_parity(problem.parameters['ddfy'],Tx,Ty,Tz,1,-1, 1)
        ju.set_parity(problem.parameters['ddfz'],Tx,Ty,Tz,1, 1,-1)
        ddf = dftf.heaviside2(math.pi*yy/HW,math.pi*zz/HH,Ny)/np.maximum(1e-10,rr)
        ju.logger.info('nry.shape {}'.format(nry.shape))
        ju.logger.info('nrz.shape {}'.format(nrz.shape))
        ju.logger.info('ddf.shape {}'.format(ddf.shape))
        ju.logger.info('wdivx.shape {}'.format(wdivx.shape))
        problem.parameters['ddfy']['g'] = nry*ddf*wdivx.reshape((Nxx,1,1))
        problem.parameters['ddfz']['g'] = nrz*ddf*wdivx.reshape((Nxx,1,1))
        avar['div'] = gx.new_field(Px=1) 
        
    if wdivxl is None and wdivxr is None: wdivx=None
    
    if noise is not None:
        wnoise = domain.new_field()  # weighting for noise
        wnoise['g']=ju.stepwf(x,nfx,xn1,xn2,xn3,xn4,L)*wr/wrmax

    if forcing is not None and (B is not None or S is not None):
        if forcing <7: 
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
      
    if forcing>6: # Plume
        # u velocity field is odd
        # g          is constantju
        # b and S are odd
        # Density is given Inlet profile between x6 and x7
        # Density is set to zero between between x4 and x5
        # Fluid is inserted between x2 and x3 and horizontal velocities forced to zero 
        # Fluid is removed  between x0 and x1

        nfd  = domain.new_field() # divergence field
        M=7
        if Tx != 'SinCos':
            ju.logger.info('Tx must be SinCos for forcing 7')
            ju.set_error('Tx must be SinCos for forcing 7')
            ju.test_error()
        (wul,wdl,wdls) = FFTAx.fdfsl(    xx,2)  # Bottom fluid source region             
        (wur,wdr,wdrs) = FFTAx.fdfsl(L - xx,2)  # Top    fluid sink   region          
        #(sdr,sdrs)    = FFTAx.deltasl(xx-xx[wdrs.start],1)  # Sampling for top velocities
        (sdr,sdrs)     = FFTAx.dfsl(L - xx,2)   # Sample same regions as sink 
        
        pwvl = pwvl or pwwl
        pwvr = pwvr or pwwr 

        X0 = xx[wdls].max()

       
        if fixI: # Fix injection profile
            pmif=np.power(np.maximum(xx,X0),-1/5)*FFTAx.heaviside(L-X0-xx,1)[0]*FFTAx.heaviside(xx-X0,1)[0]
            pmif=pmif/np.sqrt((pmif**2).sum()*dxx)

        if pwur or pwul:
            wwls=slice(0,0)
            wwrs=slice(Nxx,Nxx)
            if pwul: (wwl,wwls) = FFTAx.dfsl(      xx,4)  # Forcing for zeroing velocities
            if pwur: (wwr,wwrs) = FFTAx.dfsl(  L - xx,2)  # Forcing for zeroing velocities
            wulrs = slice(wdls.stop,wdrs.start)
            
            wuri = FFTAx.heaviside(xx[wulrs]-xx[wulrs].mean(),2)[0]
            wuli = 1-wuri
            nwu  = domain.new_field() 
            nfu  = domain.new_field()
            ju.set_parity(nwu,Tx,Ty,Tz,1,1,1)
            ju.set_parity(nfu,Tx,Ty,Tz,-1,1,1)
            nwu.set_scales(AA)
            nfu.set_scales(AA)
            if pwul: nwu['g'] += FFTAx.Ndf(  xx, 4)
            if pwur: nwu['g'] += FFTAx.Ndf(L-xx, 2)
            

        if pwvl or pwvr:
            if pwvl==pwul and pwvr==pwur:
                problem.substitutions['wv']   = 'wu'
            else:
                if pwvl: (wwl,wwls) = FFTAx.dfsl(      xx,4)  # Forcing for zeroing velocities
                if pwvr: (wwr,wwrs) = FFTAx.dfsl(  L - xx,2)  # Forcing for zeroing velocities
                nwv  = domain.new_field() 
                ju.set_parity(nwv,Tx,Ty,Tz,1,1,1)
                nwv.set_scales(AA)
                if pwvl: nwv['g'] += FFTAx.Ndf(  xx, 4)
                if pwvr: nwv['g'] += FFTAx.Ndf(L-xx, 2)
            problem.substitutions['ww']   = 'wv'


 
        if oddg: Bwxx =   FFTAx.heaviside(L/2-xx,2)[0]
        else:    Bwxx = 2*FFTAx.heaviside(    xx,2)[0]-FFTAx.heaviside(xx-L/2,2)[0]-1
        
        if nfb is None: nfb = domain.new_field()
        nfb.set_scales(AA)
        if oddg: nfb.meta['x']['parity'] =  1
        else:    nfb.meta['x']['parity'] = -1
        nfb['g'] = wrr*Bwxx # Buoyancy field
        if forcing==7:
            if fbmult:
                if B is not None: problem.substitutions['wb'] = 'wv'
                if S is not None: problem.substitutions['ws'] = 'wv'
                if B is not None: [wbr,wbrs] = FFTA.Ndfsl(L-xx,1)   # Buoyancy top forcing region
                if S is not None: [wsr,wsrs] = FFTA.Ndfsl(L-xx,1)   # Sediment top forcing region
            else:
                nwb.set_scales(AA)
                nwb['g']           = FFTAx.Ndf(xx,4)  # Bottom buoyancy forcing region
                if pmzt:  nwb['g'] = np.maximum(nwb['g'], FFTAx.Ndf(L-xx,2))
                if pmzrt: nwb['g'] = np.maximum(nwb['g'], FFTAx.Ndf(L-xx,2)*Hrr)
                if pmzr:  nwb['g'] = np.maximum(nwb['g'], Hrr)
                if not oddg and not pmzt: nwb['g'] = np.maximum(nwb['g'],FFTAx.Ndf(L-xx,1))
                nwb['g'] = ju.clip01(nwb['g'])
        
    if forcing==6: # Dynamically adjust the divergence so as to make the velocity field correct

        nfu  = domain.new_field() # Target u
        nwu  = domain.new_field() # Weight for u
        nfd  = domain.new_field() # divergence field
 
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

        IAA=ju.jsum(wrr)*dAA
        inwdr=jtheta.periodic_gaussianu(1/(param['Radius']+0.1*(L-x3)),yy,zz,W,H)
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
    sed=se.fu(fgzz,gIzz,H,A=U1S,h=hu,w=wu)
    seb=se.fb(fgzz,gIzz,H,A=B,h=hb,w=wb)
    ses=se.fb(fgzz,gIzz,H,A=S,h=hs,w=ws)
    if ddserf:
        fhu=None
        fwu=None 
    
    if U1 is None: U1=0
    if forcing is not None:
        #ju.logger.info('hu: {}'.format(hu));
        #ju.logger.info('H:  {}'.format(H));
        if hu is not None:
            hu=min(H,max(0,hu))
            RL=hu
        elif hab is not None: RL=hab
        # This is what the forcing should achieve at the Inlet
        problem.parameters['EFB']=0
        problem.parameters['EFS']=0
        problem.parameters['EFU']=-param['Velocity']
        problem.parameters['EFV']=0
        problem.parameters['EFW']=0
        problem.parameters['EFN']=Nxx-1

        U0=-1
        #ju.logger.info('H {} hu {} U {} U0 {} U1 {}'.format(H,hu,U,U0,U1))
        if  hu is None: U2 = U0
        elif   hu <= 0: U2 = U0
        elif   hu >= H: U2 = 0.0
        else:           U2 = (H*U0-hu*U1)/(H-hu);   # h*U1+(H-h)*U2=H*U0
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
        if   hb>=H: hbz  = np.ones_like(z)
        elif hb<=0: hbz  = np.zeros_like(zz)
        else:       hbz  = ju.smoothwf(z,nfz,hb+Wz/2,hb-Wz/2)
        if   hb>=H: hbzz = np.ones_like(zz)
        elif hb<=0: hbzz = np.zeros_like(zz)
        else:       hbzz = ju.smoothwf(zz,nfz,hb+Wz/2,hb-Wz/2)
    else:
        hbz  = ju.smoothwf( z,nfz,H,0)
        hbzz = ju.smoothwf(zz,nfz,H,0)
    if xb is not None:
        ifb=domain.new_field()
        if Tx == "SinCos": ifb.meta['x']['parity'] = 1
        if Ty == "SinCos": ifb.meta['y']['parity'] = 1
        if Tz == "SinCos": ifb.meta['z']['parity'] = 1
        ifb.set_scales(AA)
        if hab is None: ihbzz=hbzz
        else:           ihbzz = ju.smoothwf(zz,nfz,hab+Wz/2,hab-Wz/2)
        if xa is not None: ifb['g'] = ihbzz*(FFTAx.heaviside(xx-xa,1)[0]-FFTAx.heaviside(xx-xb,1)[0])
        else:              ifb['g'] = ihbzz* FFTAx.heaviside(xb-xx,1)[0]
    if Tz=='Cheb':
        if hu is  None:
            hu1=np.zeros_like(z)
        else:
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
        if forcing!=7 and x7 is not None and x6 is not None:
            W76=x7-x6
            nwb['g'] = ju.stepwf(x,nfx,x0,x0+Wx,x6,    x7,      L) # Weight field
            if hb is None: nfb['g'] = 0 
            else: nfb['g'] = hbz*ju.stepwf(x,nfx,x4,x5,   x7+W76,x7+2*W76,L) # Buoyancy field
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
    elif forcing==5:  # velocity forcing only in x0 x2 
 
        # Constant Velocity forcing  from x0,x1 to x2,x3 amd is constant
        # x0 0.0
        # x1 0.5
        # x2 1.0
        # x3 1.5
    
        nwu      = domain.new_field()  # Velocity forcing weight 
        if Ty=="SinCos": nwu.meta['y']['parity'] = 1
        #nwu.set_scales(AA)
        #nwb.set_scales(AA)
        #nfb.set_scales(AA)
        #nwu['g'] =      (FFTAx.heaviside(x1-xx,4)[0]-FFTAx.heaviside(x0-xx,2)[0])
        #nwb['g'] =      (FFTAx.heaviside(x3-xx,4)[0]-FFTAx.heaviside(x2-xx,2)[0])
        #nfb['g'] = hbzz*(FFTAx.heaviside(x5-xx,4)[0]-FFTAx.heaviside(x4-xx,2)[0])

        nf=2
        nb=8

        x8=2*x5-x4
        nwu['g'] = FFTAx.Hn(x,x0,x1,nf)-FFTAx.Hn(x,x2,x3,nb)
        nwb['g'] = FFTAx.Hn(x,x0,x1,nf)-hbz*FFTAx.Hn(x,x4,x5,nb)-(1-hbz)*FFTAx.Hn(x,x1,x6,nb)
        nfb['g'] = FFTAx.Hn(x,x6,x7,nf)-FFTAx.Hn(x,x5,x8,nb)
        #nwb['g'] = hbz*(FFTAx.Hn(x,x0,x1,nf)-FFTAx.Hn(x,x4,x5,nb))+(1-hbz)*(FFTAx.Hn(x,x0,x1,nf)-FFTAx.Hn(x,x1,x6,nb))
        #nwb['g'] = hbz*(FFTAx.Hn(x,x0,x1,2)-FFTAx.Hn(x,x4,x5,4))
        #nwb['g'] = hbz*(FFTAx.Hn(x,x6,x7,2)-FFTAx.Hn(x,x4,x5,4))+(FFTAx.Hn(x,x0,x1,2)-FFTAx.Hn(x,x1,x6,4))

        
        if Ny>1: problem.substitutions['wv'] = 'wu'
        if Nz>1: problem.substitutions['ww'] = 'wu'

        if PIDX is not None and not PIDG:
            nfu = domain.new_field()  # Velocity function
            nfu['g']=-1
        else: problem.substitutions['fu'] = '-U'
        
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
    elif forcing==7: # Gravity current 
        nfd = domain.new_field()  # div(u)
        ju.set_parity(nfd,Tx,Ty,Tz,1,1,1)
        nfd.set_scales(AA)
        if pwul or pwur:
            nfu = domain.new_field()  # Velocity forcing function
            nfu.set_scales(AA)
            ju.set_parity(nfu,Tx,Ty,Tz,-1,1,1)

        
        if Tx=='SinCos': # ful is defined on too small a region
            UL=0 # velocity is zero on left of domain
            if ddj is None:
                mm=2
                (fur,wdr,wdrs) = FFTAx.fdfsl(xx-L,mm)           # Forcing for right divergence
                (ful,wdl,wdls) = FFTAx.fdfsl(xx,mm)             # Forcing for left divergence
                furs=wdrs         # Forcing for right velocity
                fuls=wdls         # Forcing for left velocity
            else:
                ju.logger.info('ddj {} ddk {}'.format(ddj,ddk))
                (wdls,wdl,ful,Nl,dwdl)=FFTAx.deltajks(xx.flatten(),  ddj,0,ddk)
                (wdrs,wdr,fur,Nr,dwdr)=FFTAx.deltajks(-np.flip(xx.flatten()),ddj,0,ddk)
                Nwdl = (wdl.sum()*dxx)
                Nwdr = (wdr.sum()*dxx)
                wdl  = wdl/Nwdl
                wdr  = wdr/Nwdr
                dwdl = dwdl/Nwdl
                dwdr = dwdr/Nwdr
                dwdl = dwdl.reshape((dwdl.shape[0],)+xx.shape[1:])
                dwdr = dwdr.reshape((dwdr.shape[0],)+xx.shape[1:])
                Nful = ful.shape[0]
                Nfur = fur.shape[0]
                ful  = np.concatenate((ful,0.5+np.flip(ful)/2,0.5-ful/2))
                fur  = np.concatenate((-0.5-fur/2,-0.5+np.flip(fur)/2,fur))
                fuls = slice(0,3*Nful)
                furs = slice(Nxx-3*Nfur,Nxx)
                ful  = ful.reshape((ful.shape[0],)+xx.shape[1:])
                fur  = fur.reshape((fur.shape[0],)+xx.shape[1:])
                wdl  = wdl.reshape((wdl.shape[0],)+xx.shape[1:])
                wdr  = wdr.reshape((wdr.shape[0],)+xx.shape[1:])

            sdl  = wdl/(wdl.sum()*dxx)
            sdls = wdls
            wul  = wdl
            wur  = wdr
            wuls = wdls
            wurs = wdrs
            wudl = ful[:wdl.shape[0]]*wdl
            wudr = fur[-wdr.shape[0]:]*wdr
            if dbj is None:
                (wbl,wbls) = FFTAx.dfsl(xx,4)       # Forcing for left buoyancy
                (wbr,wbrs) = FFTAx.dfsl(L-xx,4)     # Forcing for right buoyancy
                (sbl,sbls) = FFTAx.dfsl(xx,4)       # Sampling for buoyancy
                sbl=sbl.reshape((sbl.shape[0],1))/(dxx*sbl.sum())
                (fbl,fbls) = FFTAx.fsl(   xx,4)
                #fbl=FFTAx.heaviside(L/2-xx,2)[0]
            else:
                (wbls,wbl,fbl)=FFTAx.deltajks(xx.flatten()  ,dbj,0,dbk)[0:3]
                (wbrs,wbr)=FFTAx.deltajks(-np.flip(xx.flatten()),dbj,0,dbk)[0:2]
                wbl  = wbl/wbl.max()
                wbr  = wbr/wbr.max()
                wbl  = wbl.reshape((wbl.shape[0],)+xx.shape[1:])
                wbr  = wbr.reshape((wbr.shape[0],)+xx.shape[1:])
                sbl  = wbl.reshape((wbl.shape[0],1))/(wbl.sum()*dxx)   # Sampling for buoyancy
                sbls = wbls
                #fbl  = FFTAx.heaviside(L/2-xx,2)[0]
                fbl  = np.concatenate((np.ones_like(fbl),0.5+np.flip(fbl)/2,0.5-fbl/2))
                fbls = slice(0,fbl.shape[0])
                fbl  = fbl.reshape((fbl.shape[0],)+xx.shape[1:])
            sdl  =  sdl.reshape(( sdl.shape[0],1))
        else:
            UL=-1 # Velocity is -1 on left of domain
            # ful = [__/-\]
            [ful,fuls] = FFTAx.HH(dbj,AA*dbk,(2,1,0)) # Velocity forcing function
            [wul,wuls] = FFTAx.HH(dbj,AA*dbk,(0,1,0)) # Velocity forcing weight
            nWU  = np.int64(ful.shape[0]/5)
            wdls=slice(2*nWU,3*nWU)
            tu   = gx.new_field(Px=0)
            tu['g'][fuls] = ful.copy()
            wdl  = tu.differentiate('x')['g'][wdls]          # desired  du/dx
            dwdl = tu.differentiate('x','x')['g'][wdls]      # desired ddu/dxx

            wudl = wdl*ful[wdls].copy()                      # advective force
            sdl  = wdl/(wdl.sum()*dxx)                       # sampling function
            sdls = wdls


            wdl  =  wdl.reshape(( wdl.shape[0],)+xx.shape[1:])
            sdl  =  sdl.reshape(( sdl.shape[0],)+xx.shape[1:])
            wul  =  wul.reshape(( wul.shape[0],)+xx.shape[1:])
            ful  =  ful.reshape(( ful.shape[0],)+xx.shape[1:])
            wudl = wudl.reshape((wudl.shape[0],)+xx.shape[1:])
            dwdl = dwdl.reshape((dwdl.shape[0],)+xx.shape[1:])

            wdr  = None
            wur  = None
            dwdr = None
            wur   = None

            
            [wbl,wbls] = FFTAx.HH(dbj,AA*dbk,(0,1,0)) # Buoyancy  forcing weight
            [fbl,fbls] = FFTAx.HH(dbj,AA*dbk,(1,1,0)) # Buoyancy  forcing function
            sbl  = FFTAx.Djk(dbj,AA*dbk)              # Buoyancy sampling function
            sbl  = sbl/(sbl.sum()*dxx)
            sbls = slice(2*sbl.shape[0],3*sbl.shape[0])
            
            wbl  =  wbl.reshape(( wbl.shape[0],)+xx.shape[1:])
            sbl  =  sbl.reshape(( sbl.shape[0],1))#+xx.shape[1:])
            fbl  =  fbl.reshape(( fbl.shape[0],)+xx.shape[1:])

            wbr  =  None
            sbr  =  None
            fbr  =  None


            #dtuu  = tu.differentiate('x')['g']
            #ddtuu = tu.differentiate('x','x')['g']
            #tuu=tu['g']
            #ju.write_hdf5_vars('force/tuu.hdf5',('tuu','dtuu','ddtuu','nwul'),locals())
            #quit()

        # find incoming velocity profile on each node
        if hu is not None and U1 is not None:
            if  wu is not None:
                UR=sed.f(zz)
                #(UR,URI,URW)=se.serf1(zz,hu,wu)
                #UR=((1+U1/U)*UR-(1+URI/H*U1/U))/(1-URI/H)
            elif Wz is not None:
                UR=ju.smoothwf(zz,nfz,hu-Wz/2,hu+Wz/2)
                #ju.logger.info('int(UR) {}'.format((gIzz*UR).sum()))
                UR = U1/U-H*(1+U1/U)*UR/intzz(UR).sum()#(H-hu)
        else: UR = 1-4*ju.smoothwf(zz,nfz,0,H)

        UR /= -intzz(UR).sum()/H
        #ju.logger.info('U1 {} H {} hu {} Wz {} UR {}'.format(U1,H,hu,Wz,UR))
        #quit()

        # Collect together to get complete profile on each node
        GUR=np.zeros((Nzz,), dtype=np.float64)
        if yycolor: yycomm.Gatherv(sendbuf=UR, recvbuf=(GUR, zzcounts), root=0) 
        GUR = ju.comm.bcast(GUR,0)
        
        if dd:
            state_ev['dd']     = GUR
            gdw                = np.zeros((Nzz,), dtype=np.float64)
            gww                = np.zeros((Nzz,), dtype=np.float64)
            avar['dwdz']       = np.zeros((Nzz,), dtype=np.float64)
            avar['ww']         = np.zeros((Nzz,), dtype=np.float64)
            ju.logger.info('wbr.shape {} wbrs {}'.format(wbr.shape,wbrs))
        #ju.logger.info('wbr.shape {} wbrs {}'.format(wbr.shape,wbrs))
        #ju.logger.info('wbl.shape {} wbls {}'.format(wbl.shape,wbls))
        nwb.set_scales(AA)
        if wbl is not None: nwb['g'][wbls,...]  = wbl/wbl.max()
        if wbr is not None: nwb['g'][wbrs,...] += wbr/wbr.max()
                    

        if db:
            gbp            = np.zeros((Nzz,), dtype=np.float64)
            state_ev['db'] = np.zeros((Nzz,), dtype=np.float64)
            avar['db']     = np.zeros((Nzz,), dtype=np.float64)
            if yycolor: yycomm.Gatherv(sendbuf=hbzz, recvbuf=(state_ev['db'], zzcounts), root=0)
            state_ev['db'] = ju.comm.bcast(state_ev['db'],0)
            sbp=state_ev['db']
            hbzz=((state_ev['db'][zzslices]).reshape(szzzl))
 
        else:
            nfb.set_scales(AA)
            nfb['g'] = hbzz*FFTAx.heaviside(L/2-xx,2)[0]


        # Give sampling regions a fixed width    
            
            
        if dbws:  # sampling width for dynamic buouyancy 
            sbls  = slice(0,np.max(np.argwhere(xxf<dbws)))
            sbl   = ju.smoothi(xxf[sbls],dbj,dbws,0)[0] # dbj number of vanishing derivative
            sbl   = sbl.reshape((sbl.shape[0],1))/(dxx*sbl.sum())

        if ddws:  # sampling width for dynamic buouyancy 
            sdls  = slice(0,np.max(np.argwhere(xxf<ddws)))
            sdl   = ju.smoothi(xxf[sdls],ddj,ddws,0)[0] # ddj number of vanishing derivative
            sdl   = sdl.reshape((sdl.shape[0],1))/(dxx*sdl.sum())

        if pwvl or pwvr or pwwl or pwwr or pwul or pwur:
            nW = domain.new_field()  
            ju.set_parity(nW,Tx,Ty,Tz,1,1,1)
            nW.set_scales(AA)
            nW['g']=0
            if pwul or pwvl or pwwl: nW['g'][wuls,...] += wul/wul.max()
            if pwur or pwvr or pwwr: nW['g'][wurs,...] += wur/wur.max()

        if pwul or pwur: problem.substitutions['wu'] = 'W'
        if pwvl or pwvr: problem.substitutions['wv'] = 'W'
        if pwwl or pwwr: problem.substitutions['ww'] = 'W'
        
        if wdl is not None: nfd['g'][wdls,...] = wdl*(UR-UL)
        if wdr is not None: nfd['g'][wdrs,...] = wdr

        if nfu is not None: nfu['g'][fuls,...] = UL+(UR-UL)*ful

        if pful or pfur:
            rfu = domain.new_field()  
            ju.set_parity(rfu,Tx,Ty,Tz,-1,1,1)
            rfu.set_scales(AA)
            DU=U*(UR-UL)
            if Tx=='SinCos': DUU=0
            else:            DUU=intzz(DU**2)/H
            if conservative: pfum=2
            else:            pfum=1
            if pful: rfu['g'][wdls,...] = pfum*U*wdl*DU*UL + pfum*wudl*(DU**2-DUU)-dwdl/Re*DU
            if pfur: rfu['g'][wdrs,...] =                    pfum*wudr* U**2-dwdr/Re*U

            
        # This forcing approach forces the horizontal velocity with a constant value, -U,
        # The Buoyancy must be large enough to overcome this
        # Force the velocity and buoyancy to zero in a region behind the moving current
    elif forcing==8:
        nwu = domain.new_field()  # Velocity forcing weight 
        if Ny>1: problem.substitutions['wv'] = 'wu'
        problem.substitutions['ww'] = 'wu'
        if B is not None: problem.substitutions['wb'] = 'wu'
        if S is not None: problem.substitutions['ws'] = 'ws'
        nwb=None
        nfb=None
    
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
    else:            nfd.meta['x']['parity'] = 0
    if Ny>1:
        if Ty=="SinCos": nfd.meta['y']['parity'] = 1
        else:            nfd.meta['y']['parity'] = 0
    if Tz=="SinCos": nfd.meta['z']['parity'] = 1
    else:            nfd.meta['z']['parity'] = 0

 
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
ju.set_wfun(problem, nW,     'W'     , L*H*W, Tx, Ty, Tz, AA, JA, AAJ)

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
rfu = ju.set_ffun(problem,rfu,'rfu',U**2,Tx,Ty,Tz,Uflg, AA,-1,JA,AAJ)
for a in problem.parameters:
    pp=problem.parameters[a]
    if ju.jtype(pp)=='fxyz':
        ju.logger.info('{} {} {}'.format(a,type(pp),type(pp.meta)))
        for c in pp.meta:
            #ju.logger.info('{} {} {} {}'.format(a,c,type(pp.meta[c]),pp.meta[c]))
            if not 'parity' in pp.meta[c]: pp.meta[c]['parity']=0

 
    
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
    elif forcing==8: fluxd='yz'
    else:
        ju.logger.error('forcing should be 0 to 8 for plumes'.format(forcing))
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

if clipB:   ju.logger.info('clipB: BM:{:5.2f} [{:5.2f}, {:5.2f}]'.format(BM,minB,maxB))
if clipS:   ju.logger.info('clipS: SM:{:5.2f} [{:5.2f}, {:5.2f}]'.format(SM,minS,maxS))
if clipu is not None: ju.logger.info("clipu: {:4.0f},".format(clipu))
ju.logger.info("restart:            {}".format(restart))
ju.logger.info("reset:              {}".format(reset))
ju.logger.info("start_new_files:    {}".format(start_new_files))
if dtcheck is not None:  ju.logger.info("dtcheck:            {:7.4f} (min)".format(dtcheck))
#ju.logger.info("Wall time:          {}".format(param['WallTime']))
ju.logger.info("mindt:  {:8.6f}, maxdt: {:8.6f}".format(mindt, maxdt))
ju.logger.info(ju.nns( 'dx:    {:7.5f}, ',dx)    + ju.nns( 'dy:    {:7.5f}, ',dy)    + ju.nns( 'dz:    {:7.5f}, ',dz)    + ju.nns( 'dr:    {:7.5f}',dr)) 
ju.logger.info(ju.nns( 'Wx:    {:7.5f}, ',Wx)    + ju.nns( 'Wy:    {:7.5f}, ',Wy)    + ju.nns( 'Wz:    {:7.5f}, ',Wz)    + ju.nns( 'Wr:    {:7.5f}',Wr)) 
ju.logger.info(ju.nnsd('Wx/dx: {:7.5f}, ',Wx,dx) + ju.nnsd('Wy/dy, : {:7.5f}',Wy,dy) + ju.nnsd('Wz/dz: {:7.5f}, ',Wz,dz) + ju.nnsd('Wr/dr: {:7.5f}',Wr,dr)) 
if Tdiv: ju.logger.info("Tdiv:  {:7.5e}".format(Tdiv))

#if fsmin:  ju.logger.info("Minimum sediment enforced by min")
#if fsmax:  ju.logger.info("Maximum sediment enforced by max")
#if dd:   ju.logger.info("Dynamic divergence")
#if db:     ju.logger.info("Dynamic bouyancy")
#if noise is not None:
#    ju.logger.info("noise:  {:7.4f}, noiseL: {:7.5f}, noiseT: {:7.5f}, noised {}".format(noise,noiseL,noiseT,noised))
#    ju.logger.info("xn1:    {:7.4f}, xn2:     {:7.4f}, xn3:    {:7.4f},  xn4: {:7.4f}".format(xn1,xn2,xn3,xn4))
#if U1 is not None: ju.logger.info("U1:     {:7.4f}".format(U1))
#if U2 is not None: ju.logger.info("U2:     {:7.4f}".format(U2))
#if pwul:           ju.logger.info("Forcing u velocity left")
#if pwur:           ju.logger.info("Forcing u velocity right")
#if pwvl:           ju.logger.info("Forcing v velocity left")
#if pwvr:           ju.logger.info("Forcing v velocity right")
#if pwwl:           ju.logger.info('Forcing w velocity left')
#if pwwr:           ju.logger.info('Forcing w velocity right')
#if rescale:        ju.logger.info('Rescaling velocities and pressure')

if sType=='pm':
    ju.logger.info("Inlet shape:    {:}".format(Inlet))
    ju.logger.info("Inlet radius:   {:7.4f}".format(R))
    ju.logger.info("Inlet velocity: {:7.4f}".format(U))
    ju.logger.info("Inlet area           IA: {:7.4f}, sqrt(IA/pi): {:7.4f}".format(IA,np.sqrt(IA/np.pi)))
    ju.logger.info("Inlet flux:         {:7.4f} {:7.4f}".format(IQ,IIQ))
    ju.logger.info("Inlet flux squared: {:7.4f} {:7.4f}".format(IQQ,IIQQ))
    if forcing==0: ju.logger.info("TQ1:   {:7.4f} TQ2: {:7.4f} TQ3: {:7.4e}".format(TQ1,TQ2,TQ3))
if pmzt:  ju.logger.info('Buoyancy is set to zero at top');
if pmzr:  ju.logger.info('Buoyancy is set to zero at edges');
if pmzrt: ju.logger.info('Buoyancy is set to zero at top edges');

#for s in problem.substitutions: ju.logger.info('substitutions: {:7s}: {}'.format(s,problem.substitutions[s]))
#for s in problem.parameters:    ju.logger.info('parameters:    {:7s}: {}'.format(s,problem.parameters[s]))

if 'parity' not in problem.meta['u']['x']: problem.meta['u']['x']['parity']=0

if sType=='gc' and U is not None:
    if Tz=='Cheb' and (param['lbc']=='noslip' or param['rbc']=='noslip'):
        problem.parameters['UBC']     = domain.new_field()
    
        for a in problem.meta['u']:
            if 'parity' in problem.meta['u'][a]:
                problem.parameters['UBC'].meta[a]['parity'] = problem.meta['u'][a]['parity']
                problem.parameters['UBC'].meta['z']['constant'] = True
                problem.parameters['UBC']['g'] = -U*FFTAx.odd_constant(x)
    else: problem.parameters['UBC']=-U



if sType=='bg': ju.burgers_equations(domain,problem,param)
else:           ju.boussinseq_equations(domain,problem,param,FFTAx)
#ju.print_parity(problem)

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
            problem.add_bc(         "left(uz)    = 0")
            if Ny>1: problem.add_bc("left(vz)    = 0")
        elif param['lbc']=='noslip':
            if PIDX is not None and not PIDG: problem.add_bc("left(u)     = left(fu)")
            elif U is not None:               problem.add_bc("left(u)     = UBC")
            else:                             problem.add_bc("left(u)     = 0")
            if Ny>1:                          problem.add_bc("left(v)     = 0")
        if   param['rbc']=='slip':
            problem.add_bc("right(uz)   = 0")
            if Ny>1:                          problem.add_bc("right(vz)   = 0")
        elif param['rbc']=='noslip':
            if PIDX is not None and not PIDG: problem.add_bc("right(u)    = right(fu)")
            elif U is not None:               problem.add_bc("right(u)    = UBC")
            else:                             problem.add_bc("right(u)    = 0")
            if Ny>1:                          problem.add_bc("right(v)    = 0")

        problem.add_bc("right(w)    = 0")
        if Ny>1: 
            problem.add_bc("left(w)     = 0", condition="(nx != 0) or (ny != 0)")
            problem.add_bc("integ_z(p)  = 0", condition="(nx == 0) and (ny == 0)")
        else:
            problem.add_bc("left(w)     = 0", condition="(nx != 0)")
            problem.add_bc("integ_z(p)  = 0", condition="(nx == 0)")
                
        vz=list(set(problem.variables) & set(('bz','fpidbz','fpidsz','fpiduz','fpiduz','fpidwz')))
        for v in vz:
            if v != 'sz' and v!= 'bz':
                problem.add_bc("right(" + v + ") = 0")
                problem.add_bc("left("  + v + ") = 0")
            problem.meta[v]['z']['dirichlet'] = True

        if 'sz' in problem.variables:
            if   param['srbc']=='z':  problem.add_bc("right(s)    = 0")
            elif param['srbc']=='dz': problem.add_bc("right(dzs)  = 0")
            elif param['srbc']=='dq': problem.add_bc("kappas*right(dzs)-sVz*right(s) = 0")
            if   param['slbc']=='z':  problem.add_bc("left(s)     = 0")
            elif param['slbc']=='dz': problem.add_bc("left(dzs)   = 0")
            elif param['slbc']=='dq': problem.add_bc("kappab*left(dzb)-sVz*left(s) = 0")

        if 'bz' in problem.variables:
            if   param['brbc']=='z':  problem.add_bc("right(b)    = 0")
            elif param['brbc']=='dz': problem.add_bc("right(dzb)  = 0")
            elif param['brbc']=='dq':
                if 'bVh' in param: problem.add_bc("kappab*right(dzb)=bVz*right(b*(1-bVh*b))")
                else:              problem.add_bc("kappab*right(dzb)-bVz*right(b) = 0")
            if   param['blbc']=='z':  problem.add_bc("left(b)     = 0")
            elif param['blbc']=='dz': problem.add_bc("left(dzb)   = 0")
            elif param['blbc']=='dq':
                if 'bVh' in param: problem.add_bc("kappab*left(dzb)=bVz*left(b*(1-bVh*b))")
                else:              problem.add_bc("kappab*left(dzb)-bVz*left(b) = 0")

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
    if dbreset:
        state_ev['db'] = ju.comm.bcast(sbp,0)
        ju.logger.info('Resetting density profile')
    if 'divx' in  state_ev: state_ev['dd'] = pop(state_ev,'divx')
    if 'divz' in  state_ev: state_ev['dd'] = pop(state_ev,'divz')
        
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
ju.write_hdf5_vars('force/fx.hdf5',('dudxrn','pmif','Bwxx','wmult','wbls','wbl','wbr','wbrs','fbl','fbls','sbl','sbls','sbr','wdl','wdls','sdr','sdrs','sdl','sdls','dwdr','dwdl','wdr','wdrs','wudl','wudr','wul','wuls','wur','wux1','wux2','wux3','wux4','wux5','wux6','wux7','wux8','Iux1','Iux2','Iux3','Iux4','Iux5','sminf','bminf','wwl','wwr','ful','fur','fuls','furs','wdivx','wuli','wuri','GUR','UL'),locals())
#ju.logger.info('MM {}, M {}, xx[MM] {}, hb {}'.format(MM,M,xx[MM],hb))
#ju.logger.info('wux3[0:MM] {}'.format(wux3[0:MM]))

if debug_flow: ju.logger.info('ded_gc: saved yz') 
JA.saved('force/fyz.hdf5',('inwdr','nwdr','wdrn','rr','wrr','nyz','nry','nrz','dnry','dnrz'),'yz',locals(),AA)

   
del(ifb,ifs,ifu,ifv,ifw)

if dim==3:
    for a in set(('ddf','ddfx','ddfy','ddfz')) & set(problem.parameters.keys()):
        JA.saved('force/ddf.hdf5',(a,),'xyz',(problem.parameters[a],),1,'a')
        
    if 'oddcx' in problem.parameters: JA.saven('force/oddcx.hdf5','oddcx',problem.parameters['oddcx']['g'][:,0,0].flatten(),AA,'x')
    if 'oddcy' in problem.parameters: JA.saven('force/oddcy.hdf5','oddcy',problem.parameters['oddcy']['g'][0,:,0].flatten(),AA,'y')
    if 'oddcz' in problem.parameters: JA.saven('force/oddcz.hdf5','oddcz',problem.parameters['oddcz']['g'][0,0,:].flatten(),AA,'z')
else:
    if 'oddcx' in problem.parameters: JA.saven('force/oddcx.hdf5','oddcx',problem.parameters['oddcx']['g'][:,0].flatten(),AA,'x')
    if 'oddcz' in problem.parameters: JA.saven('force/oddcz.hdf5','oddcz',problem.parameters['oddcz']['g'][0,:].flatten(),AA,'z')



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
    
if 'WallTime' in param: solver.stop_wall_time = ju.get_sec(param['WallTime'])
else: solver.stop_wall_time=np.inf
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
    #    if ju.mpirank==0:
    if 'PIDIT' in param:
        state_ev['PIDIT']=param.pop('PIDIT',0.0)
        ju.logger.info('Setting PIDIT to {:7.3f}'.format(state_ev['PIDIT']))
    if 'PIDDD' in param:
        state_ev['PIDDD']=param.pop('PIDDD',0.0)
        ju.logger.info('Setting PIDDD to {:7.3f}'.format(state_ev['PIDDD']))
        
    if 'g' in state_ev and 'Gravity' in param:
        state_ev['g']=param.pop('Gravity')
        ju.logger.info('Setting g to {:7.3f}'.format(state_ev['g']))
    if 'U' in state_ev and 'Velocity' in param:
        state_ev['U']=param.pop('Velocity')
        ju.logger.info('Setting U to {:7.3f}'.format(state_ev['U']))
        
    if state_ev['PIDIT']==0:
        if PIDG: state_ev['PIDIT']=state_ev['g']
        else:    state_ev['PIDIT']=state_ev['U']
    if PIDG and state_ev['g']>0: g=state_ev['g']
    elif state_ev['U']>0:        U=state_ev['U']
    # Must not broadcast state_ev  as then reference in jrec will not update
    #for a in ('PIDIT','g','U','PIDDD'): state_ev[a]=ju.comm.bcast(state_ev[a],0)

    if PIDG: ju.logger.info("PID controller for front position is active on g: [{:7.3f},{:8.4f},{:7.3f}]".format(param['gmin'],state_ev['g'],param['gmax']))
    else:    ju.logger.info("PID controller for front position is active on U: [{:7.3f},{:8.4f},{:7.3f}]".format(param['Umin'],state_ev['U'],param['Umax']))
    ju.logger.info("PIDX: {:7.4f},  PIDIT: {:7.4f}, PIDDD: {:7.1e}".format(PIDX, state_ev['PIDIT'], state_ev['PIDDD']))
    #ju.logger.info("PIDT: {},   PIDP: {}   PIDI: {},  PIDD: {}".format(PIDT, PIDP, PIDI, PIDD))
    ju.logger.info("PIDT: {:7.4f},   PIDP: {:7.4f}   PIDI: {:7.5f},  PIDD: {:7.4f}".format(PIDT, PIDP, PIDI, PIDD))
    if ju.mpirank==0:
        if PIDG: pid = jpid.pid(PIDP, PIDI, PIDD, PIDL, PIDX, state_ev['PIDIT'], state_ev['PIDDD'],  state_ev['g'], st0, X, param['gmin'],param['gmax'],5,-1)
        else:    pid = jpid.pid(PIDP, PIDI, PIDD, PIDL, PIDX, state_ev['PIDIT'], state_ev['PIDDD'],  state_ev['U'], st0, X, param['Umin'],param['Umax'],0.5,1)
     
UF=None
cRe=None

dforce=True

del(wnoise,nwd,nwb,nws,nwu,nwv,nww)
if PIDX is None:
    del(nfb,nfd,nfu,nfv,nfw)
else:
    if fb is None: del(nfb)
    if fs is None: del(nfs)
    if fd is None: del(nfd)
    if fu is None: del(nfu)
    if fv is None: del(nfv)
    if fw is None: del(nfw)
    
if False:
    maxu=ju.jmax(np.abs(u['g']))
    if Ny>1: maxv=ju.jmax(np.abs(v['g']))
    else:    maxv=0
    maxw=ju.jmax(np.abs(w['g']))
    ju.logger.info('maximum velocities: [{:6.3f},{:6.3f},{:6.3f}]'.format(maxu,maxv,maxw))

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

if dd:   fd = problem.parameters['fd']
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

#if 'dda' in state_ev:          state_ev['dd'][:,0] = state_ev.pop('dda')
#if 'ddb' in state_ev and dd: state_ev['dd'][:,1] = state_ev.pop('ddb')
#JA.save_state(dt,state_ev,'state3.hdf5')

if sType=='gc' and not dd: fd = None


if rescale:
    if dim==2: vsqr = de.operators.integrate(u*u+w*w,     'x',      'z')
    else:      vsqr = de.operators.integrate(u*u+v*v+w*w, 'x', 'y', 'z')

iteration_one=True


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
        fflag = ju.reduce_bool(fflag)
        if fflag:
            ju.set_error('Fields not finite')
            ju.logger.info('Fields not finite')
            JA.save_state(dt,state_ev)
        ju.check_error()
        

        

    if forcing==7 and sType=='gc':
        if db:
            b.set_scales(AA)
            #b['g']=np.maximum(0,np.minimum(1,b['g']))
            #JY.print()
            #ju.srange('bp1 [{:6.3f},{:6.3f}]',bp)
            #ju.srange('b   [{:6.3f},{:6.3f}]',b['g'])
            bp = JY.intnr(b['g'],1/(AA*W)) 
            if yycolor:
                bp = dxx*np.sum(sbl*bp[sbls,:],axis=(0,))
                yycomm.Gatherv(sendbuf=bp, recvbuf=(gbp, zzcounts), root=0)
                #print('bp:[{} {}]'.format(ju.jmin(bp,yycomm),ju.jmax(bp,yycomm)))
               
            if ju.mpirank==0:
                #print('gbp: [{:6.3f},{:6.3f}]'.format(gbp.min(),gbp.max()))
                #gbp=np.maximum(0,np.minimum(1,gbp))
                avar['db']=gbp
                if dbin: # use old inflow conditoin
                    if dd:   sf = ju.db_forcing(state_ev['dd'].flatten())
                    else:    sf = ju.db_forcing(GUR.flatten())
                    gbp = sf + (1-sf)*gbp # Make sure density for inflowing region is set to 1

                state_ev['db'] =  (dbT*state_ev['db']+dt*gbp)/(dbT+dt)

                if dbdirect: # Evolve state towards B in inflow regions
                    if dd: sf = state_ev['dd'].flatten()
                    else:  sf = GUR.flatten()
                    state_ev['db'] += dt/(dbT+dt)/sf.max()*np.maximum(sf,0)*(B-state_ev['db'])

                if dbserf:
                    seb.fit(state_ev['db'])
                    state_ev['db']=seb.f(fgzz)

                if dblogh: # logistic evolution towards hb
                    state_ev['db'] = ju.logistic_h(state_ev['db'],hb,gIzz,H,0,B,scale=False,a=dblogT/dt)
 
                if dbD: # Diffuse dynamic divergence according to degree of non-monotonicity
                    state_ev['db'] = ju.tvd_diffusion(state_ev['db'],Dgzz)

                if dbz: # Adjust density profile to zero flux by logistic evolution
                    if dd:   sf = gIzz*state_ev['dd'].flatten()
                    else:    sf = gIzz*GUR.flatten()
                    czbu = (sf*state_ev['db']).sum()/(sf*state_ev['db']*(1-state_ev['db']/B)).sum()
                    czbu = max(-1/2,min(1/2,czbu))
                    state_ev['db'] = state_ev['db']/(1+czbu*(1-state_ev['db']/B))
                    if DED_DEBUG_DBZ: ju.logger.info('czbu {:5.1e} bu {:5.1e} bh {:8.6f}'.format(czbu,(sf*state_ev['db']).sum(),(gIzz*state_ev['db']).sum()))

            state_ev['db'] = ju.comm.bcast(state_ev['db'],0)
            bps=(state_ev['db'][zzslices]).reshape(szzzl)
            problem.parameters['fb']['g'][fbls,...] = bps*fbl
          
        if dd:
            
            #gyz.intx(w,wdl,wleft,wdls)
            #wleft.differentiate('z')
            #WA=JY.intnr(w.differentiate('z')['g']+w.differentiate(x=2).integrate('z')['g'])
            #ju.srange('avar[dwdz] [{:6.3f},{:6.3f}]',avar['dwdz'])
            #ju.srange('avar[bp]   [{:6.3f},{:6.3f}]',avar['db'])
            w.set_scales(AA)    
            u.set_scales(AA)    
            wz.set_scales(AA)    
            #WA=-JY.intnr(w.differentiate('z')['g'])
            WA=-JY.intnr(solver.state['wz']['g']) 
            WW =JY.intnr(w['g']**2)
            if yycolor:
                dw  = dxx/AA*np.sum(sdl*WA[sdls,:],axis=(0,)) # make sure sdl only has two dimensions
                ww  = dxx/AA*np.sum(sdl*WW[sdls,:],axis=(0,))
                yycomm.Gatherv(sendbuf=dw, recvbuf=(gdw, zzcounts), root=0)
                yycomm.Gatherv(sendbuf=ww, recvbuf=(gww, zzcounts), root=0)
            if ju.mpirank==0:
                avar['ww'] = gww/W
                avar['dwdz'] = -gdw/W
                state_ev['dd'] += gdw*dt/(dt+ddT)
                if dbu: # project onto mean velocity zero and zero mass flux
                    Ibb = (gIzz*state_ev['db']*state_ev['db']).sum()
                    Ibu = (gIzz*state_ev['db']*state_ev['dd']).sum()
                    Ib  = (gIzz*state_ev['db']).sum()
                    Iu  = (gIzz*state_ev['dd']).sum()
                    qbu =0 # required mass flux
                    state_ev['dd'] += (Ib*Ibu - Ibb*Iu - H*Ibb -qbu*Ib-(H*Ibu - Ib*Iu - H*Ib-H*qbu)*state_ev['db'])/(H*Ibb - Ib**2)
                    Ibu2 = (gIzz*state_ev['db']*state_ev['dd']).sum()
                    Iu2  = (gIzz*state_ev['dd']).sum()
                    #state_ev['dd'] -= 1
                    #state_ev['dd'][fgzz<=hb]=np.maximum(state_ev['dd'][fgzz<=hb],0)
                    Ibu3 = (gIzz*state_ev['db']*state_ev['dd']).sum()
                      
                    ju.logger.info('hu: {:7.4f} wu: {:7.4f} Iu2: {:7.4f} Ibu: {:7.4f} Ibu2: {:7.4f} Ibu3: {:7.4f} U1: {:7.4f}'.format(hu,wu,Iu2,Ibu,Ibu2,Ibu3,state_ev['dd'][0]))
                   
                if ddserf:
                    #ju.write_hdf5('dd1.h5',('dd'),(state_ev['dd']))
                    #sedf.fit(state_ev['dd'])
                    (fhu,fwu,fau,fbu,yy)=se.fit_serf(fgzz,state_ev['dd'],gIzz,fhu,fwu)
                    # Now enforce height and width constraints and finally constant volume
                    if hu  is not None: fhu=hu
                    if wu  is not None: fwu=wu
                    fau=min(5,max(0.01,fau))
                    serff0=se.serff(0,fhu,fwu)
                    serffH=se.serff(H,fhu,fwu)
                    serfFH=se.serfF(H,fhu,fwu)/H
                    serfu0=fbu+fau*serff0
                    fbu = -(serff0 + serfFH*serfu0)/(serff0 - serfFH)
                    fau = +(serfu0 + 1)/(serff0 - serfFH)
                    state_ev['dd']=fbu+fau*se.serff(fgzz,fhu,fwu) #ju.write_hdf5('dd2.h5',('dd'),(state_ev['dd']))

                    #state_ev['hu']=fhu
                    #state_ev['wu']=fwu
                    #A+B*se.serff(0,fhu,fwu)=sserfu0
                    #A+B*se.serfF(H,fhu,fwu)=-1
                    #A=-1-B*se.serfF(H,fhu,fwu)
                    #state_ev['u0']=fbu+fau*serff0
                    #state_ev['u1']=fbu+fau*serffH
                    #ju.write_hdf5('dd2.h5',('dd'),(dd2))
                    #ju.write_hdf5('dd3.h5',('dd'),(dd3))
                    #print('dd [{:8.5} {:8.5} {:8.5}]'.format(np.min(state_ev['dd']),np.max(state_ev['dd']),fbu*H+fau*se.serfF(H,fhu,fwu)))
                    #sed.print()
                     
                if ddlogh: # logistic evolution If we use min and max to scale to between 0 and 1 then we stay within these limits
                    state_ev['dd'] = ju.logistic_h(state_ev['dd'],hu,gIzz,H,scale=True,a=ddlogT/dt)
              
                if False:
                    ju.srange('state_ev[dd] [{:6.3f},{:6.3f}]',state_ev['dd'])
                    ju.srange('gdw            [{:6.3f},{:6.3f}]',gdw)
                    ju.srange('w              [{:6.3f},{:6.3f}]',dxx*(wdl*solver.state['w']['g'][wdls,...]).sum(axis=(0,)))
                    ju.srange('wz             [{:6.3f},{:6.3f}]',dxx*(wdl*solver.state['wz']['g'][wdls,...]).sum(axis=(0,)))
                    ju.srange('w              [{:6.3f},{:6.3f}]',solver.state['w']['g'][wdls,...])
                    ju.srange('wz             [{:6.3f},{:6.3f}]',solver.state['wz']['g'][wdls,...])
                    ju.srange('wdl            [{:6.3f},{:6.3f}]',wdl)
                    ju.logger.info('{:7.4f} {:7.4f}'.format(wdl.sum()*dxx,dxx))
                    ju.write_hdf5('ww.hdf5',('w','wz','start','stop','wdl','gdw','dxx'),(solver.state['w']['g'],solver.state['wz']['g'],wdls.start,wdls.stop,wdl,gdw,dxx))
                    quit()
                    
                if False: # Fix the minimum value to U2=H/(H-hu)*U
                    U2=H/(H-hu)*U
                    # u=-U2+a*d
                    # -H*U_2 + a*I = -H U
                    # a = H (U_2-U)/I
                    state_ev['dd'] = state_ev['dd']-min(state_ev['dd'])
                    Idd=(gIzz*state_ev['dd']).sum()
                    if Idd!=0: state_ev['dd'] = state_ev['dd']*H*(U2-U)/Idd
                    state_ev['dd'] -= U2
                    #ju.logger.info('{:7.4f} {:7.4f}  {:7.4f}'.format(Idd/H,max(state_ev['dd']),min(state_ev['dd'])))
                if False: #Qu is not None:
                    Qf=state_ev['dd']                       #Qf=state_ev['dd']                       
                    QF0=1 +(           Qf   *gIzz).sum()/H    #QF0=(           Qf   *gIzz).sum()   # Deviation from satisying mean criterion    
                    QF1=Qu-(np.maximum(Qf,0)*gIzz).sum()/H    #QF1=(np.maximum(Qf,0)*gIzz).sum()   # Deviation from satisying positive flux criterion    
                    Qh= (np.float64(Qf>0)*gIzz).sum()/H       #Qh= (np.float64(Qf>0)*gIzz).sum()         
                    QQ=max(0.5,Qu-QF1-QF0*Qh+Qh)              #Qh=max(0.01*H,Qh)                         
                    QA=(QF1 - Qu*QF0)/QQ                      #QA=(QF1*H*U+QF0*Qu)/(QF0*Qh-QF1*H)        
                    QB=(QF1 + Qh*QF0)/QQ                      #QB=(U*Qh+Qu)*H/(-QF0*Qh+QF1*H)            
                    state_ev['dd']+=QA+QB*Qf                #state_ev['dd']=QA+QB*Qf                 
                    ju.logger.info('Qu {} {}'.format(Qu,(np.maximum(state_ev['dd'],0)*gIzz).sum()/H))
                if False: #U1 is not None:
                    # Fix the bottom value to U1
                    state_ev['dd'] -= state_ev['dd'][0]
                    Idd=(gIzz*state_ev['dd']).sum()
                    if Idd!=0: state_ev['dd'] *= -H*(U1/U+1)/Idd
                    state_ev['dd'] += U1
                if ddD: # Diffuse dynamic divergence according to degree of non-monotonicity
                    state_ev['dd']=ju.tvd_diffusion(state_ev['dd'],Dgzz)
                     
                maxdd=state_ev['dd'].max()
                if maxdd<0: state_ev['dd'] -= maxdd
                Idd=-(gIzz*state_ev['dd']).sum()/(H)
                state_ev['dd']=state_ev['dd']/max(0.5,Idd)

                if Tx=='Fourier': DUU=(gIzz*(state_ev['dd']-UL)**2).sum()/H
                if not np.all(np.isfinite(state_ev['dd'])) or (state_ev['dd']).max()>10 or (state_ev['dd']).min()<-10:
                    ju.logger.info('dd not finite')
                    ju.logger.info('Idd {:6.3f} {:6.3f}'.format((gIzz*state_ev['dd']).sum(),Idd))                   
                    ju.logger.info('gdw range [{:6.3f},{:6.3f}]',np.min(gdw),np.max(gdw))
                    ju.logger.info('dd  range [{:6.3f},{:6.3f}]',np.min(state_ev['dd']),np.max(state_ev['dd']))
                    ju.set_error('dd not finite')
                    #ju.srange('bp  range [{:6.3f},{:6.3f}]',state_ev['db'])
                    #ju.logger.info('Ibp {:6.3f} hb {:6.3f} min {:6.3f} max {:6.3f}'.format(Ibp,hb,state_ev['db'].min(),state_ev['db'].max()))
                #ju.logger.info('Idd {:6.3f} {:6.3f}'.format((gIzz*state_ev['dd']).sum(),Idd))
                               
            ju.check_error()
  
                
            state_ev['dd'] = ju.comm.bcast(state_ev['dd'],0)
            DUU            = U*ju.comm.bcast(DUU,0)

            fd['g'][wdls,...]  = wdl*state_ev['dd'][zzslices].reshape(szzzl)
            UR=state_ev['dd'][zzslices].reshape(szzzl)
            DU=U*(UR-UL)

            if wdl is not None: fd['g'][wdls,...]  = wdl*DU
            if wdr is not None: fd['g'][wdrs,...]  = wdr*U
            #if pful: problem.parameters['rfu']['g'][wdls,...]  = ful*state_ev['dd'][zzslices].reshape(szzzl)**2+fulv*state_ev['dd'][zzslices].reshape(szzzl)
            if pful: problem.parameters['rfu']['g'][wdls,...] = pfum*wdl*U*DU*UL + pfum*wudl*(DU**2-DUU)-dwdl/Re*DU
            if pfur: problem.parameters['rfu']['g'][wdrs,...] =                    pfum*wudr*U**2 -dwdr/Re*U
            if pwul: problem.parameters['fu']['g'][fuls,...]  = U*UL+DU*ful
            
            #ju.logger.info('fd    {:8.2e}'.format(dVV*ju.jsum(Iww[-1]*fd['g'])))
            #ju.logger.info('Iwux5 {:8.2e}'.format(dxx*wux5.sum()-1))
            #ju.logger.info('Iwux1 {:8.2e}'.format(dxx*wux1.sum()-1))
            #w['g'][1,...]=0
                
    if fbmax or fbmin or wbl is not None or wbr is not None: b.set_scales(AA)
    if fbmax:  b['g'] = np.maximum(bmaxf,b['g'])
    if fbmin:  b['g'] = np.minimum(bminf,b['g'])
 
    if fsmax or fsmin or wsl is not None or wsr is not None: s.set_scales(AA)
    if fsmax:  s['g'] = np.maximum(smaxf,s['g'])
    if fsmin:  s['g'] = np.minimum(sminf,s['g'])

    if fbmult:
        if wbr is not None: b['g'][wbrs,...] -= b['g'][wbrs,...] * wbr
        if wsr is not None: s['g'][wsrs,...] -= s['g'][wsrs,...] * wsr

    if 'pmss' in param: # moving source for the plume
        RRR=R
        if ju.mpirank==0:
            state_ev['SY']=ju.noiseint(state_ev['SY'],dt,param['pmsl'],param['pmss'])
            state_ev['SZ']=ju.noiseint(state_ev['SZ'],dt,param['pmsl'],param['pmss'])
            RRR=ju.get_radius(R,state_ev['SY'],state_ev['SZ'],Py,Pz)
        state_ev['SY'] = ju.comm.bcast(state_ev['SY'],0)
        state_ev['SZ'] = ju.comm.bcast(state_ev['SZ'],0)
        RRR            = ju.comm.bcast(RRR,0)
        avar['R']=RRR
        wrr=InletA*ju.pinlet(yy-state_ev['SY'],zz-state_ev['SZ'],Inlet,RRR,FFTAz)[0]
        if B is not None: fb['g'] = wrr*Bwxx
        #if S is not None: fs = fb

        if not dd and forcing<7: fu['g'] = wrr*wux4+mu*(1-wux4)*inwdr 
        #ju.logger.info('Inlet area {:7.4f},  R={:7.4f}, Y={:7.4f}, Z={:7.4f}'.format(dyy*dzz*ju.jsum(wrr)/IA,RRR,state_ev['SY'],state_ev['SZ']))
        #ju.srange('fb range [{:6.3f},{:6.3f}]',fb)
    if dxnoise is not None: dxnoise.int(dt)
    if dynoise is not None: dynoise.int(dt)
    if dznoise is not None: dznoise.int(dt)

    if sType=='pm':
        u.set_scales(AA)
        v.set_scales(AA)
        w.set_scales(AA)
        if ftopu    is not None: ftopu['g']  = dxx*np.sum(sdr*u['g'][sdrs,...]   ,axis=(0,)) # Top u
        if ftopv    is not None: ftopv['g']  = dxx*np.sum(sdr*v['g'][sdrs,...]   ,axis=(0,)) # Top u
        if ftopw    is not None: ftopw['g']  = dxx*np.sum(sdr*w['g'][sdrs,...]   ,axis=(0,)) # Top u
        if ftopuu   is not None: ftopuu['g'] = dxx*np.sum(sdr*u['g'][sdrs,...]**2,axis=(0,)) # Top uu
        if 'u'      in avar: gx.intyz(u['g'],dAAnwdr,avar['u'])
        if 'Iu'     in avar: avar['Iu'] = gx.integrate(avar['u'])
        if 'dudx'   in avar: gx.differentiate(avar['u'],out=avar['dudx'])
        if 'q'      in avar: gx.intyz2(v['g'],dnry,w['g'],dnrz,avar['q'])
        if 'topdiv' in avar: avar['topdiv']['g'] = ftopv.differentiate('y')['g']+ftopw.differentiate('z')['g']

        if topdivr: avar['topdivr']['g'] = nry**2*ftopv.differentiate('y')['g']+nrz**2*ftopw.differentiate('z')['g'] + (nry*ftopv['g']+nrz*ftopw['g'])/rrm

        
        if topu:  fd3 = (topT*state_ev['topd']['g']+dt*ftopu['g'])/(topT+dt)  # Set top divergence according to incoming velocity
        if topd:  fd3 = state_ev['topd']['g']+dt/topT*avar['topdiv']['g']     # Dynamic evolution of top divergence according to dv/dy+dw/dz

        # Use a gaussian for the top divergence accord to the radius defined by the velocity
        if topr: 
            Itu   = ju.comm.reduce(ftopu['g'].sum(),  op=MPI.SUM)  
            Ituu  = ju.comm.reduce(ftopuu['g'].sum(), op=MPI.SUM)  
            Itd   = ju.comm.reduce( state_ev['topd']['g'].sum(),     op=MPI.SUM)  
            Itdd  = ju.comm.reduce((state_ev['topd']['g']**2).sum(), op=MPI.SUM)  
            if ju.mpirank==0:
                avar['topIu']  = dAA*Itu
                avar['topIuu'] = dAA*Ituu
                avar['topId']  = dAA*Itd
                avar['topIdd'] = dAA*Itdd
                if avar['topIuu']>0: avar['topur']=avar['topIu']/np.sqrt(avar['topIuu']*CA)
                if avar['topIdd']>0: avar['topdr']=avar['topId']/np.sqrt(avar['topIdd']*CA)
                state_ev['topr'] = max(R,state_ev['topr']*(topT + dt*avar['topur']/avar['topdr'])/(topT + dt))
            r3=ju.comm.bcast(state_ev['topr'],0)
            fd3    = jtheta.periodic_gaussianu(np.sqrt(2)/r3,yy,zz,FW,FH)/(W*H)

        if topq:
            if ju.mpirank==0:   
                avar['topq']=dAA*(nry*ftopv['g']+nrz*ftopw['g']).sum()
                state_ev['topr'] *= np.exp(-dt/topT*avar['topq'])
                state_ev['topr'] = max(state_ev['topr'],R)
                ju.logger.info('topq {:8.3f}, topr: {:8.3f}'.format(avar['topq'],state_ev['topr']))
                r3=state_ev['topr']
            r3   = ju.comm.bcast(r3,0)
            fd3 = jtheta.periodic_gaussianu(np.sqrt(2)/r3,yy,zz,FW,FH)/(W*H)
  
        # Set the top divergence to match the incoming velocity topT can be zero
        if dd:
            fd.set_scales(AA)
            if ju.mpirank==0:
                if ddq:   state_ev['divx']['g'] *= np.exp(-dt/dT*avar['Iu']/(L*U))
                if ddq:   state_ev['divx']['g']  = (dT*state_ev['divx']['g']-dt*avar['q']['g'])/(dT+dt)
                else:       state_ev['divx']['g'] -= dt/dT*avar['dudx']['g']
                if divxI: state_ev['divx']['g'] -= dt/topT*avar['Iu']/(2*L)
                if wdivx is not None: state_ev['divx']['g']  *=wdivx
                if ddSc is not None: # Diffuse according to cubic rule
                    #state_ev['divx']['g'][ :2]=0 # divx should be even and zero at the ends
                    #state_ev['divx']['g'][-2:]=0
                    Ddx[1:-1] = np.diff(state_ev['divx']['g']) # Take derivative
                    Dc=Ddx**2
                    Dcmax=Dc.max()
                    if  Dcmax>0:
                        Dc[1:-1] = (Dc[:-2]+2*Dc[1:-1]+Dc[2:])/4
                        Dc=Dc*min(1/2,dt/(ddSc*Re*dxx**2))/Dc.max()
                        state_ev['divx']['g']  += np.diff(Dc*Ddx) # Perform nonlinear diffusion
                sevdivx = copy.copy(state_ev['divx']['g'])
            sevdivx = ju.comm.bcast(sevdivx,0)
            fd2 = sevdivx.reshape(szxx)*nwdr      # .reshape(szxx)
            Ifd2    = dVV*ju.jsum(fd2)
        else:
            fd2=0
            Ifd2=0

           
 
        if fd3 is not None:
            fd3  = np.maximum(0,fd3)
            if topdrt: fd3 = np.minimum(fd3,ju.jmax(fd3)*wdrn)
            Ifd1  = U*dAA*ju.jsum(wrr)
            Ifd3 = dAA*ju.jsum(fd3)
            if Ifd3<=0:
                fd3=np.ones((Nyyl,Nzzl))/(W*H)
            else:
                fd3*=(Ifd1+Ifd2)/Ifd3
            if 'topd' in state_ev: state_ev['topd']['g']=fd3

            fd3=fd3.reshape((1,Nyyl,Nzzl))
             
            fd['g']            = fd2
            fd['g'][wdls,...] += U*wdl*wrr
            fd['g'][wdrs,...] -= wdr*fd3
            if pwul: problem.parameters['fu']['g'][wdls,...] = wul*U*wrr
            if pwur: problem.parameters['fu']['g'][wdrs,...] = wur*fd3
            if pwul or pwur: problem.parameters['fu']['g'][wulrs,...] = U*wrr*wuli+fd3*wuri
            #if pwul: problem.parameters['rfu']['g'][wdls,...] = (U*wrr)**2*ful
            #if pwur: problem.parameters['rfu']['g'][wdrs,...] = (fd3.reshape((1,Nyyl,Nzzl)))**2*fur
            if (pwul or pwur) and  np.mod(solver.iteration,100)==0: JA.savexyz('force/fu-{:05d}.hdf5'.format(np.int(solver.iteration/100)),('fu'),(problem.parameters['fu']['g']),AA)
        
            Ifd3    = dAA*ju.jsum(fd3)
            Ifd     = dVV*ju.jsum(fd['g']) 
            #ju.logger.info('Ifd1: {:8.4f}, Ifd2: {:8.4f}, Ifd3: {:8.4f}, Ifd: {:9.2e}'.format(Ifd1/IAU,Ifd2/IAU,Ifd3/IAU,Ifd/IAU))
            #if np.mod(solver.iteration,10)==0: 
            
            #ju.logger.info('r3: {:8.4f}, IIu: {:8.4f}, Itu: {:8.4f}, Ituu: {:8.4f}, Ifd1: {:8.4f}, Ifd2: {:8.4f}, Ifd3: {:8.4f}, Ifd: {:9.2e}'.format(r3,IIu,Itu/IAU,Ituu/IAUU,Ifd1/IAU,Ifd2/IAU,Ifd3/IAU,Ifd/IAU))
            #if wwl is not None:  w['g'][wwls,...] -= dt*wwl*w['g'][wwls,...]
            #if wwr is not None:  w['g'][wwrs,...] -= dt*wwr*w['g'][wwrs,...]
            #if wwl is not None:  v['g'][wwls,...] -= dt*wwl*v['g'][wwls,...]
            #if wwr is not None:  v['g'][wwrs,...] -= dt*wwr*v['g'][wwrs,...]
            if dtjdivyz is not None:
                if solver.sim_time >= dtjdivyz*dtjdivyzn:
                    dtjdivyzn +=1
                    JA.saveyz(ddir / 'divyz/divyz-{:05d}.hdf5'.format(dtjdivyzn), 'divyz', fd3, AA)
        else: fd['g'][wdls,...] = U*wdl*wrr
    if dd and sType=='pm' and forcing==6:
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
            state_ev['qvx'] = state_ev['qvx']-dt/dT*qv1
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
        fd['g'] = fd['g'] * wux1  + np.reshape(state_ev['qvx'],szxx)*nwdr
        Iqv=dVV*ju.jsum(fd['g'])  # Now add negative divergence to make the system divergence free
        fd['g']=fd['g']-Iqv*Iux1  # This correction should be small becuase of mu imposition
        
        #ju.logger.info('uE {:9.2e}, dnoise {:8.4f}, NYZ {:8.4f}'.format(uE,dnoise,np.sqrt(NY**2+NZ**2)))
        #ju.logger.info('uE {:9.2e}, EQ {:8.4f}, s+b {:8.4f}, Fx {:8.4f}, Ifd {:9.2e}, NN {:8.4f}'.format(uE,EQ/(U*IA*H*W),Isb,state_ev['Fx'],ju.jint(fd),NN))


        #JA.savexyz('fu/fu-{:05d}.hdf5'.format(solver.iteration),'u',fu['g'],AA)
        #JA.savexyz('fd/fd-{:05d}.hdf5'.format(solver.iteration),'d',fd['g'],AA)


                
    if tX1 is not None: tX0=copy.copy(tX1)
    if sType=='gc' and np.isfinite(X):
        XO=copy.copy(X)
        if Xmax is not None:
            if X>Xmax: ju.set_abort('X exceeded Xmax {:5.2f} {:5.2f}'.format(X,Xmax),True)
        if Xmin is not None:
            if X<Xmin: ju.set_abort('Xmin exceeded X {:5.2f} {:5.2f}'.format(X,Xmin),True)
                
            
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
        
    if FB    and b is not None: b['g'] = (FB/max(2/FB,ju.jint(b)))*b['g']
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
    
    if forcing==8 and sType=='gc' and np.isfinite(X):
        problem.parameters['wu']['g']=np.sin(math.pi/L*(xx-X))**math.floor(Nx/2)
    
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
        ss2=ju.write_stats(problem,solver,domain,dt,X,U,g,state_ev,gIxx,gIyy,gIzz,tops,rmstopd)
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
    if iteration_one:
        for a in ('rfu','rfv','rfw','fu','fv','fw'):
            if a in problem.parameters: JA.savexyz('force/'+a+'.hdf5',a, problem.parameters[a]['g'],AA)
        for a in avar:     ju.print_v(a,avar[a],    'First iteration: avar ')
        for a in state_ev: ju.print_v(a,state_ev[a],'First iteration: sev  ')
        iteration_one=False
    #ju.srange('IT1 state_ev[dd] [{:6.3f},{:6.3f}]',state_ev['dd'])

    if debug_fu:
        JA.savexyz('wdu/wdu-{:05d}.hdf5'.format(solver.iteration),('fu','wu','fd'),(problem.parameters['fu']['g'],problem.parameters['wu']['g'],problem.parameters['fd']['g']),AA)

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
    ju.logger.info('Updated jrec')
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
    if PIDG: param['Gravity']  = state_ev['g']
    else:    param['Velocity'] = state_ev['U']
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
    
                #if pmIu: state_ev['div']x = np.maximum(0,wdivx*(state_ev['divx']-wux4*(dt/dT*Idudx+Anwdr*avar['Iu']*10*dt)))
                #else:    state_ev['divx'] = np.maximum(0,wdivx*(state_ev['divx']-dt/dT*Idudx))
                #        state_ev['divx'] = np.maximum(0,wux6*(state_ev['divx']-dt/dT*wux4*Idudx))
                #if fixI: state_ev['divx'] = pmif*(state_ev['divx']*pmif).sum()*dxx
                #state_ev['divx'][1:-1] += (state_ev['divx'][:-2]-2*state_ev['divx'][1:-1]+state_ev['divx'][2:])/4
                # Time advance divx
                #Idudx[1:-1] = (Idudx[:-2]+2*Idudx[1:-1]+Idudx[2:])/4
                #Idudx[-2:]  = (Idudx[-2] +  Idudx[-1])/2          
                #Idudx[ :2]  = (Idudx[ 0] +  Idudx[ 1])/2          
                #state_ev['divx'] = np.maximum(0,(dT*state_ev['divx']-dt*wdivx*Idudx)/(dT+dt))
