#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#
#
# ded_extend.py pm/150/y/y-00000.hdf5 pm/150/new.hdf5
"""
Extend fluid simulation state

First the input state is doubled along requested directions with the correct parity

Usage: 
   ded_gc.py [options] NAMEIN NAMEOUT

Options:
  -H, --Height=<>   : New Box height       
  -W, --Width=<>    : New Box width        
  -L, --Length=<>   : New Box length       
  --Nx=<>           : New Nx
  --Ny=<>           : New Ny
  --Nz=<>           : New Nz
  --scalex          : Change x scale
  --scaley          : Change y scale
  --scalez          : Change x scale
  --doublex         : Double in x direction
  --doubley         : Double in y direction
  --doublez         : Double in z direction
  --parityx=<>      : Change x parity
  --parityy=<>      : Change y parity
  --parityz=<>      : Change x parity
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore',category=FutureWarning) 
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
from docopt import docopt
from mpi4py import MPI
import os
import sys
import math
import h5dict
import jutil as ju
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools  import post
from polytrope.tools.checkpointing import Checkpoint
from pathlib import Path

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

  
namein  =  np.str(args['NAMEIN'])
nameout =  np.str(args['NAMEOUT'])

ju.logger.info('Changing {} to {}'.format(namein,nameout))

# Read the output directory from the environment variable DEDALUS_DATA
bdir = Path(os.environ['DEDALUS_DATA']) / "gc"
ddir = bdir / name

ju.set_base_dir(ddir)
ju.check_version(ver)

pfn=Path('final')
fpfn=ju.find_param_file(bdir,name,pfn)
param=ju.read_param(fpfn)
Nxa       = param['Nx']    
Nya       = param['Ny']    
Nza       = param['Nz']    
Ha        = param['H']     
Wa        = param['W']     
La        = param['L']     
average   = param['average']     

Hb  = np.float(args['--Height'])
if args['--Width']:     Wb  = np.float(args['--Width'])
else:                   Wb  = Wa
if args['--Length']:    Lb  = np.float(args['--Length'])
else:                   Lb  = La
if args['--Nx']:        Nxb = np.int(args['--Nx'])
else:                   
if args['--Ny']:        Nyb = np.int(args['--Ny'])
if args['--Nz']:        Nyb = np.int(args['--Nz'])
else:                   Nyb = np.round(Hb/Ha*Nza)
scalex   = np.bool(args['--scalex'])
scaley   = np.bool(args['--scaley'])
scalez   = np.bool(args['--scalez'])
doublex  = np.bool(args['--doublex'])
doubley  = np.bool(args['--doubley'])
doublez  = np.bool(args['--doublez'])
paritybx = np.int(args['--parityx'])
parityby = np.int(args['--parityy'])
paritybz = np.int(args['--parityz'])


if scalex is None: scalex=False
if scaley is None: scaley=False
if scalez is None: scalez=False

if doublex is None: doublex=False
if doubley is None: doubley=False
if doublez is None: doublez=False

if parityx is None: parityx=False
if parityy is None: parityy=False
if parityz is None: parityz=False

da, pa=ju.boussinseq_init(Nxa,Nya,Nza,La,Wa,Ha,average,'b')

Na=(pa.Nx,pa.Ny,pa.Nz)
Sa=(pa.L,pa.W,pa.H)
vv=sa.state.keys()
fa=sa.state[vv[0]]
paritya = np.array((dima,),dtype=np.int64)
Ea      = np.array((dima,),dtype=np.float64)
Sa      = np.array((dima,),dtype=np.float64)
for j in range(dima):
    cna[j]     = da.bases[j].name
    paritya[j] = fa.meta[cna[j]['parity']
    Sa[j]=da.bases[j].interval[0]       
    Ea[j]=da.bases[j].interval[1]
                         



if dimb==3: parityb=[parityx,parityy,parityz];
else:       parityb=[parityx,parityz];



if dima==dimb: cnb=cna
else:          cnb=[cna[0],'y',cna[1]]



if Nxb is None: Nxb = np.round(Lb/La*Nxa)
if Nyb is None: Nyb = np.round(Wb/Wa*Nya)
if Nzb is None: Nzb = np.round(Hb/Ha*Nza)

if Lxb is None and     doublex: Lxb = 2*Lxa
if Lxb is None and not doublex: Lxb =   Lxa
if Wyb is None and     doubley: Wyb = 2*Wya
if Wyb is None and not doubley: Wyb =   Wya
if Hzb is None and     doublez: Hzb = 2*Hza
if Hzb is None and not doublez: Hzb =   Hza

if Nxb is None and     doublex: Nxb = 2*Nxa
if Nxb is None and not doublex: Nxb =   Nxa
if Nyb is None and     doubley: Nyb = 2*Nya
if Nyb is None and not doubley: Nyb =   Nya
if Nzb is None and     doublez: Nzb = 2*Nza
if Nzb is None and not doublez: Nzb =   Nza




db, pb=ju.boussinseq_init(Nxb,Nyb,Nzb,Lb,Wb,Hb,average,'b')


sa=jrec.jsolve(da.variables.keys())
sb=jrec.jsolve(da.variables.keys())

fnrst=ju.find_start_file(bdir,namein,Null)  

ja=jrec.jrec(da,ddir,pa,False,None)
jb=jrec.jrec(db,ddir,pb,True,None)
dt,ja.state_ev=ja.load_state(fnrst,fn,sa,sa.state_ev)

if Nya>1: dima=3
else:     dima=2
if Nyb>1: dimb=3
else:     dimb=2




c                         
                         
for n in cna:
if parityx is None: sa.state[vv[0]].meta['x'['parity']
if parityy is None: sa.state[vv[0]].meta['x'['parity']
if parityz is None: sa.state[vv[0]].meta['x'['parity']


for v in vv:
    fa=sa.state[v]
    fb=sb.state[v]
    for j in range(dima):
        fb[cn]['scale']=1
        if double[j]:
            Na[j]=Na[j]*2
            if   fa[cn]['parity'] ==  1: np.concatenate(np.flip(+y,afais=j),y),afais=j,out=y)        
            elif fa[cn]['parity'] == -1: np.concatenate(np.flip(-y,afais=j),y),afais=j,out=y)        
            else:                       np.concatenate(                 y,y),afais=j,out=y)        
        if interp[j] and Na[j] != Nb[j] fa[cn]['scale']=Nb[j]/Na[j]
        else:                           fa[cn]['scale']=1
 

if dimb==3: y[0:Nxb,0:Nyb,0:Nzb]=x
else:       y[0:Nxb,0:Nzb]=x

for j in range(dim):





#dt = checkpoint.restart(fnrst, sa)
ba = sa.state['b']
ua = sa.state['u']
va = sa.state['v']
wa = sa.state['w']
ba.set_scales(1)
ua.set_scales(1)
va.set_scales(1)
wa.set_scales(1)

bb = sb.state['b']
ub = sb.state['u']
vb = sb.state['v']
wb = sb.state['w']
bb.set_scales(1)
ub.set_scales(1)
vb.set_scales(1)
wb.set_scales(1)


save_state(self,dt,state_ev):



