#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Load a Dedalus state and write it back out at a different resolution

Usage: 
   load_state.py [options] FNIN FNOUT 

Options:
      --Nx=<>           : x grid points    
      --Ny=<>           : y grid points   
      --Nz=<>           : z grid points   

"""

from docopt import docopt
import copy
import h5py 
import jutil as ju
import pathlib
import logging
import numpy as np
from dedalus import public as de
from pathlib import Path

ju.logger=logging.getLogger(__name__)

if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

    
def make_basis(N,I,c,T):
    if   T=='Chebyshev': b=de.Chebyshev(c, N, interval=I, dealias=1)
    elif T=='SinCos':    b=de.SinCos(   c, N, interval=I, dealias=1)
    elif T=='Fourier':   b=de.Fourier(  c, N, interval=I, dealias=1)
    else:
        ju.logger.error('No match for {} {}'.format(T,len(T)))
        quit()
    return(b)

def load_state(fn,fnout,args):
    s=dict()
    t=dict()
    g=h5py.File(str(fnout), 'w')
    try:
        f=h5py.File(str(fn), 'r')
        if 'gsize'        in f['scales']: s['gsize']        = f['scales']['gsize'][:]
        if 'interval'     in f['scales']: s['interval']     = f['scales']['interval'][:]
        if 'iteration'    in f['scales']: s['iteration']    = f['scales']['iteration'][()]
        if 'sim_time'     in f['scales']: s['sim_time']     = f['scales']['sim_time'][()]
        if 'wall_time'    in f['scales']: s['wall_time']    = f['scales']['wall_time'][()]
        if 'write_number' in f['scales']: s['write_number'] = f['scales']['write_number'][()]
        if 'timestep'     in f['scales']: s['dt']           = f['scales']['timestep'][()]
        if 'type'         in f['scales']: s['type']         = f['scales']['type']                
        #if 'x'            in f['scales']: s['x']            = f['scales']['x'][:]
        #if 'y'            in f['scales']: s['y']            = f['scales']['y'][:]
        #if 'z'            in f['scales']: s['z']            = f['scales']['z'][:]
        #for a in s.keys(): ju.logger.info('load_state: /scales/{:14s}: {}'.format(a,s[a]))
    except:
        ju.logger.info('load_state: scales failed {}'.format(str(fn)))
        quit()

    g.create_group('scales')#Coordinates will be wrong as well as size
    g.create_group('tasks')#Coordinates will be wrong as well as size

    typs=['interval','iteration','sim_time','wall_time','write_number','timestep','type']             
    for a in typs:g['scales'].create_dataset(a,data=f['scales'][a])
    
    dim=s['gsize'].size
    #ju.logger.info('load_state: dimension: {}'.format(dim))

    if ju.mpirank==0: comm0 = ju.comm.Split(0,0)
    else:             comm0 = ju.comm.Split(1,mpirank)


    if dim==2: c=('x','z')
    else:      c=('x','y','z')
    N1=s['gsize']
    N2=copy.copy(N1)
    if args['--Nx']: N2[ 0]=args['--Nx']
    if args['--Ny']: N2[ 1]=args['--Ny']
    if args['--Nz']: N2[-1]=args['--Nz']

    T=list()
    for j in range(dim): T.append(s['type'][j].decode("ASCII", "ignore").strip())
    ju.logger.info('load_state: types {}'.format(T))

    
    b1=list()
    for j in range(dim): b1.append(make_basis(N1[j],s['interval'][j],c[j],T[j]))
    d1=de.Domain(b1, grid_dtype=np.float64,comm=comm0)

    b2=list()
    for j in range(dim):
        b2.append(make_basis(N2[j],s['interval'][j],c[j],T[j]))
        g['tasks'][c[j]]=b2[j].grid(1)
    g['scales']['gsize']=N2  
    d2=de.Domain(b2, grid_dtype=np.float64,comm=comm0)

    for j in range(dim):
        if T[j]=='Fourier':
            N2[j] = N2[j]/2
            N1[j] = N1[j]/2
            
    sl=list()
    for j in range(dim): sl.append(slice(0,min(N1[j],N2[j]),None))

    print('N1:  {}'.format(N1))
    print('N2:  {}'.format(N2))
    print('sl:  {}'.format(sl))
    sl=tuple(sl)
    ju.logger.info('load_state: input  size: {}'.format(N1))
    ju.logger.info('load_state: output size: {}'.format(N2))
    
    x=d1.new_field()
    y=d2.new_field()

    for a in f['tasks'].keys():
        p=f['tasks'][a].attrs['parity']
        for k in range(dim):
            if p[k]!=0:
                x.meta[c[k]]['parity']=p[k]
                y.meta[c[k]]['parity']=p[k]
        x['g']=f['tasks'][a]
        y['c'][sl]=x['c'][sl]
        g['tasks'][a]=y['g']
        g['tasks'][a].attrs['parity']=p
        
    f.close()
    g.close()


fnin=Path(args['FNIN'])
fnout=Path(args['FNOUT'])

load_state(fnin,fnout,args)


