#!/usr/bin/env python3
#import matplotlib
#matplotlib.use("TkAgg")
from mpi4py import MPI
import dedalus.public as de
import matplotlib.pyplot as plt
import numpy as np
import sys

comm      = MPI.COMM_WORLD
mpirank   = comm.Get_rank()
mpisize   = comm.Get_size()

# Resolution 1
m1=4
m2=2
AA=3/2
Nx=64
Ny=30
Nz=48

x_basis = de.Fourier(  'x', Nx, interval = (0,4), dealias=AA)
y_basis = de.Fourier(  'y', Ny, interval = (0,1), dealias=AA)
z_basis = de.Chebyshev('z', Nz, interval = (0,1), dealias=AA)
domain = de.Domain([x_basis,y_basis,z_basis], np.float64, mesh=(m1,m2))
if mpirank==0:
    print('  scales         {}'.format(1))
    print('    global_shape {}'.format(domain.dist.grid_layout.global_shape(scales=1))) 
print('    local_shape  {}'.format(domain.dist.grid_layout.local_shape(scales=1))) 
print('    slices       {}'.format(domain.dist.grid_layout.slices(scales=1)))         
print('    start        {}'.format(domain.dist.grid_layout.start(scales=1)))           
if mpirank==0:
    print('                   ')
    print('  scales         {}'.format(AA))
    print('    global_shape {}'.format(domain.dist.grid_layout.global_shape(scales=AA))) 
print('    local_shape  {}'.format(domain.dist.grid_layout.local_shape(scales=AA))) 
print('    slices       {}'.format(domain.dist.grid_layout.slices(scales=AA)))         
print('    start        {}'.format(domain.dist.grid_layout.start(scales=AA)))           
print('                   ')

quit()

if mpirank==0:
    print('domain.dist.rank:     {}'.format(domain.dist.rank))   
    print('domain.dist.size:     {}'.format(domain.dist.size))
    print('domain.dist.coords:   {}'.format(domain.dist.coords)) 
    print('domain.dist.mesh:     {}'.format(domain.dist.mesh))  
    print('domain.grid(0).shape: {}'.format(domain.grid(0).shape))
    print('domain.grid(1).shape: {}'.format(domain.grid(1).shape))
    print('domain.grid(2).shape: {}'.format(domain.grid(2).shape))
    print('                   ')
    print('domain.dist.grid_layout')
    print('  ext_coords     {}'.format(domain.dist.grid_layout.ext_coords))  
    print('  ext_mesh       {}'.format(domain.dist.grid_layout.ext_mesh))    
    print('  index          {}'.format(domain.dist.grid_layout.index))       
    print('  local          {}'.format(domain.dist.grid_layout.local))       
    print('  grid_space     {}'.format(domain.dist.grid_layout.grid_space))  
    print('                   ')
    print('                   ')
    print('                   ')
    print("vars(domain):      {}".format(vars(domain)))
    print("vars(distributor): {}".format(vars(domain.distributor)))
    print("vars(dist):        {}".format(vars(domain.dist)))
    print('vars(domain.dist.comm):             {}'.format(dir(domain.dist.comm)))
    print('vars(domain.dist.comm_cart):        {}'.format(dir(domain.dist.comm_cart)))
    print('vars(domain.dist.layouts:           {}'.format(dir(domain.dist.layouts)))
    print('vars(domain.dist.paths:             {}'.format(dir(domain.dist.paths)))
    print('vars(domain.dist.coeff_layout:      {}'.format(dir(domain.dist.coeff_layout)))
    print('vars(domain.dist.grid_layout:       {}'.format(dir(domain.dist.grid_layout)))
    print('vars(domain.dist.grid_layout.blocks:      {}'.format(dir(domain.dist.grid_layout.blocks)))
    print('vars(domain.dist.grid_layout.buffer_size: {}'.format(dir(domain.dist.grid_layout.buffer_size)))
    print('vars(domain.dist.grid_layout.domain:      {}'.format(dir(domain.dist.grid_layout.domain)))
    print('vars(domain.dist.grid_layout.ext_coords.shape:  {}'.format(dir(domain.dist.grid_layout.ext_coords.shape)))
    print('vars(domain.dist.grid_layout.ext_mesh.shape:    {}'.format(dir(domain.dist.grid_layout.ext_mesh.shape)))
    print('vars(domain.dist.grid_layout.global_shape:{}'.format(dir(domain.dist.grid_layout.global_shape)))
    print('vars(domain.dist.grid_layout.grid_space.shape:  {}'.format(dir(domain.dist.grid_layout.grid_space.shape)))
    print('vars(domain.dist.grid_layout.local_shape: {}'.format(dir(domain.dist.grid_layout.local_shape)))
    print('vars(domain.dist.grid_layout.slices:      {}'.format(dir(domain.dist.grid_layout.slices)))
    print('vars(domain.dist.grid_layout.start:       {}'.format(dir(domain.dist.grid_layout.start)))
    print('vars(domain.dist.comm):             {}'.format(vars(domain.dist.comm)))
    print('vars(domain.dist.comm_cart):        {}'.format(vars(domain.dist.comm_cart)))
    print('vars(domain.dist.layouts:           {}'.format(vars(domain.dist.layouts)))
    print('vars(domain.dist.paths:             {}'.format(vars(domain.dist.paths)))
    print('vars(domain.dist.coeff_layout:      {}'.format(vars(domain.dist.coeff_layout)))
    print('vars(domain.dist.grid_layout:       {}'.format(vars(domain.dist.grid_layout)))
    print('vars(domain.dist.layout_references: {}'.format(vars(domain.dist.layout_reference)))
    
