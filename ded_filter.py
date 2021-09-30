#!/usr/bin/env python3
import numpy as np
import math
import matplotlib.pyplot as plt
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools  import post
from dedalus.core   import operators
plt.ion()
rand = np.random.RandomState()
L=4
Lc=.5;
Nx=128;
x_basis=de.Fourier(  'x', Nx, interval=(0, L), dealias=3/2)
domain=de.Domain([x_basis], grid_dtype=np.float64)

kc=np.float(L/Lc);

k = domain.bases[0].wavenumbers[:]/kc

#print(L*k/(Nx*math.pi))
nfx = domain.new_field(name='nfx')
nfx.set_scales(1)
#print(nfx['c'].shape)
gshape = domain.dist.grid_layout.global_shape(scales=1/2)
slices = domain.dist.grid_layout.slices(scales=1/2)
nfx['c'] = (rand.standard_normal(gshape)[slices]+1J*rand.standard_normal(gshape)[slices])*np.exp(-k**2)
nfx['c'][0]=0
nfx['g']
nfx.set_scales(3/2, keep_data=True)
print(np.mean(nfx['g']))
ax = plt.gca()
ax.plot(nfx['g'])

plt.draw()
plt.waitforbuttonpress()

#print(nfx['c'].shape)
#print(nfx.data)
