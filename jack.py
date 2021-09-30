#!/usr/bin/env python3
#import matplotlib
#matplotlib.use("TkAgg")
import dedalus.public as de
import matplotlib.pyplot as plt
import numpy as np
import sys
# Resolution 1
x_basis = de.Fourier('x', 4, interval = (0,1))
domain = de.Domain([x_basis], np.complex128)

plotting_scale=2
mode = 2*np.pi
x = domain.grid(0, scales=plotting_scale)
kx = domain.elements(0)
a = domain.new_field()
a['c'][kx == mode] = 1
print(" a: {}, {}".format(a['c'].shape, a['c']))
print("kx: {}, {}".format(kx.shape, kx))
a.set_scales(plotting_scale)
plt.figure(1)
plt.plot(x, a['g'].real, '-o')
plt.plot(x, a['g'].imag, '--')
# plt.show()

# Resolution 2
x_basis2 = de.Fourier('x2', 16, interval = (0,1))
domain2 = de.Domain([x_basis2], np.complex128)

# mode = 2
x2 = domain2.grid(0)
kx2 = domain2.elements(0)
a2 = domain2.new_field()
a2['c'][kx2 == mode] = 1
print(" a2: {}, {}".format(a2['c'].shape, a2['c']))
print("kx2: {}, {}".format(kx2.shape, kx2))

# plt.figure(2)
plt.plot(x2, a2['g'].real, '-o')
plt.plot(x2, a2['g'].imag, '--')


rand = np.random.RandomState()#(seed=42)
N=128
L=1
x_basis3 = de.Fourier('x3', N, interval = (0,L))
domain3  = de.Domain([x_basis3], np.float64)
x3       = domain3.grid(0)
kx3      = domain3.elements(0)/(2*np.pi*(N-1))
n1       = domain3.new_field()
n2       = domain3.new_field()

n1['g']=rand.standard_normal(x3.shape)
n1['c'][kx3 == 0]=0
n2['c']=n1['c']*np.exp(-(20*kx3)**2)
print(kx3)
plt.figure(2)
plt.plot(x3, n1['g'],x3, n2['g'])





plt.show()


