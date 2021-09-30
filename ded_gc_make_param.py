#!/usr/bin/env python3
import h5py
import numpy as np
import pprint
f = h5py.File("/tmp/gc-param.hdf5", "w")
bc = f.create_group("bc")  # Boundary conditions
ic = f.create_group("ic")  # initial  conditions
pp = f.create_group("pp")  # program  parameters
op = f.create_group("op")  # output   parameters
fp = f.create_group("fp")  # forcing  parameters

dt = h5py.special_dtype(vlen=unicode)

bc.create_dataset('lb',data=np.str("noslip"))
bc.create_dataset('ub',data=np.str("slip"))


#group = f.get("bc")
#group.items()
#

f.close()

f = h5py.File("/tmp/gc-param.hdf5", "r")
for name in f:
    print(name)


    
fbc=f['bc']
for name in fbc:
    print(name)

print(fbc['lb'])
print(fbc['ub'])
print(fbc['x'])

x=f['bc']['x'][...]
lb=f['bc']['lb'][...]
if lb=="noslip":
    print("noslip")
else:
    print("slip")





