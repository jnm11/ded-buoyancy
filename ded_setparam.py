#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ded_settime.py param.h5 300
ded_h5print.py param.h5
"""

import sys
import h5py
import numpy as np

n=len(sys.argv) - 1
name  = sys.argv[1]
param = sys.argv[2]
value = sys.argv[3]
f = h5py.File(name, 'r+')

if len(sys.argv)>3 :
    type=sys.argv[4]
    #print('Forcing type {} for {}'.format(type,param))
else:
    type=np.dtype(f[param])
    print('Exists type {} for {}'.format(type,param))
    
if param in f:  del f[param]
    
if type=="int64":
    f.create_dataset(param,data=np.int(value))
elif type=="float64":
    f.create_dataset(param,data=np.float64(value))
elif type=="object" or type=="str":
    f.create_dataset(param,data=np.str(value))
elif type=="bool":
    if   (value=='T' or value=='True'):  f.create_dataset(param,data=bool(True))
    elif (value=='F' or value=='False'): f.create_dataset(param,data=bool(False))
    else:  print('Unknown boolean value {}'.format(value))
elif type=="delete":
    print('Deleting {}'.format(param))
else:
    print('Unknown type {}'.format(type))
    
if param in f: print("Setting {} to {} in {} {}".format(param,f[param],name,type))
f.close();
    







