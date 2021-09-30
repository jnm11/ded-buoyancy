#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
remove uv entries in dump file and replace with uw

Usage: 
   ded_gc.py [options] NAME

Options:
  -h, --help            : show this help message

"""

import sys
import h5py
import numpy as np
import pprint
import h5dict

def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.items():
        print("    {}: {}".format(key, val))
            
                
def foo(name, obj):
   print(name, obj)
   return None

n=len(sys.argv) - 1
k=0

def print_type(x,s):
    if isinstance(x, dict):
        for a in x: print_type(x[a],s+'/'+a)
    elif   isinstance(x, np.ndarray):
        if np.size(x)<20:
            if   x.dtype=='float64':    print("{:25s}: {:11s} {} {:>11.4f}".format(s, 'afloat64',   x.shape, x))
            elif x.dtype=='int64':      print("{:25s}: {:11s} {} {}".format(s, 'aint64',     x.shape, x))
            elif x.dtype=='complex128': print("{:25s}: {:11s} {} {}".format(s, 'complex128', x.shape, x))
            else:                       print("{:25s}: {:11s} {} {}".format(s, str(x.dtype), x.shape, x))
        else:
            if   x.dtype=='float64':    print("{:25s}: {:11s} {} [{:>11.5f} {:>11.5f}]".format(s, 'afloat64',   x.shape, x.min(), x.max()))
            elif x.dtype=='int64':      print("{:25s}: {:11s} {} [{:>11d} {:>11d}]".format(    s, 'aint64',     x.shape, x.min(), x.max()))
            elif x.dtype=='complex128': print("{:25s}: {:11s} {} [{:>11d} {:>11d}]".format(    s, 'complex128', x.shape, x.min(), x.max()))
            else:                       print("{:25s}: {:11s} {} [{} {}]".format(            s, str(x.dtype), x.shape, x[0],x[-1]))
    elif isinstance(x, np.float):   print("{:25s}: {:11s} {:>11.5f}".format(              s, 'np.float', x))
    elif isinstance(x, np.int):     print("{:25s}: {:11s} {:>11d}".format(                s, 'np.int', x))
    elif isinstance(x, np.float64): print("{:25s}: {:11s} {:11.5f}".format(               s, 'np.float64', x))
    elif isinstance(x, np.int64):   print("{:25s}: {:11s} {:>11d}".format(                s, 'np.int64', x))
    elif isinstance(x, np.bool_):   print("{:25s}: {:11s} {:>11s}".format(                s, 'np.bool_', str(x)))
    elif isinstance(x, np.str_):    print("{:25s}: {:11s} {:>11s}".format(                s, 'np.str_', x))
    elif isinstance(x, bool):       print("{:25s}: {:11s} {:>11s}".format(                s, 'bool', str(x)))
    elif isinstance(x, float):      print("{:25s}: {:11s} {:>11.5f}".format(              s, 'float', x))
    elif isinstance(x, int):        print("{:25s}: {:11s} {:>11d}".format(                s, 'int', x))
    elif isinstance(x, str):        print("{:25s}: {:11s} {:>11s}".format(                s, 'str', x))
    else:                           print("{:25s}: {:11s} {} else".format(               s, str(type(x)), x))


while (k<n):
    k=k+1
    name =  sys.argv[k]
    print("")
    print("Filename   : {}".format(name))
    #    f = h5py.File(name, 'r')
    p = h5dict.load_dict_from_hdf5(name)
    print_type(p,'  ')

    #    f.visititems(foo)
    #a=list(f.keys())
    #a = ['tasks/ab','tasks/ap','tasks/au','tasks/aw','tasks/auu', 'tasks/avv','tasks/aww', 'tasks/auw']
    #for b in a:       
    #    pprint(b,":",f[b])
    #
    #    f.close()        

    






