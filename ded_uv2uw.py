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


def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.iteritems():
        print("    {}: {}".format(key, val))
            
                
def foo(name, obj):
   print(name, obj)
   return None

n=len(sys.argv) - 1
k=0
while (k<n):
    k=k+1
    name =  sys.argv[k]
    f = h5py.File(name, 'r+')
    #print("File {} opened".format(name))

    if "tasks/av" in f:
        print("File {} contains av ".format(name))
        del f["tasks/av"]

    if "tasks/auv" in f:
        print("File {} contains auv ".format(name))
        del f["tasks/auv"]

    a = ['tasks/ab','tasks/ap','tasks/au','tasks/aw','tasks/auu', 'tasks/avv','tasks/aww', 'tasks/auw']
    for b in a:       
        if not b in f:
            print("File {} does not contains {}".format(name,b))
            f[b] = f["tasks/u"]
            data=f[b]
            data[...]=0*data[...]
    f.close()        

    






