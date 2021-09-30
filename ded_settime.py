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
name = sys.argv[1]
time = np.float(sys.argv[2])
print("Setting time to ",str(time)," in ",name)
f = h5py.File(name, 'r+')
T = f['T']
T[...]=time
f.close();

    






