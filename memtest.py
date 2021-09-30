#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test memory useage

"""
Test memory useage

VMEM size of mapped memory which the process occupies
DATA stack and heap limit doesn't work on OSX

Usage: 
   memtest.py [options]

Options:
  -C, --collect=<>  : delete variable each iteration       [default: 0]
  -D, --delete=<>   : garbage collect each iteration       [default: 0]
  -S, --steps=<>    : fractional step in memory power of 2 [default: 1]
  -M, --start=<>    : starting           memory power of 2 [default: 10]
  -N, --end=<>      : ending             memory power of 2 [default: 20]
  -R, --RSS=<>      : set maximum DATA (=0 unlimited)      [default: 0]
  -V, --VMEM=<>     : set maximum VMEM (=0 unlimited)      [default: 0]
"""

# on asahi
# with varaiable deletion works DATA segment size control
# memtest.py -M 22 -N 24 -S 4 -R 0.7 -D 1
# without varaiable deletion fails
# memtest.py -M 22 -N 24 -S 4 -R 0.7 -D 0
#
# with varaiable deletion works VM size control
# memtest.py -M 22 -N 24 -S 4 -D 1 -V 1.14
# memtest.py -M 22 -N 24 -S 4 -D 0 -V 1.14
#
# on tokachi
# DATA segment  and VM size control don't do anything
#
# with varaiable deletion works VM size control
# memtest.py -M 22 -N 24 -S 4 -D 1 -V 1.14
# memtest.py -M 22 -N 24 -S 4 -D 0 -V 1.14
#
#
#memtest.py -M 24 -N 26 -S 8 -D 0
#memtest.py -M 24 -N 26 -S 8 -D 1
#memtest.py -M 10 -N 30 -S 1 -D 0
#memtest.py -M 10 -N 30 -S 1 -D 1 -R 0.2 -V 6

import sys
import os
import psutil
import resource
import numpy as np
import gc
from mpi4py import MPI
from docopt import docopt
if __name__ == '__main__':
    args = docopt(__doc__)  # parse arguments based on docstring above

#print(args)

OSX   = (sys.platform == "darwin")
LINUX = (sys.platform == "linux")

if LINUX: MS=1024
if OSX:   MS=1

D = np.bool(int(args['--delete']))
C = np.bool(int(args['--collect']))
S = np.int(args['--steps'])
M = np.int(args['--start'])
N = np.int(args['--end'])
R = np.int(np.float(args['--RSS'])*1024**3)
V = np.int(np.float(args['--VMEM'])*1024**3)

comm      = MPI.COMM_WORLD
mpirank   = comm.Get_rank()
mpisize   = comm.Get_size()

rand = np.random.RandomState()#(seed=42)

pid=os.getpid()
process = psutil.Process(pid)

print('OSX   {}'.format(OSX))
print('LINUX {}'.format(LINUX))

if R>0: print('DATA  {:6.3f} GB'.format(float(R)/1024**3))
if V>0: print('VMEM  {:6.3f} GB'.format(float(V)/1024**3))

#https://docs.python.org/3/library/resource.html
if R>0: resource.setrlimit(resource.RLIMIT_DATA,(R, resource.RLIM_INFINITY)) # in bytes
if V>0: resource.setrlimit(resource.RLIMIT_AS,  (V, resource.RLIM_INFINITY)) # in bytes
#RLIMIT_RSS  This has no effect

def get_mss():
    m=MS*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return(m)

m=process.memory_info()
rss0=m.rss
vms0=m.vms
mss0=get_mss()

    
print("Deleting:   {}".format(D))
print("Collecting: {}".format(C))
print("                          RSS0: {:8.3f} MB, VMS0: {:8.3f} MB, MSS0: {:8.3f} MB".format(rss0/1024.**2,vms0/1024.**2,mss0/1024.**2))

def print_mem():
    m=process.memory_info()
    mss=get_mss()
    print("i: {:4.1f}, a: {:9.3f} MB, RSS: {:9.3f} MB, VMS: {:9.3f} MB, MSS: {:9.3f} MB".format(float(i)/S, float(n)*8/1024**2,(m.rss-rss0)/1024.**2,(m.vms-vms0)/1024.**2,(mss-mss0)/1024**2))

for i in range(S*M, S*N):
    n=int(pow(2,float(i)/S))
    #a = np.zeros((n))
    a = rand.standard_normal(n)
    print_mem()
    if D: del a
    if C: gc.collect(2)
   
#print("Max resident set size m.ru_maxrss {}".format(m.ru_maxrss))



#uss (Linux, OSX, Windows): aka “Unique Set Size”, this is the memory which is unique to a process and which would be freed if the process was terminated right now.
#pss (Linux): aka “Proportional Set Size”, is the amount of memory shared with other processes, accounted in a way that the amount is divided evenly between the processes that share it. I.e. if a process has 10 MBs all to itself and 10 MBs shared with another process its PSS will be 15 MBs.

#def get_proc_status(keys = None):
#    with open('/proc/self/status') as f:
#        data = dict(map(str.strip, line.split(':', 1)) for line in f)
#    return tuple(data[k] for k in keys) if keys else data

#VmHWM,VmRSS,VmPeak,VmSize,VmData = get_proc_status(('VmHWM', 'VmRSS','VmPeak', 'VmSize','VmData'))
#print('VmHWM:{:>11}\nVmRSS:{:>11}\nVmPeak:{:>10}\nVmSize:{:>10}\nVmData:{:>10}\n'.format(VmHWM,VmRSS,VmPeak,VmSize,VmData))

