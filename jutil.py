
import numpy as np
import scipy.special as sp
import scipy.optimize as opt
import scipy.interpolate as it
import scipy.fftpack as fft
#import pyfftw.interfaces.scipy_fftpack as fft
import os
import sys
import traceback
from dedalus import public as de
from dedalus import core as co
import signal
import time
import operator
import platform
import h5dict
import logging
import resource
import h5py
import math
import copy
import warnings
import psutil
import re
import ded_types as dp
from sys import stdout
#from pympler import asizeof

logger=logging.getLogger(__name__)

from mpi4py import MPI
from pathlib import Path
import jgrid

global nodenum     # Number of nodes
global nodeid      # which node we are on
global noderank    # Which process we are on this node
global nodesize    # How many processes on this node
global nodecomm    # communicator for processes on this node
global rootcomm    # communicator for root process on each node
global RSSPNODE
global VMSPNODE

global terminate
global abort
global jerror
global abort_s
global jerror_s
global dump 
global comm
global mpirank
global mpisize
global basedir
global JOB_ID
global JOB_QUEUE
global HOSTNAME 
global PID
global NODENAME
global pid
global fnmem
global fnlog
global fnlock
global fnstatus
global fnnode
global fntime
global mm
global param
global VmHWM
global VmPeak
global xxf 
global xx 
global yy 
global zz 
global Nx
global Ny
global Nz
global Nxx
global Nyy
global Nzz
global AA
global Kx
global Sx
global nstats
global fpstats
global fnsnm
global mem
global process


global gx
global gyz
global gyzJ

DED_DEBUG_LOGH = os.getenv('DED_DEBUG_LOGH') is not None
DED_DEBUG_TVDD = os.getenv('DED_DEBUG_TVDD') is not None

gyz=None
gyzJ=None
gx=None

realmin=np.finfo(np.float64).tiny
realeps=np.finfo(np.float64).eps

mem=dict()
mem['RSSHwm'] = np.int64(0)
mem['VMSHwm'] = np.int64(0)

basedir   = None
JOB_ID    = None
JOB_QUEUE = None
HOSTNAME  = None
NODENAME  = None
PID       = None
pid       = -1
VmHWM     = 0
VmPeak    = 0
mm        = 0
fnmem     = None
fnnode    = None
fnlog     = None
fnnode    = None
fnlock    = None
fnstatus  = None
fntime    = None
comm      = MPI.COMM_WORLD

mpirank   = comm.Get_rank()
mpisize   = comm.Get_size()

    
pi=math.pi
sqrtpi=np.sqrt(math.pi)

def jintenv(s):
    i=None
    if mpirank==0:
        v=os.getenv(s)
        if v is not None:
            try:
                i = int(v)
            except ValueError:
                logger.info('Environment variable "{}" is not an integer but {}'.format(s,v))
    i = comm.bcast(i,root=0) 
    return(i)

HOSTNAME           = os.getenv('HOSTNAME', 'Unknown')
JOB_ID             = os.getenv('JOB_ID',   'Unknown')
JOB_QUEUE          = os.getenv('JOB_QUEUE','Unknown')
JOB_NAME           = os.getenv('JOB_NAME', 'Unknown')

JOB_NUM_NODES      = jintenv('JOB_NUM_NODES')
TASKS_PER_NODE     = jintenv('TASKS_PER_NODE')
NTASKS_PER_NODE    = jintenv('NTASKS_PER_NODE')
NTASKS             = jintenv('NTASKS')
NPROCS             = jintenv('NPROCS')
NNODES             = jintenv('NNODES')

#JOB_NUM_TASKS         = jintenv('JOB_NUM_TASKS')
#PBS_NODENUM           = jintenv('PBS_NODENUM') 
#PBS_VNODENUM          = jintenv('PBS_VNODENUM') 
#PBS_TASKNUM           = jintenv('PBS_TASKNUM') 
#SLURM_NODEID          = jintenv('SLURM_NODEID') 
#SLURM_PROCID          = jintenv('SLURM_PROCID') 
#SLURM_LOCALID         = jintenv('SLURM_LOCALID') 
#SLURM_NTASKS_PER_NODE = jintenv('SLURM_NTASKS_PER_NODE') 
#SLURM_NNODES          = jintenv('SLURM_NNODES') 
#SLURM_NTASKS          = jintenv('SLURM_NTASKS') 
#SLURM_NPROCS          = jintenv('SLURM_NPROCS') 

#TASKS_PER_NODE        = int(re.search('[0-9]+', os.getenv('SLURM_TASKS_PER_NODE','-1')).group())
#JOB_TASKS_PER_NODE2   = int(re.search('[0-9]+', os.getenv('JOB_TASKS_PER_NODE','-1')).group())
#SLURM_TASKS_PER_NODE  = int(re.search('[0-9]+', os.getenv('SLURM_TASKS_PER_NODE','-1')).group())

pid=os.getpid()
PID=str(pid)
process = psutil.Process(pid)
NODENAME=str(platform.node())
OSX   = (sys.platform == "darwin")
LINUX = (sys.platform == "linux")


def node_info():
    global fnlog
    global nodenum     # Number of nodes
    global nodeid      # which node we are on
    global noderank    # Which process we are on this node
    global nodesize    # How many processes on this node
    global nodecomm    # communicator for processes on this node
    global rootcomm
    
    r = comm.gather(NODENAME,root=0)
    u,b=np.unique(r,return_inverse=True)
    nodenum   = comm.bcast(len(u),root=0)
    nodeid    = comm.scatter(b,root=0)
    nodecomm  = comm.Split(nodeid, mpirank)
    noderank  = nodecomm.Get_rank()
    nodesize  = nodecomm.Get_size()
    rootcomm  = comm.Split(noderank, mpirank)
    append_file(fnlog,'nodename  {}\n'.format(NODENAME))
    append_file(fnlog,'nodenum   {}\n'.format(nodenum))
    append_file(fnlog,'nodeid    {}\n'.format(nodeid))
    append_file(fnlog,'nodesize  {}\n'.format(nodesize))
    append_file(fnlog,'noderank  {}\n'.format(noderank))
 

#if mpisize>0:
#    MPI.COMM_SELF.Set_errhandler(MPI.ERRORS_ARE_FATAL)
#    MPI.COMM_WORLD.Set_errhandler(MPI.ERRORS_ARE_FATAL)
#    sys_excepthook = sys.excepthook
#    def mpi_excepthook(v, t, tb):
#        sys_excepthook(v, t, tb)
#        MPI.COMM_WORLD.Abort()
#    sys.excepthook = mpi_excepthook

terminate = False
abort = False
jerror = False
abort_s = ''
jerror_s = ''
dump  = False

#res4 = zeta**2
#res4 = res4.evaluate()
#value = res4.integrate('x','y')

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

#warnings.filterwarnings('ignore',category=FutureWarning) 
#warnings.filterwarnings('ignore',category=DeprecationWarning) 
warnings.filterwarnings('ignore')
#warnings.filterwarnings('ignore',category=h5py.H5pyDeprecationWarning) 

#warnings.showwarning = warn_with_traceback
#warnings.simplefilter("always")
#warnings.simplefilter("ignore", DeprecationWarning)
#python -W error myprogram.py  This makes all warnings fatal, see here for more information

def jdot():
    if mpirank==0:
        sys.stdout.write('.')
        sys.stdout.flush()
        
def jup():
    #if mpirank==0: sys.stdout.write("\033[F")
    if mpirank==0: sys.stdout.write("\r")
        
def jnl():
    if mpirank==0: print('')

def abool(a):
    if   a=="None":  x=None
    elif a=="False": x=np.int(0)
    elif a=="0":     x=np.int(0)
    elif a=="F":     x=np.int(0)
    elif a=="True":  x=np.int(1)
    elif a=="1":     x=np.int(1)
    elif a=="T":     x=np.int(1)
    else:
        logger.error("abool: Boolean argument must be True, False, 0 or 1")
        set_abort('abool')
    return(x)

def bool(a):
    if a: x=np.int(1)
    else: x=np.int(0)
    return(x)

def Barrier():
    comm.Barrier()

def set_base_dir(sType,name):
    global basedir
    global fnlock
    global fnstatus
    global fntime
    DD   = Path(os.environ['DEDALUS_DATA'])
    bdir = DD / sType
    basedir  = bdir / name
    jfinal=basedir / 'final' / 'jfinal.hdf5'
    if mpirank==0:
        fnlock   = basedir / "lock"
        fnstatus = basedir / "status"
        fntime   = basedir / "time"
        if not bdir.is_dir():    bdir.mkdir(parents=True)
        if not basedir.is_dir(): basedir.mkdir(parents=True)
        if not 'JOB_ID' in  os.environ:
            fnout    = basedir / '{}.out'.format(PID)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            fh = logging.FileHandler(str(fnout))
            fh.setLevel(logging.INFO)
            fh.setFormatter(formatter)
            logger.addHandler(fh)
            logger.info("fnout:    {}".format(fnout))
        logger.info("DEDALUS_DATA: {}".format(DD))
        logger.info("bdir:         {}".format(bdir))
        logger.info("ddir:         {}".format(basedir))
        logger.info("fnstatus:     {}".format(fnstatus))
        logger.info("fnlock:       {}".format(fnlock))
        logger.info("fntime:       {}".format(fntime))
        logger.info("version:      {}".format(sys.version))
    comm.Barrier()
    return(bdir,basedir)

def check_version(v):
    if mpirank==0:
        fn=basedir / "version"
        if not fn.is_file(): write_file(fn,"{:.2f}".format(v))
        else:
            f=fn.open('r')
            v2=f.read()
            f.close()
            v2=np.float(v2)
            if v2>v:
                logger.error("Read version {:.2f} but program is version {:.2f} ".format(v2,v))
                set_error('version')
            else:
                logger.info("Program version {:.2f} restarting from {:.2f} ".format(v,v2))
    check_error()

def write_status(s):
    global fnstatus
    if mpirank==0 and basedir is not None: write_file(fnstatus,s  + '\n')
        
def write_time(t):
    global fntime
    if mpirank==0 and basedir is not None: write_file(fntime,"{:7.3f}\n".format(t))
    
def write_param(p):
    if mpirank==0:
        if type(p) != type(dict()):
            logger.info('write_param: p must be a dict not "{}"'.format(type(p)))
            set_abort('write_param')
        else:
            fnp=basedir / 'param.h5'
            h5dict.save_dict_to_hdf5(p,fnp)
            dp.write_options(p,basedir / 'options.txt')

def write_dict(fn,p):
    if mpirank==0:
        if type(p) != type(dict()):
            logger.info('write_dict: p must be a dict not "{}"'.format(type(p)))
            set_abort('write_dict')
        else:
            fnp=basedir / fn
            h5dict.save_dict_to_hdf5(p,fnp)

def read_param(fnp):
    global basedir
    p=None
    if mpirank==0:
        if not fnp:
            logger.error("read_param called with empty argument")
            set_error('read_param')
        try:
            logger.info("Reading parameters from {}".format(fnp))
            p = h5dict.load_dict_from_hdf5(fnp)
        except:
            logger.error("Cannot read parameters from {}".format(fnp))
            write_file(basedir / "error","read_param" + '\n')
            set_error('read_param')
    check_error()
    p = comm.bcast(p, root=0)
    return(p)

def jsum(x):
    y=x.sum(axis=None)
    y=comm.allreduce(y, op=MPI.SUM)
    return(y)

def MemPhys():
    global mem
    pmem = psutil.virtual_memory()
    mem_min  = comm.allreduce(pmem.total, op=MPI.MIN)         # min  physical memory available
    mem_mean = comm.allreduce(pmem.total, op=MPI.SUM)/mpisize # mean physical memory available
    mem_max  = comm.allreduce(pmem.total, op=MPI.MAX)         # max  physical memory available
    if mpirank==0: logger.info('Physical memory: [{:} {:} {:}]'.format(MemStr(mem_min),MemStr(mem_mean),MemStr(mem_max)))
    mem['PHYS']=mem_min
    #resource.setrlimit(resource.RLIMIT_VMEM, (mem_mean, mem_mean))
    
def jmax(x,c=comm):
    y=np.NINF
    if x is not None:
        jt=jtype(x)
        if jt[0]=='f':   x=x['g']
        if jt=='scalar': y=x
        elif x.size>0:   y=x.max()
    y=comm.allreduce(y, op=MPI.MAX)
    return(y)

def jmin(x,c=comm):
    y=np.PINF
    if x is not None:
        jt=jtype(x)
        if jt[0]=='f':   x=x['g']
        if jt=='scalar': y=x
        elif x.size>0:   y=x.min()
    y=comm.allreduce(y, op=MPI.MIN)
    return(y)

def jrange(x):
    return(jmin(x),jmax(x))

def prange(s,x):
    if x is None:
        r=''
    else:
        if type(x) is np.ndarray: minx,maxx=jrange(x)
        else:                     minx,maxx=jrange(x['g'])
        if not np.isfinite(minx) or not np.isfinite(maxx):
            set_error('prange: ' + s + ' not finite')
            r = ', {}: not finite:'.format(s)
        else : r = ', {}:[{:6.3f},{:6.3f}]'.format(s,minx,maxx)
    return(r,minx,maxx)

def srange(s,x):
    n,x=jrange(x)
    logger.info(s.format(n,x))
    

def touch(fname):
    fname.touch()
    
def abort_run():
    logger.info("Abort triggered")
    write_status("Exception")
    delete_lock()

    #if mpisize>1: comm.Abort()
    exit(1)


def write_files():
    global fnnode
    global fnmem
    global fnlog
    debug = False
    dnode = basedir / "nodes"
    dmem  = basedir / "mem"
    dlog  = basedir / "log"
    if mpirank==0:
        fnt=basedir / "time"
        fnt.touch()
        write_file(basedir / "ntasks",str(mpisize))
        write_file(basedir / "queue",JOB_QUEUE + '\n')
        write_file(basedir / "jobid",JOB_ID + '\n')
        fenv=basedir / "env"
        f=fenv.open('w')
        for key in os.environ.keys():
            f.write(key + "= " + os.environ[key] +"\n")
        f.close()
        if not dnode.is_dir(): dnode.mkdir()
        if not dmem.is_dir():  dmem.mkdir()
        if not dlog.is_dir():  dlog.mkdir()
        #cmd=os.path.basename(__file__)
        #write_file(basedir/ "cmd",cmd + '\n')
    comm.Barrier() # make sure directories are created
   
    hn=NODENAME + " " + PID
    fnnode = dnode / '{:04d}'.format(mpirank)
    fnmem  = dmem  / '{:04d}'.format(mpirank)
    fnlog  = dlog  / '{:04d}'.format(mpirank)
    write_file(fnnode,'NODENAME {}\n'.format(NODENAME))
    write_file( fnlog, 'mpirank/mpisize {}/{}\n'.format(mpirank,mpisize))
    append_file(fnlog, 'HOSTNAME              {}\n'.format(HOSTNAME))
    append_file(fnlog, 'NODENAME              {}\n'.format(NODENAME))
    #append_file(fnlog, 'PBS_NODENUM           {}\n'.format(PBS_NODENUM))
    #append_file(fnlog, 'PBS_VNODENUM          {}\n'.format(PBS_VNODENUM))
    #append_file(fnlog, 'PBS_TASKNUM           {}\n'.format(PBS_TASKNUM))
    #append_file(fnlog, 'SLURM_NODEID          {}\n'.format(SLURM_NODEID))          
    #append_file(fnlog, 'SLURM_PROCID          {}\n'.format(SLURM_PROCID))         
    #append_file(fnlog, 'SLURM_LOCALID         {}\n'.format(SLURM_LOCALID))         
    #append_file(fnlog, 'SLURM_NTASKS_PER_NODE {}\n'.format(SLURM_NTASKS_PER_NODE)) 
    #append_file(fnlog, 'SLURM_NNODES          {}\n'.format(SLURM_NNODES))          
    #append_file(fnlog, 'SLURM_NTASKS          {}\n'.format(SLURM_NTASKS))         
    #append_file(fnlog, 'SLURM_NPROCS          {}\n'.format(SLURM_NPROCS))          
    logger.info('')
    logger.info('Environment variables')
    logger.info('  JOB_ID:            {}'.format(JOB_ID))         
    logger.info('  JOB_QUEUE:         {}'.format(JOB_QUEUE))       
    logger.info('  HOSTNAME:          {}'.format(HOSTNAME))         
    logger.info('  NODENAME:          {}'.format(NODENAME))         
    logger.info('  JOB_NUM_NODES      {}'.format(JOB_NUM_NODES))     
    logger.info('  TASKS_PER_NODE     {}'.format(TASKS_PER_NODE))    
    logger.info('  NTASKS_PER_NODE    {}'.format(NTASKS_PER_NODE))   
    logger.info('  NTASKS             {}'.format(NTASKS))            
    logger.info('  NPROCS             {}'.format(NPROCS))            
    logger.info('  NNODES             {}'.format(NNODES))            
    
    MemoryUsage()
    check_error()
    comm.Barrier()
    
def write_file(fn,msg):
    if fn.is_dir():
        print('write_file: "{}" is a directory \n "{}"'.format(fn,msg))
        set_error('write file: "' + str(fn) + '"')
        return
    try:
        f=fn.open('w')
        f.write(msg)
        f.close()
    except:
        print('write_file: Failed to write to "{}" \n "{}"'.format(fn,msg))
        set_error('write file: "' + str(fn) + '"')
     
def append_file(fn,msg):
    if fn.is_dir():
        print('write_file: "{}" is a directory \n "{}"'.format(fn,msg))
        set_error('write file: "' + str(fn) + '"')
        return
    try:
        f=fn.open('a')
        f.write(msg)
        f.close()
    except:
        print('append_file: Failed to append to {} \n {}'.format(fn,msg))
        set_error('append file: "' + str(fn) + '"')

    
def make_lock():
    global fnlock
    if mpirank==0:
        if fnlock.is_file():
            set_error('make_lock: lock file exists')
            logger.info("make_lock: lock file exists {}".format(fnlock))
        try:
            f=fnlock.open('w')
            f.write(HOSTNAME + ' ' + str(PID))
            f.close()
            logger.info("make_lock: {}".format(fnlock))
        except:
            set_error('make_lock: create')
            logger.info("make_lock: failure {}".format(fnlock))
    test_error()
    if jerror: exit(1)
    
def delete_lock():
    global fnlock
    if mpirank==0:
        if not fnlock.is_file(): logger.info("Lock file does not exist {}".format(fnlock))
        else:
            fnlock.unlink()
            logger.info("Lock file deleted {}".format(fnlock))
    
def get_sec(time_str):
    try:
        h, m, s = time_str.split(':')
    except:
         d,h,m,s= time_str.split(':')
         h=h+24*d
         
    t = np.float(int(h) * 3600 + int(m) * 60 + int(s))
    if t==0:
        t=np.inf
    return(t)



# Find index of nearest value in an array using bisection
def find_nearest(array, value):
    idx_sorted = np.argsort(array)
    sorted_array = np.array(array[idx_sorted])
    idx = np.searchsorted(sorted_array, value, side="left")
    if idx >= len(array):
        idx_nearest = idx_sorted[len(array)-1]
    elif idx == 0:
        idx_nearest = idx_sorted[0]
    else:
        if abs(value - sorted_array[idx-1]) < abs(value - sorted_array[idx]):
            idx_nearest = idx_sorted[idx-1]
        else:
            idx_nearest = idx_sorted[idx]
    return idx_nearest


def linear_fit_calc(x,y,X): # Find a linear fit to y=f(x)=a*x+b then calculate f(X)
    Ix1 = (x).sum()
    Iy1 = (y).sum()
    Ixy = (x*y).sum()
    Ixx = (x*x).sum()
    I11 = x.size
    d   =  I11*Ixx-Ix1*Ix1
    a   = (I11*Ixy-Ix1*Iy1)/d
    b   = (Ixx*Iy1-Ix1*Ixy)/d
    Y   = b+a*X
    return(Y)

def linear_fit_calc_test():
    a =  np.random.normal()
    b =  np.random.normal()
    x =  np.random.normal(size=[10,])
    y = a+b*x
    ny=linear_fit_calc(x,y,x)
    print('{:}'.format(y-ny))
    

              
 
def switch(x,D,A1,A2):
    if D==0:
        y = (A1+A2+(A1-A2)*np.sign(x))/2
        y = np.round(y)
    else:
        y = (A1+A2+(A1-A2)*sp.erf(x*D))/2
    return(y)

def indicator(x,x1,x2,D):
    if D==0:
        z = (np.sign(x-x1)-np.sign(x-x2))/2
        z = np.round(z)
    else:
        z = (sp.erf((x-x1)*D)-sp.erf((x-x2)*D))/2
    return(z)

def indicatorp(x,x1,x2,D,L):
    z=np.maximum(0,np.minimum(1,indicator(x,x1,x2,D) + indicator(x-L,x1,x2,D) + indicator(x+L,x1,x2,D)))
    return(z)

def heaviside(x,D):
    if D==0:
        z = (1+np.sign(x))/2
        z = np.round(z)
    else:
        z = (1+sp.erf(x*D))/2
    return(z)


def int_weights2(X,N,T,j,domain):
    
    if T=="Cheb":
        ls=domain.dist.grid_layout.slices(scales=1)[j]
        I=np.ndarray((N,),dtype=np.float)
        x=de.Chebyshev('x', N, interval=(0,X), dealias=1)
        if mpirank==0: comm0 = comm.Split(0,0)
        else:          comm0 = comm.Split(1,mpirank)
        if mpirank==0:
            d=de.Domain([x], grid_dtype=np.float64,comm=comm0)
            f=d.new_field()
            for k in range(N):
                f['g']=0
                f['g'][k]=1
                I[k] = np.float(f.integrate('x')['g'][0])

        if mpirank==0: logger.info('Integration weight: {}, interval: {}'.format(sum(I),X))
        I = comm.bcast(I, root=0)
        I = I[ls]        
    else: I=X/N
    j=j+1
    return(I,j)

def int_weights(AA,domain):
    dim=range(domain.dim)
    Il=list(dim)
    Ig=list(dim)
    for i in dim:
        N=domain.dist.grid_layout.global_shape(scales=AA)[i]
        b=domain.bases[i]
        X=b.interval[1]-b.interval[0]
        T=get_bases_type(b)
        #logger.info('{},Length {}, type: {}'.format(i,X,T))
        if T=="Chebyshev":
            ls=domain.dist.grid_layout.slices(scales=AA)[i]
            I=np.ndarray((N,),dtype=np.float)
            if mpirank==0:
                k=np.arange(1,1+np.int(N/2))
                for j in range(N): I[j] = 1-2*np.sum(np.cos(k*(2*j+1)*np.pi/N)/(4*k**2-1))
                I=I*X/N
                logger.info('Integration weight: {}, interval: {}'.format(sum(I),X))
            Ig[i] = comm.bcast(I, root=0)
            Il[i] = Ig[i][ls]        
        elif T=="Fourier":
            Il[i]=X/N
            Ig[i]=X/N
        elif T=="SinCos":
            Il[i]=X/N
            Ig[i]=X/N
        else:
            logger.info('int_weights: unknown coordinate type {}'.format(T))
            set_error('int_weights: unknown coordinate type {}'.format(T))
            check_error()
    return(Il,Ig)

#def basis_init(I,T,N,A):
#    n=I.shape((0,))
#    b=list()
#    if n==3: c=('x','y,'z'):
#    else:    c=('x','z'):
#    for k in range(n):
#        
#    if   T[:,k]=="SinCos":    b.append=de.SinCos(   c[k], N[k], interval=I[k], dealias=A[k])
#    elif T[:,k]=="Fourier":   b.append=de.Fourier(  c[k], N[k], interval=I[k], dealias=A[k])
#    elif T[:,k]=="Chebyshev": b.append=de.Chebyshev(c[k], N[k], interval=I[k], dealias=A[k])
#    elif T[:,k]=="Cheb":      b.append=de.Chebyshev(c[k], N[k], interval=I[k], dealias=A[k])
#    logger.info('basis_init: unknown coordinate type {}'.format(T:,k)
def slice_intersect(a,b):
    n=len(a)
    c=list()
    for k in range(n): c.append(slice(max(a[k].start,b[k].start),min(a[k].stop,b[k].stop)))
    return(c)

def boussinseq_init(p):
    global xxf 
    global xx 
    global yy 
    global zz 
    global Nx
    global Ny
    global Nz
    global Nxx
    global Nyy
    global Nzz
    global AA
    global param
    global Kx
    global Sx
    global gyz
    global gyzJ
    global gx
    
    L     = p['Length']
    W     = p['Width']
    H     = p['Height']
    Nx    = p['Nx']
    Ny    = p['Ny']
    Nz    = p['Nz']
    Tx    = dp.get_param(p,'Tx')
    Ty    = dp.get_param(p,'Ty')
    Tz    = dp.get_param(p,'Tz')
    AA    = p['AA']
    sType = p['sType']

    U     = dp.get_param(p,'Velocity')
    oddg  = dp.get_bool_param(p,'oddg')
    fpidu = dp.get_bool_param(p,'fpidu')
    fpidv = dp.get_bool_param(p,'fpidv')
    fpidw = dp.get_bool_param(p,'fpidw')
    fpidb = dp.get_bool_param(p,'fpidb')
    fpids = dp.get_bool_param(p,'fpids')
    
    Conservative = dp.get_bool_param(p,'Conservative')
    average      = dp.get_bool_param(p,'avrg')
    intdx        = dp.get_bool_param(p,'intdx')

    B     = dp.get_param(p,'Buoyancy')
    S     = dp.get_param(p,'Sediment')
    noise = dp.get_param(p,'Noise')
    VMEM  = dp.get_param(p,'VMEM')

    LL=(0,L)
    WW=(0,W)
    HH=(0,H)
    if sType=='pm':
        if Ty != "SinCos": WW=(-W/2,W/2)
        if Tz != "SinCos": HH=(-H/2,H/2)
     
    param={}   
    param['cm']=dp.get_bool_param(p,'cm')
    param['U']=U
    param['L']=L
    param['W']=W
    param['H']=H
    param['Lx']=LL
    param['Ly']=WW
    param['Lz']=HH
    param['Nx']=Nx
    param['Ny']=Ny
    param['Nz']=Nz
    param['Rx']=LL[1]-LL[0]
    param['Ry']=WW[1]-WW[0]
    param['Rz']=HH[1]-HH[0]
    param['Volume']=(LL[1]-LL[0])*(WW[1]-WW[0])*(HH[1]-HH[0])
    param['NN']=Nx*Ny*Nz
    # Create bases and domain

    
    if   Tx=="SinCos":  x_basis=de.SinCos(   'x', Nx, interval=LL, dealias=AA)
    elif Tx=="Fourier": x_basis=de.Fourier(  'x', Nx, interval=LL, dealias=AA)
    elif Tx=="Cheb":    x_basis=de.Chebyshev('x', Nx, interval=LL, dealias=AA)
    if   Ty=="SinCos":  y_basis=de.SinCos(   'y', Ny, interval=WW, dealias=AA)
    elif Ty=="Fourier": y_basis=de.Fourier(  'y', Ny, interval=WW, dealias=AA)
    elif Ty=="Cheb":    y_basis=de.Chebyshev('y', Ny, interval=WW, dealias=AA)
    if   Tz=="SinCos":  z_basis=de.SinCos(   'z', Nz, interval=HH, dealias=AA)
    elif Tz=="Fourier": z_basis=de.Fourier(  'z', Nz, interval=HH, dealias=AA)
    elif Tz=="Cheb":    z_basis=de.Chebyshev('z', Nz, interval=HH, dealias=AA)
       
       
    # Dedalus will only work with 2d mesh decompositions
       
    param['NN']=Nx*Ny*Nz
    if Ny>1: N1=Ny
    else:    N1=Nx
    N2=Nz
    m1=None
    m2=None
    if Ny==1:
        m1=mpisize
        m2=1
        p.pop('m1', None)
        p.pop('m2', None)
    else:
        if 'm1' in p: m1=p['m1']
        if 'm2' in p: m2=p['m2']
        if m1 is None and m2 is None: m1,m2=meshxy(mpisize,N1,N2)
        if m1 is None and not (m2 is None): m1=mpirank/m2
        if m2 is None and not (m1 is None): m2=mpirank/m1
        if m1*m2 != mpisize and mpirank==0:
            logger.info('boussinseq_init: grid size [{} {}] does not match mpi size {}'.format(m1,m2,mpisize))
            set_error('boussinseq_init: grid size does not match mpi size')
        check_error()


    dx=np.mean(np.diff(np.asarray(x_basis.grid(scale=1))))
    if Ny>1: dy=np.mean(np.diff(np.asarray(y_basis.grid(scale=1))))
    else:    dy=1
    dz=np.diff(np.asarray(z_basis.grid(scale=1)))
    mindz=np.min(dz)
    maxdz=np.max(dz)
    dz=np.mean(dz)
    
    xx = x_basis.grid(scale=AA)
    if Ny>1: yy = y_basis.grid(scale=AA)
    else: yy=None
    if round(AA*Nz) != AA*Nz:
        logger.error('L: {} W: {} H: {}'.format(L,W,H))
        logger.error('AA: {} Nx: {} Ny: {} Nz: {} AA*Nz; {}. Not an integer'.format(AA,Nx,Ny,Nz,AA*Nz))
        set_error('AA*Nz is not an integer')
    check_error()
 
    zz = z_basis.grid(scale=AA)
    Nxx=len(xx)
    if Ny>1: Nyy=len(yy)
    else:    Nyy=1
    Nzz=len(zz)
    xxf=xx
    param['NNN']=Nxx*Nyy*Nzz
    if Ny>1:
        xx=xx.reshape([Nxx,1,1])
        yy=yy.reshape([1,Nyy,1])
        zz=zz.reshape([1,1,Nzz])
    else:
        xx=xx.reshape([Nxx,1])
        zz=zz.reshape([1,Nzz])

    # Calculate integration weights
    if Ny>1: dim=3
    else: dim=2



            
    if Ny>1:
        domain=de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64, mesh=(m1,m2))
        logger.info('Mesh set to [{} {}], actual [{} {}]'.format(m1,m2,domain.dist.grid_layout.ext_mesh[1],domain.dist.grid_layout.ext_mesh[2]))
        
    else:    domain=de.Domain([x_basis, z_basis], grid_dtype=np.float64)

    #print_domain_info(domain,1)
    
    Kx = domain.elements(0)*L/(np.pi*(Nx-1)) # scale it to [0 1]
    Sx = np.exp(-(5*Kx)**2) # Damps highest mode by exp(-25)=1.4e-11 
    
    def BFF(*args):
        t = args[0].data[()] # this is a scalar; we use .data[()] to get its value
        x = args[1].data  # this is an array; we use .data to get its values
        y = args[2].data  # this is an array; we use .data to get its values
        h = args[3].data[()]
        A = args[4].data[()]
        rs=(x*x+y*y)/h**2
        return(2*A*np.exp(-2*rs))
    
    def BF(*args, domain=domain, F=BFF):
        return(de.operators.GeneralFunction(domain, layout='g', func=F, args=args))

    de.operators.parseables['BF'] = BF


    if Ny>1:   s=['u','v','w']
    else:      s=['u','w']
    if B is not None:   s=s + ['b']             # Add on buoyancy field
    if S is not None:   s=s + ['s']             # Add on sediment field
    v=['p']

    
    if average:
        if sType=='gc':
            v = v + ['ap'] + ['a' + ss for ss in s] + ['auu','aww','auw','abu','abw']
            if Ny>1: v = v + ['avv']
        elif  sType=='pm':
            v = v + ['ap','aur'] + ['a' + ss for ss in s] + ['auu','arr','aru']
            v.remove('aw')
    if noise is not None:
        s = s + ['noisef1']        
        v = v + ['noisef12']        
    if 'av' in v: v.remove('av')

    if fpidu: s = s + ['fpidu']
    if fpidv: s = s + ['fpidv']
    if fpidw: s = s + ['fpidw']
    if fpidb: s = s + ['fpidb']
    if fpids: s = s + ['fpids']
    if Tz=="Cheb": v = v + [ss + 'z' for ss in s]

    v = v+s
    
    problem = de.IVP(domain, variables=v, time='t')
    v=problem.variables
    
    if Tx=="SinCos":
        ve=set(('fpidv','fpidw','p','v','vz','w','wz','av','aw','auu','avv','aww','avw','abw','ap','noisef1','noisef2','noisef1z','wu'))
        vo=set(('fpidu','u','fu','uz','au','auw','auv','abu'))
        if sType=='pm' and p['Forcing']>6 and not oddg:
            vo=vo | set(('fpidb','fpids','b','s','sz','bz','ab','as'))
            ve=ve | set(('abu',))
        else:
            ve=ve | set(('fpidb','fpids','b','s','sz','bz','ab','as'))
            vo=vo | set(('abu',))
        ve=list(set(v) & ve)
        vo=list(set(v) & vo)
        
        for vv in ve: problem.meta[vv]['x']['parity'] = 1
        for vv in vo: problem.meta[vv]['x']['parity'] = -1
        
        #ve=list(set(v) & set(('fpidv','fpidw','fpidb','fpids','p','b','s','sz','bz','v','vz','w','wz','ap','ab','av','aw','auu','avv','aww','avw','abw','noisef1','noisef2','noisef1z','wu')))
        #vo=list(set(v) & set(('fpidu','u','fu','uz','au','auw','auv','abu')))
        #        if mpirank==0:
        #            logger.info('Tx ve {}'.format(ve))
        #            logger.info('Tx vo {}'.format(vo))
        for vv in ve: problem.meta[vv]['x']['parity'] = 1
        for vv in vo: problem.meta[vv]['x']['parity'] = -1
    if Ty=="SinCos" and Tz=="SinCos":  
        vye=list(set(v) & set(('fpidu','fpidw','fpidb','fpids','p','b','s','sz','bz','u','w','fu','uz','ap','ab','au','aw','auu','avv','aww','abu','noisef1','w','fw','wz','aw','auw','abw')))
        vze=list(set(v) & set(('fpidv','fpidu','fpidb','fpids','p','b','s','sz','bz','u','v','fu','uz','ap','ab','au','aw','auu','avv','aww','abu','noisef1','v','fv','vz','av','auv','abw')))
        vyo=list(set(v) & set(('fpidv','v','fv','vz','av','auv','auv','abv','noisef1z','wv')))
        vzo=list(set(v) & set(('fpidw','w','fw','wz','aw','auw','auw','abw','noisef1z','wv')))
        for vv in vye: problem.meta[vv]['y']['parity'] = 1
        for vv in vze: problem.meta[vv]['z']['parity'] = 1
        for vv in vyo: problem.meta[vv]['y']['parity'] = -1
        for vv in vzo: problem.meta[vv]['z']['parity'] = -1
    elif Ty=="SinCos": 
        ve=list(set(v) & set(('fpidu','fpidw','fpidb','fpids','p','b','s','sz','bz','u','fu','uz','w','wz','ap','ab','au','aw','auu','avv','aww','auw','abw','abu','noisef1','noisef1z','wv')))
        vo=list(set(v) & set(('fpidv','v','fv','vz','av','avw','auv','abv')))
        #        if mpirank==0:
        #            logger.info('Ty ve {}'.format(ve))
        #            logger.info('Ty vo {}'.format(vo))
        for vv in ve: problem.meta[vv]['y']['parity'] = 1
        for vv in vo: problem.meta[vv]['y']['parity'] = -1
    elif Tz=="SinCos": 
        ve=list(set(v) & set(('fpidu','fpidv','fpidb','fpids','p','b','s','sz','bz','u','fu','uz','v','vz','ap','ab','au','av','auu','avv','aww','auv','abv','abu','noisef1','noisef1z','ww')))
        vo=list(set(v) & set(('fpidv','w','fw','wz','aw','aww','auw','abw')))
        for vv in ve: problem.meta[vv]['z']['parity'] = 1
        for vv in vo: problem.meta[vv]['z']['parity'] = -1
  #  for vv in v: print("{:5s} x {} z {}".format(vv,problem.meta[vv]['x']['parity'],problem.meta[vv]['y']['parity']))
 
    problem.substitutions['rho'] = '1'   # Boussinesq simulation
    if sType=='pm': problem.substitutions['ur'] = '((x*u+y*v)/sqrt(x*x+y*y))'       # r velocity component
    # Enstrophy

    #    problem.substitutions['EN']   = 'dx(u)**2+dx(v)**2+dx(w)**2+dy(u)**2+dy(v)**2+dy(w)**2+dz(u)**2+dz(v)**2+dz(w)**2'
    if Tz != "Cheb":
        for ss in ('u','v','w','b','s'):
            if ss in problem.variables: problem.substitutions[ss + 'z'] = '0'

    if not 'u' in problem.variables: problem.substitutions['u'] = '0'
    if not 'v' in problem.variables: problem.substitutions['v'] = '0'
    if not 'w' in problem.variables: problem.substitutions['w'] = '0'

    if not 'u' in problem.variables: problem.substitutions['dx(X)']='0'
    if not 'v' in problem.variables: problem.substitutions['dy(X)']='0'
    if not 'w' in problem.variables: problem.substitutions['dw(X)']='0'
    
    
    problem.substitutions['vorx'] = 'dy(w)-dz(v)'       # x vorticity
    problem.substitutions['vory'] = 'dz(u)-dx(w)'       # y vorticity
    problem.substitutions['vorz'] = 'dx(v)-dy(u)'       # z vorticity
    problem.substitutions['SR']   = '2*dx(u)**2+2*dy(v)+2*dz(w)**2+(dx(v)+dy(u))**2+(dx(w)+dz(u))**2+(dy(w)+dz(v))**2' #Total shear
    problem.substitutions['SS']   = 'u*u + v*v + w*w'
    problem.substitutions['E']    = 'vorx*vorx+vory*vory+vorz*vorz'
    problem.substitutions['P']    = 'p+(v*v+u*u+w*w)/2' # Total pressure
    problem.substitutions['WO']   = 'dx(u)**2+dy(v)**2+dz(w)**2+2*dy(w)*dz(v)+2*dx(w)*dz(u)+2*dx(v)*dy(u)' # Weiss-Okubo

    problem.substitutions['POS(X)'] = '(X+abs(X))/2'
    problem.substitutions['SGN(X)'] = '(X/max(1e-50,abs(X)))'
    problem.substitutions['PIN(X)'] = '(0.5+X/(2*abs(X)))'


    if Conservative:
        cdm           = '(  dx(u*X)'
        if Ny>1: cdm += ' + dy(v*X)' 
        cdm           += ' + dz(w*X)' 
        cdm+=' )'
        problem.substitutions['CDM(X,Y)'] = cdm
        
    cd                 ='(   u*dx(X)'
    if Ny>1:       cd += ' + v*dy(X)' 
    if Tz=="Cheb": cd += ' + w*Y    '
    else:          cd += ' + w*dz(X)'
    cd+=' )'
    problem.substitutions['CD(X,Y)'] = cd
    
    lp                 = '(  dx(dx(X))'
    if Ny>1:       lp += ' + dy(dy(X))'
    if Tz=="Cheb": lp += ' + dz(Y)    '
    else:          lp += ' + dz(dz(X))'
    lp+=' )'
    problem.substitutions['LP(X,Y)'] = lp
    
    check_error()

    #logger.info("variables:    {}".format(problem.variables))
    #logger.info("variables: {}".format(v))
    MemPhys()
    logger.info("AA:    {:7.2f}, AAS:    {:7.2f}, AAJ:    {:7.2f}".format(AA,p['AAS'],p['AAJ']))
    if Ny>1:
        logger.info("Tx:   {:8s},    Ty: {:8s}, Tz: {:8s}".format(Tx,Ty,Tz))
        logger.info("x  rng: [{:7.2f},{:7.2f}] {}".format(LL[0],LL[1],Tx))
        logger.info("y  rng: [{:7.2f},{:7.2f}] {}".format(WW[0],WW[1],Ty))
        logger.info("z  rng: [{:7.2f},{:7.2f}] {}".format(HH[0],HH[1],Tz))
        logger.info("dx:     {:8.5f},  dy: {:8.5f}".format(dx,dy))
        logger.info("dz range: [{:8.5f} {:8.5f}]".format(mindz,maxdz))
        logger.info("Nx:    {:7d},    Ny: {:7d},    Nz: {:7d}".format(Nx, Ny, Nz))
        logger.info("Nxx:   {:7d},   Nyy: {:7d},   Nzz: {:7d}".format(Nxx,Nyy,Nzz))
        logger.info("m1:    {:7d},    m2: {:7d}".format(m1, m2))
        logger.info("N1/m1: {:7.2f}, N2/m2: {:7.2f}".format(N1/m1, N2/m2))
    else:
        logger.info("Tx:   {:8s},    Tz: {:8s}".format(Tx,Tz))
        logger.info("x  rng: [{:7.2f},{:7.2f}] {}".format(LL[0],LL[1],Tx))
        logger.info("z  rng: [{:7.2f},{:7.2f}] {}".format(HH[0],HH[1],Tz))
        logger.info("dx:     {:8.5f},  dz range: [{:8.5f} {:8.5f}]".format(dx,mindz,maxdz))
        logger.info("Nx:    {:7d},    Nz: {:7d}".format(Nx, Nz))
        logger.info("Nxx:   {:7d},   Nzz: {:7d}".format(Nxx,Nzz))
        logger.info("m1:    {:7d},   Nz/m1: {:7.2f}".format(m1, Nz/m1))
    logger.info("Grid cells                 {:8.3f} {:8.3f} (million)".format(param['NN']/1e6,param['NNN']/1e6))
    N1=np.prod(domain.dist.grid_layout.local_shape(scales=1))
    N2=np.prod(domain.dist.grid_layout.local_shape(scales=AA))
    MNM1=comm.allreduce(N1, op=MPI.SUM)/mpisize/1e3
    MIN1=comm.allreduce(N1, op=MPI.MIN)/1e3
    MAX1=comm.allreduce(N1, op=MPI.MAX)/1e3
    MNM2=comm.allreduce(N2, op=MPI.SUM)/mpisize/1e3
    MIN2=comm.allreduce(N2, op=MPI.MIN)/1e3
    MAX2=comm.allreduce(N2, op=MPI.MAX)/1e3
    logger.info('Cells/core {:3.1f}: min={:8.3f} mean={:8.3f} max={:8.3f} (thousand)'.format( 1,MIN1,MNM1,MAX1))
    logger.info('Cells/core {:3.1f}: min={:8.3f} mean={:8.3f} max={:8.3f} (thousand)'.format(AA,MIN2,MNM2,MAX2))
    if VMEM is not None: logger.info("VMEM:           {:6.1f}GB".format(VMEM))
    append_file(fnlog,'domain.dist.rank:     {}\n'.format(domain.dist.rank))   
    append_file(fnlog,'domain.dist.size:     {}\n'.format(domain.dist.size))
    append_file(fnlog,'domain.dist.coords:   {}\n'.format(domain.dist.coords)) 
    append_file(fnlog,'domain.dist.mesh:     {}\n'.format(domain.dist.mesh))  
    append_file(fnlog,'domain.grid(0).shape: {}\n'.format(domain.grid(0).shape))
    append_file(fnlog,'domain.grid(1).shape: {}\n'.format(domain.grid(1).shape))
    if Ny>1: append_file(fnlog,'domain.grid(2).shape: {}\n'.format(domain.grid(2).shape))
    append_file(fnlog,'                   \n')
    append_file(fnlog,'domain.dist.grid_layout\n')
    append_file(fnlog,'  ext_coords     {}\n'.format(domain.dist.grid_layout.ext_coords))  
    append_file(fnlog,'  ext_mesh       {}\n'.format(domain.dist.grid_layout.ext_mesh))    
    append_file(fnlog,'  index          {}\n'.format(domain.dist.grid_layout.index))       
    append_file(fnlog,'  local          {}\n'.format(domain.dist.grid_layout.local))       
    append_file(fnlog,'  grid_space     {}\n'.format(domain.dist.grid_layout.grid_space))  
    append_file(fnlog,'                   \n')
    append_file(fnlog,'  scales         {}\n'.format(1))
    append_file(fnlog,'    global_shape {}\n'.format(domain.dist.grid_layout.global_shape(scales=1))) 
    append_file(fnlog,'    local_shape  {}\n'.format(domain.dist.grid_layout.local_shape(scales=1))) 
    append_file(fnlog,'    slices       {}\n'.format(domain.dist.grid_layout.slices(scales=1)))         
    append_file(fnlog,'    start        {}\n'.format(domain.dist.grid_layout.start(scales=1)))           
    append_file(fnlog,'                   \n')
    append_file(fnlog,'  scales         {}\n'.format(AA))
    append_file(fnlog,'    global_shape {}\n'.format(domain.dist.grid_layout.global_shape(scales=AA))) 
    append_file(fnlog,'    local_shape  {}\n'.format(domain.dist.grid_layout.local_shape(scales=AA))) 
    append_file(fnlog,'    slices       {}\n'.format(domain.dist.grid_layout.slices(scales=AA)))         
    append_file(fnlog,'    start        {}\n'.format(domain.dist.grid_layout.start(scales=AA)))           
    append_file(fnlog,'                   \n')
    append_file(fnlog,'domain.dim:      {}\n'.format(domain.dim))
#    append_file(fnlog,'domain.grid:     {}\n'.format(domain.grid))
#    append_file(fnlog,'domain.grids:    {}\n'.format(domain.grids))
#    append_file(fnlog,'vars(domain.grid): {}\n'.format(vars(domain.grid)))
#    append_file(fnlog,"domain.grid['axis']: {}\n".format(vars(domain.grid['axis'])))
#    append_file(fnlog,'vars(domain.grids):{}\n'.format(vars(domain.grids)))
#    append_file(fnlog,'domain.grids.axis:       {}\n'.format(domain.grids.axis))
#    append_file(fnlog,'vars(domain):      {}\n'.format(vars(domain)))
#    append_file(fnlog,'vars(distributor): {}\n'.format(vars(domain.distributor)))
#    append_file(fnlog,'vars(dist):        {}\n'.format(vars(domain.dist)))
 
    #   logger.info('vars(domain.dist.layout_references: {}'.format(dir(domain.dist.layout_reference)))
    if False:
     #    logger.info("vars(domain):      {}".format(vars(domain)))
    #    logger.info("vars(distributor): {}".format(vars(domain.distributor)))
    #    logger.info("vars(dist):        {}".format(vars(domain.dist)))
    #    logger.info('vars(domain.dist.comm):             {}'.format(dir(domain.dist.comm)))
    #    logger.info('vars(domain.dist.comm_cart):        {}'.format(dir(domain.dist.comm_cart)))
    #    logger.info('vars(domain.dist.layouts:           {}'.format(dir(domain.dist.layouts)))
    #    logger.info('vars(domain.dist.paths:             {}'.format(dir(domain.dist.paths)))
    #    logger.info('vars(domain.dist.coeff_layout:      {}'.format(dir(domain.dist.coeff_layout)))
    #    logger.info('vars(domain.dist.grid_layout:       {}'.format(dir(domain.dist.grid_layout)))
    #    logger.info('vars(domain.dist.grid_layout.blocks:      {}'.format(dir(domain.dist.grid_layout.blocks)))
    #    logger.info('vars(domain.dist.grid_layout.buffer_size: {}'.format(dir(domain.dist.grid_layout.buffer_size)))
    #    logger.info('vars(domain.dist.grid_layout.domain:      {}'.format(dir(domain.dist.grid_layout.domain)))
    #    logger.info('vars(domain.dist.grid_layout.ext_coords.shape:  {}'.format(dir(domain.dist.grid_layout.ext_coords.shape)))
    #    logger.info('vars(domain.dist.grid_layout.ext_mesh.shape:    {}'.format(dir(domain.dist.grid_layout.ext_mesh.shape)))
    #    logger.info('vars(domain.dist.grid_layout.global_shape:{}'.format(dir(domain.dist.grid_layout.global_shape)))
    #    logger.info('vars(domain.dist.grid_layout.grid_space.shape:  {}'.format(dir(domain.dist.grid_layout.grid_space.shape)))
    #    logger.info('vars(domain.dist.grid_layout.local_shape: {}'.format(dir(domain.dist.grid_layout.local_shape)))
    #    logger.info('vars(domain.dist.grid_layout.slices:      {}'.format(dir(domain.dist.grid_layout.slices)))
    #    logger.info('vars(domain.dist.grid_layout.start:       {}'.format(dir(domain.dist.grid_layout.start)))
        #        logger.info('vars(domain.dist.comm):             {}'.format(vars(domain.dist.comm)))
        #        logger.info('vars(domain.dist.comm_cart):        {}'.format(vars(domain.dist.comm_cart)))
        #        logger.info('vars(domain.dist.layouts:           {}'.format(vars(domain.dist.layouts)))
        #        logger.info('vars(domain.dist.paths:             {}'.format(vars(domain.dist.paths)))
        #        logger.info('vars(domain.dist.coeff_layout:      {}'.format(vars(domain.dist.coeff_layout)))
        #        logger.info('vars(domain.dist.grid_layout:       {}'.format(vars(domain.dist.grid_layout)))
        #        logger.info('vars(domain.dist.layout_references: {}'.format(vars(domain.dist.layout_reference)))
    
        logger.info("dir(domain):       {}".format(dir(domain)))
        logger.info("vars(distributor): {}".format(vars(domain.distributor)))
        logger.info("vars(dist):        {}".format(vars(domain.dist)))
        logger.info("variables:         {}".format(v))
        logger.info("variables:         {}".format(problem.variables))
        logger.info("substitutions:     {}".format(problem.substitutions))

    #  for j in range(domain.dim):
    #      basis=domain.bases[j]
    #      logger.info("basis: {} {}".format(j,basis.name))
    #      logger.info(grid_type(basis.grid))
    #print_vars(basis)
    #logger.info("basis.grid:")
    #print_vars(basis.grid)
    
    #if isinstance(basis,de.core.basis.Fourier):     logger.info("basis: {} {} Fourier".format(j,basis.name))
    #elif isinstance(basis,de.core.basis.SinCos):    logger.info("basis: {} {} SinCos".format(j,basis.name))
    #elif isinstance(basis,de.core.basis.Chebychev): logger.info("basis: {} {} Chebychev".format(j,basis.name))
    #logger.info("basis: {} {} Unknown  {}".format(j,basis.name,basis.__class__))
    # logger.info("vars(basis) {}".format(vars(basis)))
    gyz  = jgrid.gridyz(p,p['AA'], comm, domain.dist.grid_layout.ext_mesh[1:])
    gyzJ = jgrid.gridyz(p,p['AAJ'],comm, domain.dist.grid_layout.ext_mesh[1:])
    gx  =  jgrid.gridx(p, p['AA'], comm)
    return(domain,problem)


def boussinseq_equations(domain,problem,p,FFTAx):
    
    v=problem.variables

    ps = list(set(problem.parameters) | set(problem.substitutions))
    pv = list(set(('u','v','w','b','s')) & set(v))
    ps = comm.bcast(ps, root=0)
    pv = comm.bcast(pv, root=0)
    
    diffdom = dp.get_bool_param(p,'diffdom')
    average = dp.get_bool_param(p,'avrg')
    oddgx   = dp.get_bool_param(p,'oddg')
    Tx      = dp.get_param(p,'Tx')
    Ty      = dp.get_param(p,'Ty')
    Tz      = dp.get_param(p,'Tz')
    bSV     = dp.get_param(p,'bV')
    sSV     = dp.get_param(p,'SV')
    bVh     = dp.get_param(p,'bVh')
    sVh     = dp.get_param(p,'SVh')
    B       = dp.get_param(p,'Buoyancy')
    S       = dp.get_param(p,'Sediment')
    noise   = dp.get_param(p,'Noise')
    forcing = dp.get_param(p,'Forcing')
    force   = dp.get_param(p,'Force')
    forceb  = dp.get_param(p,'forceb')
    forces  = dp.get_param(p,'forces')
    meanu   = dp.get_param(p,'meanu')
    Skp     = dp.get_param(p,'Skp')
    Ski     = dp.get_param(p,'Ski')
    Skg     = dp.get_param(p,'Skg')
    Conservative = dp.get_bool_param(p,'Conservative')
    intdx   = dp.get_bool_param(p,'intdx')
    
    oddsx = False
    oddbx = False

    if meanu is not None: problem.parameters['meanu'] = meanu

    if 's' in v:
        if 'parity' in problem.meta['s']['x']: oddsx = problem.meta['s']['x']['parity']==-1

    if 'b' in v:
        if 'parity' in problem.meta['b']['x']: oddbx = problem.meta['b']['x']['parity']==-1

    
    if oddgx or oddsx or oddbx:
        nfx=7
        L     = domain.bases[0].interval[1]
        maxkx = domain.bases[0].wavenumbers[-1]
        maxkx = max(domain.bases[0].wavenumbers) 
        xx=domain.bases[0].grid(scale=AA).reshape((Nxx,1,1))
        problem.parameters['oddcx'] = domain.new_field()
        set_parity(problem.parameters['oddcx'],Tx,Ty,Tz,-1,1,1)
        problem.parameters['oddcx'].meta['y']['constant'] = True
        problem.parameters['oddcx'].meta['z']['constant'] = True
        problem.parameters['oddcx'].set_scales(AA)
        problem.parameters['oddcx']['g']=FFTAx.constant(xx)
        
    if p['Ny'] > 1: dim=3
    else:           dim=2
 
    if dim==3: cc=('x','y','z')
    else:      cc=('x','z')
    if dim==3: uu=('u','v','w')
    else:      uu=('u','w')
     
    if force is not None:
        problem.parameters['Fu']=force
        problem.parameters['Fv']=force
        problem.parameters['Fw']=force
        if B is not None: problem.parameters['Fb']=force*forceb
        if S is not None: problem.parameters['Fs']=force*forces

    if dim==3:      problem.substitutions['dyv']='dy(v)'
    else:           problem.substitutions['dyv']='0'
    if Tz=='Cheb':  problem.substitutions['dzw']='wz'
    else:           problem.substitutions['dzw']='dz(w)'
    if 's' in problem.variables:
        if Tz=='Cheb':  problem.substitutions['dzs']='sz'
        else:           problem.substitutions['dzs']='dz(s)'
        problem.substitutions['dxs']='dx(s)'
        if dim==3: problem.substitutions['dys']='dy(s)'
    if 'b' in problem.variables:
        if Tz=='Cheb':  problem.substitutions['dzb']='bz'
        else:           problem.substitutions['dzb']='dz(b)'
        problem.substitutions['dxb']='dx(b)'
        if dim==3: problem.substitutions['dyb']='dy(b)'

    cond1 = '(nx != 0)'
    cond2 = '(nx == 0)'
    if dim==3:    cond1 = cond1 + '  or (ny != 0)'
    if dim==3:    cond2 = cond2 + ' and (ny == 0)'
    if not Tz=='Cheb': cond1 = cond1 + '  or (nz != 0)'
    if not Tz=='Cheb': cond2 = cond2 + ' and (nz == 0)'

    if   S is not None and B is not None: problem.substitutions['rho']='(B*b+S*s)'
    elif S is     None and B is not None: problem.substitutions['rho']='B*b'
    elif S is not None and B is     None: problem.substitutions['rho']='S*s'
    elif S is     None and B is     None: problem.substitutions['rho']='0'
    
    problem.parameters['kappau']    = 1/p['Re']
    #problem.parameters['fnu']   = 0/p['Re']/1000
    if S is not None: problem.parameters['kappas'] = 1/(p['Re']*p['Scs'])
    if B is not None: problem.parameters['kappab'] = 1/(p['Re']*p['Scb'])
    if S is not None: problem.parameters['S'] = S
    if B is not None: problem.parameters['B'] = B

    g=dict()
    sV=dict()
    bV=dict()
    lhs=dict()
    rhs=dict()
    eq=dict()

    lhs['u'] = 'dt(u) - kappau*LP(u,uz) + dx(p)'
    lhs['v'] = 'dt(v) - kappau*LP(v,vz) + dy(p)'
    lhs['w'] = 'dt(w) - kappau*LP(w,wz) + dz(p)'
    lhs['b'] = 'dt(b) - kappab*LP(b,bz)'
    lhs['s'] = 'dt(s) - kappas*LP(s,sz)'

    rhs['u'] = ''
    rhs['v'] = ''
    rhs['w'] = ''
    rhs['b'] = ''
    rhs['s'] = ''

    if 'Angle' in p:
        g['x'] = math.sin(math.pi*p['Angle']/180)
        g['z'] =-math.cos(math.pi*p['Angle']/180)
        if p['Angle']==180 or p['Angle']==-180 or p['Angle']==0: del(g['x'])
        if p['Angle']== 90 or p['Angle']== -90: del(g['z'])
    #logger.info('g: {}'.format(g));

    
    if S is not None and sSV is not None:
        if  oddgx or oddsx:
            rhs['s'] += '-sSV*oddcx*dx(s)'
            problem.parameters['sVz']=sSV
        else:    
           for a in g: sV[a] = g[a]*sSV

    if B is not None and bSV is not None:
        if  oddgx or oddbx:
            if bVh is None: rhs['b'] += '-bSV*oddcx*dx(b)'
            else:           rhs['b'] += '-bSV*oddcx*dx(b*(1-bVh*b))'
            problem.parameters['bVz']=bSV
        else:    
           for a in g: bV[a] = g[a]*bSV

    #if B is not None and diffdom:  rhs['b'] += '+b*(u*ddfx+v*ddfy+w*ddfz)'
    #if diffdom:
    #    rhs['u'] += '-u*(u*ddfx+v*ddfy+w*ddfz)'#-ddfx*dx(p)
    #    rhs['v'] += '-v*(u*ddfx+v*ddfy+w*ddfz)'#-ddfy*dy(p)
    #    rhs['w'] += '-w*(u*ddfx+v*ddfy+w*ddfz)'#-ddfz*dz(p)
    # if 'ddfx' in ps: divu += '-u*ddfx'
    # if 'ddfx' in ps: rhs['u'] += '-max(u,0)**2*ddfx'    
    if bVh is not None: problem.parameters['bVh'] = bVh
    if sVh is not None: problem.parameters['sVh'] = sVh
    for a in g:  problem.parameters[ 'g'+a] =  g[a]  
    for a in sV: problem.parameters['sV'+a] = sV[a]  
    for a in bV: problem.parameters['bV'+a] = bV[a]  

    if sVh is None:
        for a in sV: lhs['s'] += '+sV'+a+'*d'+a+'s'
    else:
        for a in sV: rhs['s'] += '-sV'+a+'*d'+a+'s*(1-2*s*bVh)'
    if bVh is None:
        for a in bV: lhs['b'] += '+bV'+a+'*d'+a+'b'
    else:
        for a in bV: rhs['b'] += '-bV'+a+'*d'+a+'b*(1-2*b*bVh)'
    
    #logger.info('oddgx: {}'.format(oddgx))
    #logger.info('oddbx: {}'.format(oddbx))
    #logger.info('B:     {}'.format(B))
    #logger.info('bV:    {}'.format(bV))
    #logger.info('bSV:   {}'.format(bSV))
    #for j in problem.parameters: print_v(j,problem.parameters[j],s='parameters:')
    #quit(1)
    
    
    qsx=''
    qsy=''
    qsz=''
    if 'PIDG' in p: PIDG=p['PIDG']
    else:           PIDG=False
    if oddgx: rhs['u']+='+ g*oddcx*rho'
    else:
        for (a,a2) in zip(('u','v','w'),('x','y','z')):
            if a2 in g:
                if g[a2] !=0:
                    if PIDG: rhs[a] += '+ g*g'+a2+'*rho'
                    else:
                        lhs[a] += '- g*g'+a2+'*rho'
                        if Skg is not None and a2=='x': qsx += '+Skg*(' + 'g*g'+a2+'*rho' + ')' 
    for a in pv:
        if 'fpid'+a  in ps: lhs[a] += '- fpid'+a
        if 'rf'+a    in ps: rhs[a] += '+ rf'+a
        if len(rhs[a])>0:
            if rhs[a][0]=='+':  rhs[a]=rhs[a][1:]
 

    if Skp is not None:
        problem.parameters['Skp']=Skp
        qsx += '-Skp*dx(p)'
        qsy += '-Skp*dy(p)'
        qsz += '-Skp*dz(p)'
    if Ski is not None:
        problem.parameters['Ski']=Ski
        qsx += '+Ski*kappau*LP(u,uz)'
        qsy += '+Ski*kappau*LP(v,vz)'
        qsz += '+Ski*kappau*LP(w,wz)'
    if Skg is not None:   # Need to add in some lhs terms
        problem.parameters['Skg']=Skg
        if len(rhs['u'])>0: qsx += '+Skg*(' + rhs['u'] + ')' 
        if len(rhs['v'])>0: qsy += '+Skg*(' + rhs['v'] + ')' 
        if len(rhs['w'])>0: qsz += '+Skg*(' + rhs['w'] + ')' 

    if len(qsx)>0:
        if qsx[0]=='+':  qsx=qsx[1:]
        if qsy[0]=='+':  qsy=qsy[1:]
        if qsz[0]=='+':  qsz=qsz[1:]
        if 'u' in v and len(qsx)>1: rhs['s'] += '+dx(s*(' + qsx + '))'
        if 'v' in v and len(qsy)>1: rhs['s'] += '+dy(s*(' + qsy + '))'
        if 'w' in v and len(qsz)>1: rhs['s'] += '+dz(s*(' + qsz + '))'


    for a in pv:
        if 'w'+a     in ps:
            if 'f'+a in ps: rhs[a] += '+ F'+a+'*w'+a+'*(f'+a+'-'+a+')'  # '+Fu*wu*(fu-u)'
            else:           rhs[a] += '- F'+a+'*w'+a+'*'+a  # '-Fu*wu*u'
        if 'noise'+a in ps: rhs[a] += '+ noise'+a
        if len(rhs[a])>0:
            if rhs[a][0]=='+':  rhs[a]=rhs[a][1:]

    if Conservative:
        rhs['u'] += '-CDM(u,uz)'
        rhs['v'] += '-CDM(v,vz)'
        rhs['w'] += '-CDM(w,wz)'
    else:
        rhs['u'] += '-CD(u,uz)'
        rhs['v'] += '-CD(v,vz)'
        rhs['w'] += '-CD(w,wz)'

    if intdx: rhs['u'] += ' + d*d.integrate("x")'

    rhs['b'] += '-CD(b,bz)'
    rhs['s'] += '-CD(s,sz)'

    #    if ('meanu' in problem.parameters):
    #        if Tz=='Cheb':   problem.add_equation("integ_z(u)  = meanu", condition=cond2)
    #        else:            problem.add_equation("u  = meanu",          condition=cond2)
    #        if Tz=='Cheb':   logger.info('integ_z(u) = meanu, condition = {}'.format(cond2))
    #        else:            logger.info('u          = meanu, condition = {}'.format(cond2))
    #        conu=cond1
    #    else: conu=None
    conu=None
    if meanu==0:
        lhs['u']+='+g*gx*integ(rho)/(L*H*W)'
        #lhs[a]+='+integ(u)/(L*H*W)'
    for a in pv:
        eq[a] = lhs[a] + ' = ' + rhs[a]
        #logger.info(eq[a])

    #dt(s) - kappas*LP(s,sz) = -sSV*oddc*dx(s)
    #+dx(s*(Skp*dx(p)+Ski*kappau*LP(u,uz)+Skg*( g*oddc*rho)))
    #+dy(s*(Skp*dy(p)+Ski*kappau*LP(v,vz)))
    #+dz(s*(Skp*dz(p)+Ski*kappau*LP(w,wz)))+ Fs*ws*(fs-s)-CD(s,sz)

    # dt(s) - kappas*LP(s,sz) = Fs*ws*(fs-s) + s*(Skg*g*oddc*dx(rho)-Skp*LP(p))
    #                           dx(s)*(-u-Skp*dx(p)+kappau*Ski*LP(u,uz)+oddc*(Skg*g*rho-sSV))
    #                           dy(s)*(-v-Skp*dy(p)+kappau*Ski*LP(v,vz))
    #                           dz(s)*(-w-Skp*dz(p)+kappau*Ski*LP(w,wz))
    #logger.info("conu {}".format(conu))
    #logger.info("equ {}".format(eq['u']))
    #eq['u']="dt(u) - kappau*LP(u,uz) + dx(p) = -CDM(u,uz)"
    #eq['u']="dt(u)  =0"
    #for j in problem.substitutions:                   logger.info('substitutions: {:10s} = {:71s}'.format(j,problem.substitutions[j]))
    #logger.info(eq['u'])
    
    if conu is None:   problem.add_equation(eq['u'])
    else:
        problem.add_equation(eq['u'], condition=conu)
        logger.info('eq['u']:  {}, condition = {}'.format(eq['u'],conu))
    if 'v'  in v:      problem.add_equation(eq['v'])
    if 'w'  in v:      problem.add_equation(eq['w'])
    if 'b'  in v:      problem.add_equation(eq['b'])
    if 's'  in v:      problem.add_equation(eq['s'])
    if 'noisef2' in v: problem.add_equation("noiseT*   dt(noisef2)        + noisef2 - noise*noisef1 = 0") 
    if 'noisef1'  in v: problem.add_equation("noiseL**2*LP(noisef1,noisef1z) + noisef1                 = fnoise") 

    #if 'fpidu' in v: problem.add_equation("dt(fpidu) + fpid*fpidu - fnu*LP(fpidu,fpiduz) = fpid*wu*(fpidu+fu-u)")
    #if 'fpidv' in v: problem.add_equation("dt(fpidv) + fpid*fpidv - fnu*LP(fpidv,fpidvz) = fpid*wv*(fpidv+fv-v)")
    #if 'fpidw' in v: problem.add_equation("dt(fpidw) + fpid*fpidw - fnu*LP(fpidw,fpidwz) = fpid*ww*(fpidw+fw-w)")
    #if 'fpidb' in v: problem.add_equation("dt(fpidb) + fpid*fpidb - fnu*LP(fpidb,fpidbz) = fpid*wb*(fpidb+fb-b)")
    #if 'fpids' in v: problem.add_equation("dt(fpids) + fpid*fpids - fnu*LP(fpids,fpidsz) = fpid*ws*(fpids+fs-s)")

    divu='dx(u)'
    if 'v'  in  v: divu += ' + dy(v)'
    if 'w'  in  v: divu += ' + dzw'
    divu += ' = 0'
    if 'fd'   in ps: divu += '+fd'
    if 'ddfx' in ps: divu += '-u*ddfx'
    if 'ddfy' in ps: divu += '-v*ddfy'
    if 'ddfz' in ps: divu += '-w*ddfz'
    
    if Tz=='Cheb':
        problem.add_equation(divu)
    else:
        problem.add_equation(divu,    condition=cond1)
        problem.add_equation('p = 0', condition=cond2)
 
    for vz in v:
      if vz.endswith('z'): problem.add_equation(vz + '-dz(' + vz[:-1] + ')=0')
           
    if 'ap'  in v: problem.add_equation('dt(ap) - p = 0')  
    if 'ab'  in v: problem.add_equation('dt(ab) - b = 0')  
    if 'as'  in v: problem.add_equation('dt(as) - s = 0')  
    if 'au'  in v: problem.add_equation('dt(au) - u = 0') 
    if 'av'  in v: problem.add_equation('dt(av) - v = 0') 
    if 'aw'  in v: problem.add_equation('dt(aw) - w = 0') 
    if 'abu' in v: problem.add_equation('dt(abu) = b*u') 
    if 'abw' in v: problem.add_equation('dt(abw) = b*w') 
    if 'asu' in v: problem.add_equation('dt(asu) = s*u') 
    if 'asw' in v: problem.add_equation('dt(asw) = s*w') 
    if 'auu' in v: problem.add_equation('dt(auu) = u*u') 
    if 'avv' in v: problem.add_equation('dt(avv) = v*v') 
    if 'aww' in v: problem.add_equation('dt(aww) = w*w') 
    if 'auv' in v: problem.add_equation('dt(auv) = u*v') 
    if 'auw' in v: problem.add_equation('dt(auw) = u*w') 
    if 'avw' in v: problem.add_equation('dt(avw) = v*w')
    if 'aur' in v: problem.add_equation('dt(aur) = ur')
    if 'arr' in v: problem.add_equation('dt(arr) = ur*ur')
    if 'aru' in v: problem.add_equation('dt(aru) = u*ur')

 
def burgers_init(p):
    global xxf 
    global xx 
    global yy 
    global zz 
    global Nxx
    global Nyy
    global Nzz
    global AA
    global param
    global Kx
    global Sx
    
    L     = p['L']
    Nx    = p['Nx']
    AA    = p['AA']
    noise = p['noise']
    average = p['avrg']
    
    if p['sType']!='bg':
        logger.error('burgers_init: sType is not Burgers'.format(p['sType']))
        set_error('burgers_init: wrong type')
    
    check_error()
    param={}   
    param['U']=p['Velocity']
    param['L']=L
    param['Lx']=[0,L]
    param['Nx']=Nx
    param['Rx']=L
    param['Volume']=L
    param['NN']=Nx
 
    
    x_basis=de.Fourier('x', Nx, interval=[0,L], dealias=AA)
    x = x_basis.grid(scale=1)
    print('rank/size: {}/{}, Nx: {}, L: {}, AA: {}, x:[{},{}]'.format(mpirank,mpisize,Nx,L,AA,x[0],x[-1]))
    domain=de.Domain([x_basis], grid_dtype=np.float64,mesh=(mpisize,))
    dx=np.mean(np.diff(np.asarray(x_basis.grid(scale=1))))
    xx = x_basis.grid(scale=AA)
    Nxx=len(xx)
    xxf=xx
    param['NNN']=Nxx


    v=['u']
    
    if average: v = v + ['au']
    if noise>0:
        v = v + ['noisef1']        
        v = v + ['noisef2']        
    if p['fpidu']: v = v + ['fpidu']
    problem = de.IVP(domain, variables=v, time='t')
    v=problem.variables

    problem.substitutions['SR']   = 'dx(u)**2'
 
    #logger.info("variables:    {}".format(problem.variables))
    #logger.info("variables: {}".format(v))
    logger.info("AA:  {:7.2f}, AB:    {:7.2f}".format(AA,p['AB']))
    logger.info("L:   {:7.2f}, dx:    {:8.5f}, Nx:    {:7d}, ".format(L,dx,Nx))
    logger.info("Grid cells                 {:8.3f} {:8.3f} million".format(param['NN']/1e6,param['NNN']/1e6))
    logger.info("Grid cells/core            {:8.3f} {:8.3f} thousand".format(np.float(param['NN'])/np.float(mpisize)/1e3,np.float(param['NNN'])/np.float(mpisize)/1e3))
    if VMEM is not None:logger.info("VMEM:           {:6.1f}GB".format(VMEM))
    append_file(fnlog,'domain.dist.rank:     {}\n'.format(domain.dist.rank))   
    append_file(fnlog,'domain.dist.size:     {}\n'.format(domain.dist.size))
    append_file(fnlog,'domain.dist.coords:   {}\n'.format(domain.dist.coords)) 
    append_file(fnlog,'domain.dist.mesh:     {}\n'.format(domain.dist.mesh))  
    append_file(fnlog,'domain.grid(0).shape: {}\n'.format(domain.grid(0).shape))
    append_file(fnlog,'                   \n')
    append_file(fnlog,'domain.dist.grid_layout\n')
    append_file(fnlog,'  ext_coords     {}\n'.format(domain.dist.grid_layout.ext_coords))  
    append_file(fnlog,'  ext_mesh       {}\n'.format(domain.dist.grid_layout.ext_mesh))    
    append_file(fnlog,'  index          {}\n'.format(domain.dist.grid_layout.index))       
    append_file(fnlog,'  local          {}\n'.format(domain.dist.grid_layout.local))       
    append_file(fnlog,'  grid_space     {}\n'.format(domain.dist.grid_layout.grid_space))  
    append_file(fnlog,'                   \n')
    append_file(fnlog,'  scales         {}\n'.format(1))
    append_file(fnlog,'    global_shape {}\n'.format(domain.dist.grid_layout.global_shape(scales=1))) 
    append_file(fnlog,'    local_shape  {}\n'.format(domain.dist.grid_layout.local_shape(scales=1))) 
    append_file(fnlog,'    slices       {}\n'.format(domain.dist.grid_layout.slices(scales=1)))         
    append_file(fnlog,'    start        {}\n'.format(domain.dist.grid_layout.start(scales=1)))           
    append_file(fnlog,'                   \n')
    append_file(fnlog,'  scales         {}\n'.format(AA))
    append_file(fnlog,'    global_shape {}\n'.format(domain.dist.grid_layout.global_shape(scales=AA))) 
    append_file(fnlog,'    local_shape  {}\n'.format(domain.dist.grid_layout.local_shape(scales=AA))) 
    append_file(fnlog,'    slices       {}\n'.format(domain.dist.grid_layout.slices(scales=AA)))         
    append_file(fnlog,'    start        {}\n'.format(domain.dist.grid_layout.start(scales=AA)))           
    append_file(fnlog,'                   \n')
    append_file(fnlog,'domain.dim:      {}\n'.format(domain.dim))
    return(domain,problem)

def burgers_equations(domain,problem,p):
    
    if p['sType']!='bg':
        logger.error('burgers_equations: sType is not Burgers'.format(p['sType']))
        set_error('burgers_equations: wrong type')
    check_error()
    
    v=problem.variables
    
    problem.parameters['Fu']=p['Force']
    
    #    problem.add_equation("dt(u) - nu*dx(dx(u))  - fpidu = g + noiseu + Fu*wu*(fu-u) - u*dx(u)")
    problem.add_equation("dt(u) - kappau*dx(dx(u))  - fpidu = g + noiseu + Fu*wu*(fu-u) - dx(u**2/2)")
    if 'fpidu'  in v: problem.add_equation("dt(fpidu) + fpid*fpidu - fnu*dx(dx(fpidu)) = fpid*wu*(fpidu+fu-u)")
    if 'au'     in v: problem.add_equation('dt(au)     - u     = 0')
    if 'afpidu' in v: problem.add_equation('dt(afpidu) - fpidu = 0') 
   
def print_problem(p):
    for j in range(len(p.equations)):           logger.info('equation  {:2d}: {}'.format(j, p.equations[j]['raw_equation']))
    if 'boundary_conditions' in dir(p):
        for j in range(len(p.boundary_conditions)): logger.info('boundary  {:2d}: {}'.format(j, p.boundary_conditions[j]['raw_equation']))
    for j in p.substitutions:                   logger.info('substitutions: {:10s} = {:71s}'.format(j,p.substitutions[j]))
    for j in p.parameters:                      print_v(j,p.parameters[j],s='parameters:')
    
def cubic(x,d):  # cubic spline centre on zero with step size d 
    if d<=0:
        f=(1+np.sign(x))/2
 #       g=f
    else:
        z=1-abs(x/d)
        f=(1+np.sign(x)*np.minimum(1,1-z*z*z))/2
    return(f)

def cubich(y,x1,x2): # cubic hump interpolating between 0 at x1 and x2 and 1 at (x1+x2)/2
    if x1==x2:
        f=0*y
    else:
        x=2*((2*y-x1-x2)/(x2-x1))**2
        z=2-x
        f=np.maximum(0,np.minimum(1-3*x**2/2*(1+x/2),z**3/4));

    return(f)

def cubici(y,x1,x2):  # cubic spline interpolating between 0 and 1 between x1 and x2
    if x1==x2:
        f = (1-np.sign(y-x1))/2
    else:
        x=(2*y-x1-x2)/(x2-x1)
        z=np.maximum(0,(-1+2*abs(x)))
        f=0.5+np.sign(x)*np.minimum(3, 2*abs(x)*(-2*x*x+3)+z*z*z)/6
    return(f)

def dcubici(y,x1,x2): # derivative of cubic spline interpolating between 0 and 1 between x1 and x2
    x=2*abs((2*y-x1-x2)/(x2-x1))
    z1=np.maximum(0,1-x)
    z2=np.maximum(0,2-x)
    f=(z2*z2-2*z1*z1)/(x2-x1);
    return(f)

def handler(signum, frame):
    debug =  False
    
    global dump 
    if signum == signal.SIGTERM:
        if debug: print("SIGTERM {:d}".format(mpirank))
        set_abort('SIGTERM')
    elif signum == signal.SIGINT:    
        if debug: print("SIGINT {:d}".format(mpirank))
        set_abort('SIGINT')
    elif signum == signal.SIGHUP:    
        if debug: print("SIGHUP {:d}".format(mpirank))
        dump = True
    elif signum == signal.SIGUSR1: 
        if debug: print("SIGUSR1 {:d}".format(mpirank))
    elif signum == signal.SIGUSR2: 
         if debug: print("SIGUSR2 {:d}".format(mpirank))
    elif signum == signal.SIGQUIT: 
        if debug: print("SIGQUIT {:d}".format(mpirank))
        set_abort('SIGQUIT')
    elif signum == signal.SIGALRM: 
        if debug: print("SIGALRM {:d}".format(mpirank))
        dump = True
    elif signum == signal.SIGVTALRM:
        if debug: print("VTSIGALRM {:d}".format(mpirank))
        dump = True
    else:
        if debug: print('Signal handler called with signal {:d}'.format(mpiranksignum))


def init_signals():
    signal.signal(signal.SIGTERM,  handler)
    signal.signal(signal.SIGINT,   handler)  
    signal.signal(signal.SIGQUIT,  handler)
    signal.signal(signal.SIGHUP,   handler)  
    signal.signal(signal.SIGUSR1,  handler)
    signal.signal(signal.SIGUSR2,  handler)
    signal.signal(signal.SIGALRM,  handler)
    signal.signal(signal.SIGVTALRM,handler) 


# Convert seconds to dy:hr:min:sec string         
def sec2dhms(x):
    y=np.int(x/365/3600/24)
    x=x-24*3600*365*y
    d=np.int(x/3600/24)
    x=x-24*3600*d
    h=np.int(x/3600)
    x=x-3600*h
    m=np.int(x/60)
    s=int(x-60*m)

    r=''
    
    if y>0:               r=r+'{:04d}-'.format(1970+y)
    if y>0 or d>0:        r=r+'{:02d} '.format(d)
    if y>0 or d>0 or h>0: r=r+'{:02d}:'.format(h)
    r=r+'{:02d}:'.format(m)
    r=r+'{:02d}'.format(s)
    
    return(r)

def run_info():
    logger.info("MPI rank: {:d}, size: {:d}".format(mpirank,mpisize))
    if JOB_ID and JOB_QUEUE: logger.info("JOB_ID: {:s}, JOB_QUEUE: {:s}".format(JOB_ID,JOB_QUEUE))
    logger.info("NODENAME: {:s}, PID: {:s}".format(NODENAME,PID))

# Check if an abort file has been created or an abort signal has been received
# propagate this abort signal to all nodes
def check_abort_file():
    global abort
    if mpirank==0:
        afn=basedir / "abort"
        if afn.is_file():
            set_abort('abort file')
            os.remove(str(afn))
            logger.info('Abort file detected')
    
def reduce_bool(b):
    if  b: a=np.int(1)
    else : a=np.int(0)
    a=comm.allreduce(a, op=MPI.SUM)
    if a>0: b = True
    else:   b = False
    #if a>0 : logger.info("reduce_bool: received by {:d} tasks".format(a))             
    return(b)

def set_abort(s,t=None):
    global abort
    global abort_s
    global terminate
    if t is not None: terminate=t
    abort=True
    abort_s=s

def set_error(s):
    global jerror
    global jerror_s
    jerror=True
    jerror_s=s

def test_abort():
    global abort
    global abort_s
    abort = reduce_bool(abort)
    if abort : logger.info("Abort received: {}".format(abort_s))             
    return(abort)

def test_error():
    global jerror
    global jerror_s
    jerror = reduce_bool(jerror)
    if jerror: logger.info("Error received: {}".format(jerror_s))             
    return(jerror)
    
def check_abort():
    global abort
    global jerror
    test_abort()
    test_error()
    if abort or jerror:
        logger.info("check_abort: Aborting")
        if jerror: write_status("Exception")
        else:      write_status("Aborted")
        delete_lock()
        #if mpisize>1: comm.Abort()
        exit(1)

def check_error(s=None):
    global jerror
    if s is None: s=''
    else: jerror=True
    test_error()
    if jerror:
        logger.info("check_error: Aborting"+s)
        delete_lock()
        #if mpisize>1: comm.Abort()
        exit(1)
      
def _getmtime(entry):
    return entry.stat().st_mtime
   
def find_start_files(d):
    debug = False
    fns=[]
    if mpirank==0:
        if debug: logger.info('find_start_file d={}'.format(d))

        d=Path(d)
             
        if d.is_dir():
            ckfn=d / "final"
            files = sorted(ckfn.glob("state*.hdf5"))
            for f in files:
                if f.is_file(): fns.append(f)
            G = [d / "jfinal.hdf5",
                 d / "final/jfinal.hdf5",
                 d / "final_s1.h5",       
                 d / "final_s1.hdf5",     
                 d / "final/final_s1.h5", 
                 d / "final/final_s1.hdf5"]
            for g in G:
                #logger.info('Trying {}'.format(g))
                if g.is_file(): fns.append(g)

            ckfn=d / "checkpoint"
            files = sorted(ckfn.glob("checkpoint_s*.h*5"))
            for f in files:
                if f.is_file(): fns.append(f)
        if len(fns)==0:
            logger.info('No starting files found with {}'.format(d))
            set_error('find_start_files: No starting files found')
        else: fns=sorted(fns,key=_getmtime, reverse=True)
        #fns.sort(key=os.path.getmtime,reverse=True)        
        #logger.info('find_start_files: fns {}'.format(fns))
    check_error()
    fns = comm.bcast(fns, root=0)
    return(fns)


def find_param_file(d1,nm,fn):
    debug = False
    f=None
    if mpirank==0:
        found = False
        d1 = Path(d1)
        nm = Path(nm)
        if debug: logger.info('d1={}'.format(d1))
        if debug: logger.info('nm={}'.format(nm))
        if debug: logger.info('fn={}'.format(fn))

        f1 = d1
        f3 = d1 / nm
        if fn is not None:
            fn = Path(fn)
            if debug: logger.info('fn={}'.format(fn))
            f2 = d1 / fn
            f4 = d1 / nm / fn
            F  = [f1, f1 / '.h5', f1 / '.hdf5', f1 / 'param.h5', f1 / 'param.hdf5',
                  f2, f2 / '.h5', f2 / '.hdf5', f2 / 'param.h5', f2 / 'param.hdf5',
                  f3, f3 / '.h5', f3 / '.hdf5', f3 / 'param.h5', f3 / 'param.hdf5',
                  f4, f4 / '.h5', f4 / '.hdf5', f4 / 'param.h5', f4 / 'param.hdf5']
        else:
            F  = [f1, f1 / '.h5', f1 / '.hdf5', f1 / 'param.h5', f1 / 'param.hdf5',
                  f3, f3 / '.h5', f3 / '.hdf5', f3 / 'param.h5', f3 / 'param.hdf5']
 
        if debug: logger.info('F={}'.format(F));

        for f in F:
            if debug: logger.info('Checking for {}'.format(f))
            if f.is_file():
                found=True
                break
        if not found:
            logger.info('No starting parameter file found with "{}" "{}" "{}"'.format(d1,nm,fn))
            set_error('find_param_file: none')
        if debug: logger.info('find_param_file={}'.format(f))
    check_error()
    return(f)
        
def get_oldest_file(files):
    """ Find and return the oldest file of input file names.
    Only one wins tie.   """
    if not files:  return(None) # Check for empty list.
    now = time.time()           # Raw epoch distance.
    oldest = files[0]
    old=oldest.stat().st_ctime-now
    # Iterate over all remaining files.
    for f in files[1:]:
        #        age = now - os.path.getctime(f)
        age = f.stat().st_ctime-now
        if age<old:
            old=age
            oldest=f
    return(oldest)
    
def get_youngest_file(files):
    """ Find and return the oldest file of input file names.
    Only one wins tie.   """

    if not files:  return(None)# Check for empty list.
    now = time.time()# Raw epoch distance.
    newest = files[0]
    new=newest.stat().st_ctime-now
    # Iterate over all remaining files.
    for f in files[1:]:
        age = f.stat().st_ctime-now
        if age>new:
            new=age
            newest=f
    return(newest)
    
def sort_files_by_age(files):
    """ Find and return the files in order of newest to oldest  """
    files.sort(key=os.path.getmtime)
    return(files)
    

def end_stats(timea,sim_time_start,solver,domain):
    global param
    if mpirank==0:
        NN  = param['NN']
        NNN = param['NNN']
        
        timea.append(np.float(time.time()))
        sim_time=solver.sim_time-sim_time_start
        dt10 = timea[1]-timea[0]
        dt21 = timea[2]-timea[1]
        n_iter_loop = solver.iteration
        logger.info('Iterations:                {}'.format(n_iter_loop))
        logger.info("Grid cells:                {:10.3f} {:10.3f} million".format(param['NN']/1e6,param['NNN']/1e6))
        logger.info("Grid cells/core:           {:10.3f} {:10.3f} thousand".format(np.float(param['NN'])/np.float(mpisize)/1e3,np.float(param['NNN'])/np.float(mpisize)/1e3))
        logger.info('Sim time total:            {:6.2f}, from {:6.2f} to {:6.2f}'.format(sim_time,sim_time_start,solver.sim_time))
        logger.info('Wall time init:            {}'.format(sec2dhms(dt10)))
        logger.info('Wall time loop:            {}'.format(sec2dhms(dt21)))
        logger.info('Efficiency: it/sec         {:10.4f} it/s nodes {:4.0f}'.format(n_iter_loop/dt21,mpisize))
        logger.info('Efficiency: it/day         {:10.1f} it/d nodes {:4.0f}'.format(n_iter_loop/dt21*24*60*60,mpisize))
        logger.info('Efficiency: it/day/node    {:10.6f} it/d nodes {:4.0f}'.format(n_iter_loop/dt21*24*60*60/mpisize,mpisize))
        if  n_iter_loop>0: logger.info('Seconds per iteration      {:10.4f} s/it'.format(dt21/n_iter_loop))
        if sim_time>0:
            logger.info('Efficiency: wall/sim       {:10.1f} s/sim'.format(dt21/sim_time))
            logger.info('Efficiency: cpu*wall/sim   {:10.1f} s/sim'.format(dt21/sim_time*mpisize))
        logger.info('Efficiency: sim/wall       {:10.8f} sim/s {:10.2f} sim/d {:d} cpus'.format(sim_time/(dt21)        ,sim_time/dt21*24*60*60,        mpisize))
        logger.info('Efficiency: sim/wall/cpu   {:10.8f} sim/s {:10.2f} sim/d {:d} cpus'.format(sim_time/(dt21*mpisize),sim_time/dt21*24*60*60/mpisize,mpisize))
        logger.info('Efficiency: N*sim/wall     {:10.3f} gp*sim/s'.format(NN*sim_time/dt21))
        logger.info('Efficiency: N*sim/cpu/wall {:10.3f} gp*sim/s/cpu, {:d} cpus'.format(NN*sim_time/(dt21*mpisize),mpisize))

#Print information about the total runtime
#Write the status file
#And delete the lock file
#Delete the file with the nodename and pid for each task
def finish(time1,solver):
    global fnnode
    global terminate
    try: fnnode.unlink()
    except: print('Node file {} on {:04d}/{:04d} cannot be deleted'.format(fnnode,mpirank,mpisize))
    delete_lock()
    dt = time.time()-time1
    terminate = terminate or (solver.sim_time >= solver.stop_sim_time)
    terminate=reduce_bool(terminate)
    logger.info('Wall time: total           {}'.format(sec2dhms(dt)))
    if terminate: write_status("Terminated")
    else:         write_status("Aborted")
          
def jint(f):
    if f is None: return(np.float(0))
    global Nx
    #if 'scale' in f.meta['x']: sc=f.meta['x']['scale']
    sc=f['g'].shape[0]/Nx
    if 'y' in f.meta:
        x = f.integrate('x','y','z')['g']
        if mpirank==0: x=np.float(x[0,0,0])
        else:          x=np.float(0)
    else:
        x = f.integrate('x','z')['g']
        if mpirank==0: x=np.float(x[0,0])
        else:          x=np.float(0)
    f.set_scales(sc)
    x=comm.bcast(x, root=0)
    return(x)

def get_proc_status(keys = None):
    with open('/proc/self/status') as f:
        data = dict(map(str.strip, line.split(':', 1)) for line in f)
    return tuple(data[k] for k in keys) if keys else data

def MemAll():
    a=('VmPeak','VmSize','VmLck','VmPin','VmHWM','VmRSS','VmData','VmStk','VmExe','VmLib','VmPTE','VmSwap')
    x=np.zeros(len(a))
    if LINUX:
        b=get_proc_status(a)
        for k in range(len(a)): x[k]=np.float64(b[k][0:-3])/1024**2
        
    d=process.memory_info()
    c=d._fields
    y=np.zeros(len(c))
    for k in range(len(c)): y[k]=np.float64(d[k])/1024**3
    
    x = nodecomm.allreduce(x, op=MPI.SUM)/nodenum
    y = nodecomm.allreduce(y, op=MPI.SUM)/nodenum
         
    if mpirank==0:
        if LINUX:
            for k in range(len(a)): logger.info('Proc:    {:6s} {:11.4f} GB'.format(a[k],x[k]))
        for k in range(len(c)): logger.info('Process: {:6s} {:11.4f} GB'.format(c[k],y[k]))
     
def MemTrack():
    global mem
    global process
    global nodenum     # Number of nodes
    global nodeid      # which node we are on
    global noderank    # Which process we are on this node
    global nodesize    # How many processes on this node
    global nodecomm    # communicator for processes on this node
    global rootcomm
    global RSSPNODE
    global VMSPNODE

    if OSX:
        a=process.memory_info()
        RSS = np.int64(a.rss)
        VMS = np.int64(a.vms)
    else:
        RSS,VMS = get_proc_status(('VmRSS','VmSize'))
        RSS  = 1024*np.int64(RSS[0:-3])
        VMS  = 1024*np.int64(VMS[0:-3])

    RSSPNODE = comm.allreduce(RSS, op=MPI.SUM)/nodenum
    VMSPNODE = comm.allreduce(VMS, op=MPI.SUM)/nodenum

    mem['RSS']     = nodecomm.reduce(RSS, op=MPI.SUM)
    mem['VMS']     = nodecomm.reduce(VMS, op=MPI.SUM)

    if noderank==0:
        mem['RSSMax']  = rootcomm.reduce(mem['RSS'], op=MPI.MAX)
        mem['VMSMax']  = rootcomm.reduce(mem['VMS'], op=MPI.MAX)
        mem['RSSMin']  = rootcomm.reduce(mem['RSS'], op=MPI.MIN)
        mem['VMSMin']  = rootcomm.reduce(mem['VMS'], op=MPI.MIN)
        mem['RSSMean'] = rootcomm.reduce(mem['RSS'], op=MPI.SUM)
        mem['VMSMean'] = rootcomm.reduce(mem['VMS'], op=MPI.SUM)

    if mpirank==0:
        mem['RSSHwm']  = max(mem['RSSHwm'],mem['RSSMax'])
        mem['VMSHwm']  = max(mem['VMSHwm'],mem['VMSMax'])
        mem['RSSMean'] = mem['RSSMean']/nodenum
        mem['VMSMean'] = mem['VMSMean']/nodenum

        
def MemPrint():
    global mem
    if mpirank==0: logger.info('MemTrack: VMS: [{} {} {} {} {}]'.format(MemStr(mem['VMSMin']),MemStr(mem['VMSMean']),MemStr(mem['VMSMax']),MemStr(mem['VMSHwm']),MemStr(mem['PHYS'])))

def MemShort():
    global mem
    global RSSPNODE
    global VMSPNODE
    #if mpirank==0: return(MemStr(mem['VMSMax']))
    #else: return('')
    return(MemStr(VMSPNODE))
    
def MemStats():
    global VmHWM
    global VmPeak
    MaxRSS=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if LINUX: MaxRSS=MaxRSS*1024
    if OSX:
        VmRSS  = np.int64(process.memory_info().rss)
        VmSize = np.int64(process.memory_info().vms)
        if VmRSS  > VmHWM:  VmHWM=VmRSS
        if VmSize > VmPeak: VmPeak=VmSize
    else:
        VmHWM,VmRSS,VmPeak,VmSize = get_proc_status(('VmHWM', 'VmRSS','VmPeak', 'VmSize'))
        VmHWM  = 1024*np.int64(VmHWM[0:-3])
        VmRSS  = 1024*np.int64(VmRSS[0:-3])
        VmPeak = 1024*np.int64(VmPeak[0:-3])
        VmSize = 1024*np.int64(VmSize[0:-3])
    MaxRSS    = MemStr(comm.allreduce(MaxRSS, op=MPI.MAX))
    HwmRSS    = MemStr(comm.allreduce(VmHWM,  op=MPI.MAX))
    MaxPeakVM = MemStr(comm.allreduce(VmPeak, op=MPI.MAX))
    SumRSS    = MemStr(comm.allreduce(VmRSS,  op=MPI.SUM))
    SumVM     = MemStr(comm.allreduce(VmSize, op=MPI.SUM))
    MaxVM     = MemStr(comm.allreduce(VmSize, op=MPI.MAX))
    #s='MaxRSS: {}, HwmRSS: {}, MaxVM: {}, SumRSS: {}, SumVM: {}'.format(MaxRSS,HwmRSS,MaxPeakVM,SumRSS,SumVM)
    s='MaxPeakVM: {}, MaxVM: {}, SumVM: {}'.format(MaxPeakVM,MaxVM,SumVM)
    return(s)

def MemoryUsage():
    ''' Memory usage in bytes '''
    global mm
    global fnmem
    if OSX:
        process = psutil.Process(pid)
        m=np.int64(process.memory_info().rss)
        if m>mm:
            mm=m
            write_file(fnmem,MemStr(mm) + '\n')
            append_file(fnlog,MemStr(mm) + '\n')
    else:
        VmHWM,VmRSS,VmPeak,VmSize,VmData = get_proc_status(('VmHWM', 'VmRSS','VmPeak', 'VmSize','VmData'))
        m=1024*np.int64(VmRSS[0:-3])
        VmHWM  = MemStr(1024*np.int64(VmHWM[0:-3]))
        VmRSS  = MemStr(m)
        VmPeak = MemStr(1024*np.int64(VmPeak[0:-3]))
        VmSize = MemStr(1024*np.int64(VmSize[0:-3]))
        VmData = MemStr(1024*np.int64(VmData[0:-3]))
        write_file(fnmem,'VmHWM:{:>15}\nVmRSS:{:>15}\nVmPeak:{:>14}\nVmSize:{:>14}\nVmData:{:>14}\n'.format(VmHWM,VmRSS,VmPeak,VmSize,VmData))
        
    return(m)

def MemStr(m):
    m=np.float(m)
    if   m<1024**1: s='{:7.0f} B'.format(m)
    elif m<1024**2: s='{:7.3f} kB'.format(m/1024)
    elif m<1024**3: s='{:7.3f} MB'.format(m/1024**2)
    else:           s='{:7.3f} GB'.format(m/1024**3)
    return(s)

def MaxMem():
    x=MemoryUsage()
    x=comm.allreduce(x, op=MPI.MAX)
    return(x)


def MPImax(x,d):
    m=x.max(d);
    m=comm.allreduce(m, op=MPI.MAX)
    return(m)

def MPImin(x,d):
    m=x.min(d);
    m=comm.allreduce(m, op=MPI.MIN)
    return(m)

def MemInfo():
    m     = MemoryUsage()
    minm  = MemStr(comm.allreduce(m, op=MPI.MIN))
    maxm  = MemStr(comm.allreduce(m, op=MPI.MAX))
    summ  = comm.allreduce(m, op=MPI.SUM)
    meanm = MemStr(summ/mpisize)
    summ  = MemStr(summ)
    logger.info("Memory min: {}, max: {}, mean: {}, total: {}".format(minm,maxm,meanm,summ))

def factors(n):
    gaps = [1,2,2,4,2,4,2,4,6,2,6]
    length, cycle = 11, 3
    f, fs, next = 2, [], 0
    while f * f <= n:
        while n % f == 0:
            fs.append(f)
            n /= f
        f += gaps[next]
        next += 1
        if next == length:
            next = cycle
    if n > 1: fs.append(n)
    return fs
 
def factor2(n):
    ''' Return all pairs of integers that multiply to give n '''
    f1=[]
    f2=[]
    f=np.int(1)
    while f*f <=n:
      if n % f == 0:
          f1.append(f)
          f2.append(np.int(n/f))
      f +=1
    f3=f1
    f1=f1+f2
    f2=f2+f3
    return(f1,f2)

def factor3(n):
    ''' Return all triples of integers that multiply to give n '''
    f1=[]
    f2=[]
    f3=[]
    fa=np.int(1)
    while fa*fa <=n:
        if n % fa == 0:
            m=np.int(n/fa)
            fb=np.int(1)
            while fb*fb <=m:
                if m % fb == 0:             
                    f1.append(fa)
                    f2.append(fb)
                    f3.append(np.int(m/fb))
                fb += 1
        fa += 1
    fa = f1+f1+f2+f2+f3+f3
    fb = f2+f3+f1+f3+f1+f2
    fc=  f3+f2+f3+f1+f2+f1
    return(fa,fb,fc)

def meshxyz(n,Nx,Ny,Nz):
    mx=1;
    my=1;
    mz=1;
    if Nz==0:
        mx,my=meshxy(n,Nx,Ny)
    elif n <= min(16,Nx):
        mx = n
    elif n <= min(256,Nx*Ny):
        f1, f2 = factor2(n)
        fa=np.array(f1,dtype=np.float)/Nx
        fb=np.array(f2,dtype=np.float)/Ny
        j = np.argmin((fa-fb)**2)
        mx=f1[j]
        my=f2[j]
    else: 
        f1, f2, f3 = factor3(n)
        fa=np.array(f1,dtype=np.float)/Nx
        fb=np.array(f2,dtype=np.float)/Ny
        fc=np.array(f3,dtype=np.float)/Nz
        j = np.argmin((fa-fb)**2+(fa-fc)**2+(fb-fc)**2)
        mx=f1[j]
        my=f2[j]
        mz=f3[j]
    if mx*my*mz != n:
        logger.error("meshsz: mx*my*mx !=n")
        logger.error('n={:4d} {:4d}, mx={:4d}, my={:4d}, mz={:4d}'.format(i,mx*my*mz,mx,my,mz))
        set_error('meshxyz')

    return(mx,my,mz)

def meshxy(n,N1,N2):
    mx=1;
    my=1;
    
    if n <= min(16,N1):
        mx = n
    else:
        import multiprocessing
        nc=multiprocessing.cpu_count()
        if N1>N2:
            if nc*np.int(n/nc)==n: return(np.int(nc),np.int(n/nc))
        else:
            if nc*np.int(n/nc)==n: return(np.int(n/nc),np.int(nc))
            
        f1, f2 = factor2(n)
        fa=np.array(f1,dtype=np.float)/N1
        fb=np.array(f2,dtype=np.float)/N2
        j = np.argmin((fa-fb)**2)
        mx=f1[j]
        my=f2[j]
    if mx*my != n:
        logger.error("meshxy: mx*my !=n")
        logger.error('n={:4d} {:4d}, mx={:4d}, my={:4d}'.format(i,mx*my,mx,my))
        set_error('meshxy')

    return(mx,my)

def meshxy2(n,Nx,Ny):
    mx=1;
    my=1;
    
    if n <= min(16,Nx):
        mx = n
    else:
        N=Nx*Ny/n
        
        f1, f2 = factor2(n)
        fa=np.array(f1,dtype=np.float)/Nx
        fb=np.array(f2,dtype=np.float)/Ny
        j = np.argmin((fa-fb)**2)
        mx=f1[j]
        my=f2[j]
    if mx*my != n:
        logger.error("meshxy: mx*my !=n")
        logger.error('n={:4d} {:4d}, mx={:4d}, my={:4d}'.format(i,mx*my,mx,my))
        set_error('meshxy')

    return(mx,my)

#import jutil_131 as ju
#ju.meshsz_test(1,16)
#ju.meshsz_test(16,64)
def meshsz_test(m,n):
    i=m
    while i<=n:
        mx,my,mz=meshsz(i, 800,500,500)
        logger.info('n={:4d} {:4d}, mx={:4d}, my={:4d}, mz={:4d}'.format(i,mx*my*mz,mx,my,mz))
        i += 1

        
def write_hdf5(fn,aa,xx):
    if mpirank==0:
        p=basedir / fn
        if not p.parent.is_dir(): p.parent.mkdir()
        f=h5py.File(str(p), mode='w')
        if type(aa)==list or type(aa)==tuple:
            for (a,x) in zip(aa, xx): f[a]=x
        elif type(aa)==dict:
            for     a in aa.keys():   f[a]=x[a]
        else:
            f[aa]=xx
        f.close()

def write_hdf5_vars(fn,aa,v):
    if mpirank==0:
        p=basedir / fn
        if not p.parent.is_dir(): p.parent.mkdir()
        with h5py.File(str(p), mode='w') as f:
            for a in set(aa) & set([*v]):
                #logger.info('{} {}'.format(a,v[a]))
                if   isinstance(v[a],slice):      f[a]=np.asarray(range(v[a].start,v[a].stop))
                elif isinstance(v[a],np.ndarray): f[a]=v[a].flatten()
                elif v[a] is not None:            f[a]=v[a]
    
def restartf(rfn, solver, index):
    dt=np.float(0)
    if not rfn.is_file():
        logger.info('Restart file "{}" does not exist'.format(f))
        return(dt)
    logger.info("Loading solver state from: {}".format(rfn))
    with h5py.File(str(rfn), mode='r') as file:
    # Load solver attributes
        try:    write = file['scales']['write_number'][index]
        except: write=0
        try:             dt = file['scales']['timestep'][index]
        except KeyError: dt = None
        solver.iteration = solver.initial_iteration = file['scales']['iteration'][index]
        solver.sim_time  = solver.initial_sim_time = file['scales']['sim_time'][index]
        # Log restart info
        logger.info("Loading iteration: {}".format(solver.iteration))
        logger.info("Loading write:     {}".format(write))
        logger.info("Loading sim time:  {:8.3f}".format(solver.sim_time))
        logger.info("Loading timestep:  {:8.5f}".format(dt))
        # Load fields
        for field in solver.state.fields:
            #logger.info('restart: Looking for field "{}" in restart file'.format(field.name))
            try:
                dset = file['tasks'][field.name]
                # Find matching layout
                for layout in solver.domain.dist.layouts:
                    if np.allclose(layout.grid_space, dset.attrs['grid_space']):
                        break
                    else:
                        raise ValueError("No matching layout")
                # Set scales to match saved data
                scales = dset.shape[1:] / layout.global_shape(scales=1)
                scales[~layout.grid_space] = 1
                # Extract local data from global dset
                dset_slices = (index,) + layout.slices(tuple(scales))
                local_dset = dset[dset_slices]
                # Copy to field
                field_slices = tuple(slice(n) for n in local_dset.shape)
                field.set_scales(scales, keep_data=False)
                field[layout][field_slices] = local_dset
                field.set_scales(solver.domain.dealias, keep_data=True)
            except KeyError:
                logger.info('restart: field "{}" does not exist'.format(field.name))
                
    return(dt)


def calc_maxe(flg,ss,N,X,f):
    if X in ss.field_names:
        if  flg:
            x=ss[X]
            x.set_scales(AA)
            x=x['g']
            if np.ndim(x)==2: x=abs(x[N,:  ]-f)
            else:             x=abs(x[N,:,:]-f)
            if x.size==0: x=0
            else:          x=x.max()
        else:
            x=0
        x=comm.allreduce(x, op=MPI.MAX)
        se=', '+X+'e:{:6.4f} '.format(x)
    else:
        x=None
        se=''
    return(x,se)

def calc_uerror(problem,solver,domain):
    global AA
    if 'EFN' in problem.parameters:
        N=problem.parameters['EFN']
        s=domain.dist.grid_layout.slices(scales=AA)[0]
        flg=N>=s.start and N<s.stop
        if 'EFB' in problem.parameters: be,bs=calc_maxe(flg,solver.state,N,'b',problem.parameters['EFB'])
        if 'EFS' in problem.parameters: se,ss=calc_maxe(flg,solver.state,N,'s',problem.parameters['EFS'])
        if 'EFU' in problem.parameters: ue,us=calc_maxe(flg,solver.state,N,'u',problem.parameters['EFU'])
        if 'EFV' in problem.parameters: ve,vs=calc_maxe(flg,solver.state,N,'v',problem.parameters['EFV'])
        if 'EFW' in problem.parameters: we,ws=calc_maxe(flg,solver.state,N,'w',problem.parameters['EFW'])
        s=bs+ss+us+vs+ws
    else:
        s=''
        be=None
        se=None
        ue=None
        ve=None
        we=None
    #logger.info('s:"{}"'.format(s))
    #logger.info('bs:"{}" ss:"{}" us:"{}" vs:"{}" ws:"{}"'.format(bs,ss,us,vs,ws))
    #logger.info('be {} se {} ue {} ve {} we {} '.format(be,se,ue,ve,we))
        
    return(s,be,se,ue,ve,we)
 
def calc_noise(b):
    global Nxx
    global Nyy
    global Nzz
    bb=b['g']**2
    N=bb.sum(axis=None) # Sum all the elements
    N=comm.allreduce(N, op=MPI.SUM)
    N=np.sqrt(N/(Nxx*Nyy*Nzz))
    return(N)

def calc_cm(b):
    global param
    X=np.nan
    if param['Ny']>1: c=b.integrate('y','z')
    else:             c=b.integrate('z')
    c=c['c']
    if mpirank==0:
        if param['Ny']>1: q=np.angle(c[1,0,0])
        else:             q=np.angle(c[1,0])    
        X= param['L']/(2*math.pi)*np.mod(-q,2*math.pi)
    X=comm.bcast(X,0)
    logger.info('xa: {} xb: {} X: {}'.format(param['xa'],param['xb'],X))
    quit(1)
    return(X)

   

def calc_front(b):
    global xxf 
    global Nxx
    global param
    global Sx
   
    X=np.nan
    if param['Ny']>1: c=b.integrate('y','z')
    else:             c=b.integrate('z')
    if param['cm']:
        c=c['c']
        if mpirank==0: # This is correct for Fourier but not SinCos
            if param['Ny']>1: q=np.angle(c[1,0,0])
            else:             q=np.angle(c[1,0])    
            X= param['L']/(2*math.pi)*np.mod(-q,2*math.pi)
    else:
        c['c']=c['c']*Sx
        d=c['g']
        if mpirank==0:
            if d.ndim==3: e=d[:,0,0]
            else:         e=d[:,0]
            X=front1(e,xxf)
    X=comm.bcast(X,0)
    return(X)

def front1(e,x):
    btol=np.max(e)/20
    N=len(x)
    if e[-1]>btol: return(x[-1])
    for i in range(N-1,1,-1):
        b1=e[i-1]- btol
        b2=e[i]  - btol
        if e[i-1]>btol and e[i]<btol and e[i]>0: return((b1*x[i]-b2*x[i-1])/(b1-b2))
    return(np.nan)

# Look for an inflection point in the cubic and a strictly decreasing region
# This is only going to work well for a fully resolved simulations
def front2(e,x):
    N=len(x)
    for i in range(N-1,1,-1):
        e1=max(0,e[min(N-1,i+0)])
        e2=max(0,e[min(N-1,i+1)])
        e3=max(0,e[min(N-1,i+2)])
        e4=max(0,e[min(N-1,i+3)])
        x1=x[min(N-1,i+1)]
        x2=x[min(N-1,i+2)]
        dd1 = e1-2*e2+e3
        dd2 = e2-2*e3+e4
        if dd1*dd2<0: return((dd1*x2-dd2*x1)/(dd1-dd2))
    return(np.nan)

# Use pchip interpolation to find the front
def front3(f,x):
    ftol=np.max(e)/200
    e=it.pchip_interpolate(x, np.maximum(ftol,f), x, der=2)
    N=len(x)
    if e[-1]<0: return(x[-1])
    for i in range(N-1,1,-1):
        e1=e[i-1]
        e2=e[i] 
        if e1<0 and e2>=0: return((e1*x[i]-e2*x[i-1])/(e1-e2))
    return(np.nan)

def Cn(n):
    return(2*sp.poch(n+1,1/2)/sqrtpi)

def Dn(n):
    return(4**n/sp.binom(2*n,n)/(2*n+1)-1/(2*n+2))

def clip01(x):
    x=np.maximum(0,np.minimum(1,x))
    return(x)

# odd polynomials f(-1)=-1, f(1)=1, sith n vanishing derivative at 1 and -1
#for n in range(0,20): a,b,c=ju.zero_gradf(n,n)
# We clip the range for large values of n
# So the gradient is df=0 for |x|>rg
#                     f=-1 for x<-rg, f=1 for x>rg
#                     F=0 C-x for x<-rg F =C+x for x>rg
# rg =1 for small n and decreases for larger n
def zero_gradf(x,n):
    if n>100:
        n=100
        logger.info("zero_gradf: called with n={}, clipped to 100".format(n))
    
    C=Cn(n)
    D=Dn(n)
    #rg=np.minimum(1,1-np.power(np.finfo(np.float).eps/C,1/np.float(n+1)))
    #rg=np.minimum(1,np.sqrt(1-np.power(0.01,1/np.float(n))))
    s=np.sign(x)
    x=abs(x)
    z=np.minimum(1,x)
    z1=1-z
    z2=1+z
    y=z*z
    
    p=np.ndarray(n+1,dtype=np.float)
    q=np.ndarray(n+2,dtype=np.float)
    p2=np.ndarray(n+1,dtype=np.float)
    q2=np.ndarray(n+1,dtype=np.float)
    q[n+1]=-C*D
    for j in range(0,n+1):  p[j] = C*sp.binom(n, j)*(-1)**(n-j)/(1+2*(n-j))
    for j in range(0,n+1): p2[j] = -C*(-1)**n*sp.binom(n, j)*(-2)**j/(1+2*n-j)
    for j in range(0,n+1): q[j]  =  p[j]/(2*(n+1-j))
    for j in range(0,n+1): q2[j] =  (-1)**n*C*sp.binom(n, j)*(-2)**j/(1+2*n-j)/(2+2*n-j)
    
    df = C*(z1*z2)**n
    f  = s*np.maximum(0,np.minimum(1,z*np.polyval(p,y)))
    g  = s*np.maximum(0,np.minimum(1,1+z1**(n+1)*np.polyval(p2,z1)))
    F  = np.polyval(q,y)
    G  = -z1+z1**(n+2)*np.polyval(q2,z1)+np.maximum(x-1,0)
    rg1=0.5
    rg2=0.5

    #f=C*s*z**(2*n+1)*sp.hyp2f1(-n, -n-0.5,0.5-n,1/z**2)/(2*n+1)
    #F=  C*(z**(2*n+2)*sp.hyp2f1(-n-1, -n-0.5, .5-n, 1/z**2)/((2*n+1)*(2*n+2))-D-1/(2*n+2))

    #if len(x)>1:
    #try:
    #    rg1=min(x[np.where(np.diff(f,1)<0)])
    #except:
    #    rg=0.5
    #try:
    #    rg2=max(x[np.where(np.diff(g,1)<0)])
    #except:
    #    rg2=0.5
    #rg=(rg1+rg2)/2
    rg=0.5
    if hasattr(x, "__len__"):
        f[np.where(x>=rg)]=g[np.where(x>=rg)]
        F[np.where(x>=rg)]=G[np.where(x>=rg)]
    elif x>rg:
        f=g
        F=G
        
    #logger.info('n {} rg {} rg1 {} rg2 {}'.format(n,rg,rg1,rg2))
    if n==0: df[np.where(x>=1)]=0

    
    return(f,F,df)

# calculate entropy
def xlogx(x):
    if np.isscalar(x):
        if x<=0: y = 0
        else:    y = x*np.log(x)
    else:
        y=np.zeros_like(x)
        y[np.where(x>0)]=x[np.where(x>0)]*np.log(x[np.where(x>0)])
    return(y)
    
def mix_entropy(x):
    e=np.maximum(0,-xlogx(x)-xlogx(1-x))
    return(e)

def zero_gradf2(x,n):
    C=Cn(n)
    rg=np.minimum(1,np.maximum(0,1-np.power(np.finfo(np.float).eps/C,1/np.float(n))))
    rg=1
    z=np.minimum(1,np.maximum(-1,x))
    z1=1-z
    z2=1+z
    y=z*z
    
    p=np.ndarray(n+1,dtype=np.float)
    p2=p
    C2=0
    for j in range(0,n+1): C2=C2+4**n*(-1)**n/(n+j+1)*sp.binom(n,j)
    for j in range(0,n+1):  p[j]= C*sp.binom(n, j)*(-1)**(n-j)/(1+2*(n-j))
    for j in range(0,n+1): p2[j]=  -sp.binom(n, j)*(-2)**j/(1+2*n-j)/C2
    q=np.ndarray(n+2,dtype=np.float)
    q[n+1]=-C*Dn(n)
    for j in range(0,n+1):q[j]=p[j]/(2*(n+1-j))
    df = C*(z1*z2)**n
    f  = z*np.polyval(p,y)
    f2 = 1-z1**np.polyval(p2,z1)
                  
    #rg = np.sqrt(np.min(y[np.where(abs(f) >= 1)]))
    #f[np.where(x>=rg)]=np.float(1)
    #f[np.where(x<=-rg)]=np.float(-1)
    
    y=np.minimum(y,rg**2)
    F  = np.polyval(q,y)+np.maximum(x-rg,0)-np.minimum(rg+x,0)

    #logger.info("n {} rg {} 1-max(abs(f)) {} ".format(n,rg,1-np.max(abs(f))))
    
    return(f,F,df,f2)


# These functions are particularly good with Chebychev since they are exactly representable 
# With this scaling they get sharper and sharper with n
def smoothi(y,n,x1,x2): # interpolate between between 0 and 1 at x1 and x2 wih n vanishing derivatives at end points
    if x1==x2:
        f  = (1-np.sign(y-x1))/2
        df = 0*f
        F  = np.maximum(0,y-x1)
    else:
        S=2/(x2-x1)
        x=np.minimum(1,np.maximum(-1,(2*y-x1-x2)/(x2-x1)))
        [f,F,df]=zero_gradf(x,n)
        f=(1+f)/2
        F=(1+x+F)/(2*S)+np.maximum(y-max(x1,x2),0)
        df=df*S/2
    return(f,F,df)

def stepw(y,n,x1,x2,x3,x4,L): 
    f1a,F1a,df1a=smoothw(y,n,x1,x2)
    f2a,F2a,df2a=smoothw(y,n,x3,x4)
    f1b,F1b,df1b=smoothw(y,n,x1+L,x2+L)
    f2b,F2b,df2b=smoothw(y,n,x3+L,x4+L)
    f1c,F1c,df1c=smoothw(y,n,x1-L,x2-L)
    f2c,F2c,df2c=smoothw(y,n,x3-L,x4-L)
    return(f1a+f1b+f1c-f2a-f2b-f2c,F1a+F1b+F1c-F2a-F2b-F2c,df1a+df1b+df1c-df2a-df2b-df2c)

def stepwf(y,n,x1,x2,x3,x4,L): 
    f1a=smoothwf(y,n,x1,x2)
    f2a=smoothwf(y,n,x3,x4)
    f1b=smoothwf(y,n,x1+L,x2+L)
    f2b=smoothwf(y,n,x3+L,x4+L)
    f1c=smoothwf(y,n,x1-L,x2-L)
    f2c=smoothwf(y,n,x3-L,x4-L)
    return(f1a+f1b+f1c-f2a-f2b-f2c)

def step6(y,n,x1,x2,x3,x4,x5,x6,U0,U1,U2): 
    if   np.isfinite(x1) and np.isfinite(x6):
        f1,F1,df1=smoothw(y,n,x1,x2)
    else: 
        f1=1
        F1=y
    df1=0
    f2,F2,df2=smoothw(y,n,x3,x4)
    if   np.isfinite(x5) and np.isfinite(x6):
        f3,F3,df3=smoothw(y,n,x5,x6)
    else:
        f3=0
        F3=0
        df3=0
        
    f  = U0*( f3- f1)+U1*( f1- f2)+U2*( f3- f2)
    df = U0*(df3-df1)+U1*(df1-df2)+U2*(df3-df2)
    F  = U0*( F3- F1)+U1*( F1- F2)+U2*( F3- F2)
    
    return(f,F,df)

# Change this to use incomplete beta functions
def smoothw(y,n,x1,x2): # interpolate between between approximately 0 and 1 at x1 and x2 wih n vanishing derivatives at end points
                        # But with a fixed width, that is df(0)=1
                        # so that actual point of vanishing is further apart but
                        # is much more spectrally compact
    if x1==x2:
        f  = (1-np.sign(y-x1))/2
        df = 0*f
        F  = np.maximum(0,y-x1)
    else:
        [cf,cF,cdf]=zero_gradf(0,n)
        dx=(x2-x1)*cdf/4
        mx=(x1+x2)/2
        x=np.minimum(1,np.maximum(-1,(y-mx)/dx))
        [f,F,df]=zero_gradf(x,n)
        f=(1+f)/2
        F=(1+x+F)*dx/2+np.maximum(y-mx-np.abs(dx),0)
        df=df/(2*dx)
    return(f,F,df)

#Same as above but just return the function
def smoothwf(y,n,x1,x2): # interpolate between between approximately 0 and 1 at x1 and x2 wih n vanishing derivatives at end points
    f,F,df=smoothw(y,n,x1,x2)
    return(f)

def PGF_w(s):
    cs=np.cos(s/2)**6
    ss=np.sin(s/2)**2
    y=cs*(1+3*ss+6*ss**2);
    return(y)

def PGF_WW(s):
    ss=np.sin(s)
    y=s/2+ss*(120+20*ss**2+9*ss**4)/240;
    return(y)

def PGF(s,a):
    z=intexpcosabx(s,a,0)/2+75/128*intexpcosabx(s,a,1)-25/256*intexpcosabx(s,a,3)+3/256*intexpcosabx(s,a,5)
    return(z)
    
def periodic_gaussian(d,y,z,h,W,H):
    
    pi=math.pi
    s2=np.sqrt(2)*pi

    Y=pi*2*y/W 
    Z=pi*2*z/H 

    d0=1/h
    Q0=h**2*pi
    
    E   =  d**2*np.exp(-(H**2+W**2)* d**2/2)
    E0  = d0**2*np.exp(-(H**2+W**2)*d0**2/2)
    EY  = (PGF_WW(Y)*PGF_w(Z)-Y)*E
    EZ  = (PGF_WW(Z)*PGF_w(Y)-Z)*E
    HW0  = 3/8*pi**3*(H+W)*E0
    HW   = 3/8*pi**3*(H+W)*E
    dH0  = d0*H/s2
    dW0  = d0*W/s2
    dH   =  d*H/s2
    dW   =  d*W/s2
    
    GZ   = W/pi*d**2*PGF(Y,dW)*PGF_w(Z)*np.exp(-2*z**2*d**2)
    GY   = H/pi*d**2*PGF(Z,dH)*PGF_w(Y)*np.exp(-2*y**2*d**2)
    
    FF0 = H*W*d0**2*PGF(pi,dH0)*PGF(pi,dW0)
    FF  = H*W* d**2*PGF(pi,dH )*PGF(pi,dW )
    
    C= Q0 * pi**2/4                /(HW0+FF0)
    c= Q0/(H*W) * (HW0-HW + FF0-FF)/(HW0+FF0)
    
    psiy = y*c/2 +C*(GZ-EY)
    psiz = z*c/2 +C*(GY-EZ)
    return(psiy,psiz,c)

    #logger.info('d:{}, y:{}, z:{}, h:{}, W:{}, H:{}'.format(d,y,z,h,W,H))
    #logger.info('    E    = {:e}'.format( E    ))
    #logger.info('    E0   = {:e}'.format( E0   ))
    #logger.info('    EY   = {:e}'.format( EY   ))
    #logger.info('    EZ   = {:e}'.format( EZ   ))
    #logger.info('    HW0  = {:e}'.format( HW0  ))
    #logger.info('    HW   = {:e}'.format( HW   ))
    #logger.info('    dH0  = {:e}'.format( dH0  ))
    #logger.info('    dW0  = {:e}'.format( dW0  ))
    #logger.info('    dH   = {:e}'.format( dH   ))
    #logger.info('    dW   = {:e}'.format( dW   ))
    #logger.info('    GZ   = {:e}'.format( GZ   ))
    #logger.info('    GY   = {:e}'.format( GY   ))
    #logger.info('    FF0  = {:e}'.format( FF0  ))
    #logger.info('    FF   = {:e}'.format( FF   ))
    #logger.info('    C    = {:e}'.format( C    ))
    #logger.info('    c    = {:e}'.format( c    ))
    #logger.info('    psiy = {:e}'.format( psiy ))
    #logger.info('    psiz = {:e}'.format( psiz ))
  


def periodic_gaussian_uvw(d,y,z,h,W,H):
    
    pi=math.pi
    s2=np.sqrt(2)*pi

    Y=pi*2*y/W 
    Z=pi*2*z/H 

    d0=1/h
    Q0=h**2*pi

    rs=(y**2+z**2)*d**2/2
    Rs=(H**2+W**2)*d**2/2
    R0=(H**2+W**2)*d0**2/2
    E   =  d**2*np.exp(-Rs)
    E0  = d0**2*np.exp(-R0)
    Er  =  d**2*np.exp(-rs)
    EY  = (PGF_WW(Y)*PGF_w(Z)-Y)*E
    EZ  = (PGF_WW(Z)*PGF_w(Y)-Z)*E
    HW0  = 3/8*pi**3*(H+W)*E0
    HW   = 3/8*pi**3*(H+W)*E
    dH0  = d0*H/s2
    dW0  = d0*W/s2
    dH   =  d*H/s2
    dW   =  d*W/s2

    wZ   =  PGF_w(Z)
    wY   =  PGF_w(Y)
    
    GZ   = W/pi*d**2*PGF(Y,dW)*PGF_w(Z)*np.exp(-2*z**2*d**2)
    GY   = H/pi*d**2*PGF(Z,dH)*PGF_w(Y)*np.exp(-2*y**2*d**2)
    
    FF0 = H*W*d0**2*PGF(pi,dH0)*PGF(pi,dW0)
    FF  = H*W* d**2*PGF(pi,dH )*PGF(pi,dW )
    
    C= Q0 * pi**2/4                /(HW0+FF0)
    c= Q0/(H*W) * (HW0-HW + FF0-FF)/(HW0+FF0)

    u=c+4*C*Er*wZ*wY+2*pi*C*(1-wZ*wY)*(1/H+1/W)

    WdF=np.sqrt(2)*d*W*dFda(pi,d*W/s2)
    HdF=np.sqrt(2)*d*H*dFda(pi,d*H/s2)
    WF = PGF(pi,d*W/s2)
    HF = PGF(pi,d*H/s2)
 #   dcdx=3*h^2*((        -2+(H^2+W^2)*pi^4*(H+W)*E        -HF/3*(WdF+pi*WF)       +WF/3*(HdF*WF*HdF*W*H)*(D(d))(x)*d(x)/(8*W*d(0)^2*(3*pi^3*(H+W)*exp(-(1/2)*(H^2+W^2)*d(0)^2)*(1/8)+H*W*F(pi, d(0)*H*sqrt(2)/(2*pi))*F(pi, d(0)*W*sqrt(2)/(2*pi)))*H)
    

    
    psiy = y*c/2 +C*(GZ-EY)
    psiz = z*c/2 +C*(GY-EZ)
    return(psiy,psiz,c)

# d/da intexpcosabx(x,a,b) =  -2a*int(x^2*exp(-a^2*x^2)*cos(b*x),x);
def intexpcosabxda(x,a,b):
    return(((b**2-2*a**2)*iabx(s,a,b)+exp(-(a*s)**2)*(2*s*cos(b*s)*a**2-b*np.sin(b*s)))/(2*a**3))

# int(exp(-a^2*x^2)*cos(b*x),x);
# This only works for b O(1)
def intexpcosabx(x,a,b):
    atol=0.1
    return(np.where(a>atol*b,intexpcosabx1(x,np.maximum(a,atol*b),b),intexpcosabx2(x,np.minimum(a,atol*b),b)))

def intexpcosabx1(x,a,b): # small values of b/a
    a=np.maximum(1e-300,a)
    c=b/(2*a)
    z=sqrtpi/(2*a)*np.exp(-c**2)*np.real(sp.erf(a*x-1j*c))
    return(z)

def intexpcosabx2(x,a,b): # large values of b/a 
    c=b/np.maximum(1e-300,2*a)
    z=np.where(a>0,sqrtpi/(2*np.maximum(1e-300,a))*np.exp(-(x*a)**2)*np.real(np.exp(-1j*x*b)*sp.wofz(c-1j*a*x)),np.where(b>0,np.sin(b*x)/np.maximum(1e-300,b),x))
    return(z)

def intexpcosabx3(x,a,b):
    y=(b*x)
    c=(a/b)**2
    z=np.zeros_like(x)
    C=1
    for n in range(10):
        if n>0: C=-C*c/n
        z=z+C*intxncos(y,2*n)
    return(z)

# series expansion of cos(b*x) before integration good for large y
def intexpcosabx4(x,a,b):
    erfy=sp.erf(a*x)
    e=(b/a)**2
    y=(a*x)**2
    z=(e**6-24*e**5+480*e**4-7680*e**3+92160*e**2-737280*e+2949120)*sqrtpi*erfy/(5898240*a)-(1/958003200)*e*x*(-239500800+(y**5+(11/2)*y**4+(99/4)*y**3+(693/8)*y**2+(3465/16)*y+10395/32)*e**5+(-132*y**4-594*y**3-2079*y**2-(10395/2)*y-31185/4)*e**4+(11880*y**3+41580*y**2+103950*y+155925)*e**3+(-665280*y**2-1663200*y-2494800)*e**2+(19958400*y+29937600)*e)*np.exp(-y)
    return(z)

#integral of x^n cos x
#fn := n -> subs({cos(x) = cx, sin(x) = sx}, int(f(x, n), x))
#for i from 0 to 60 do
#CoefficientList(simplify(1/x^(((-1)^i+1)/4)*subs(x=sqrt(x),diff(fn(i),cx)),symbolic),x,termorder=reverse);
#end do
#for i from 0 to 60 do
#CoefficientList(simplify(1/x^((1-(-1)^i)/4)*subs(x=sqrt(x),diff(fn(i),sx)),symbolic),x,termorder=reverse);
#end do
def intxncos(x,n):
    if   n==0:  pc = []
    elif n==1:  pc = [1]
    elif n==2:  pc = [2]
    elif n==3:  pc = [3,-6]
    elif n==4:  pc = [4,-24]
    elif n==5:  pc = [5,-60,120]
    elif n==6:  pc = [6,-120,720]
    elif n==7:  pc = [7,-210,2520,-5040]
    elif n==8:  pc = [8,-336,6720,-40320]
    elif n==9:  pc = [9,-504,15120,-181440,362880]
    elif n==10: pc = [10,-720,30240,-604800,3628800]
    elif n==11: pc = [11,-990,55440,-1663200,19958400,-39916800]
    elif n==12: pc = [12,-1320,95040,-3991680,79833600,-479001600]
    elif n==13: pc = [13,-1716,154440,-8648640,259459200,-3113510400,6227020800]
    elif n==14: pc = [14,-2184,240240,-17297280,726485760,-14529715200,87178291200]
    elif n==15: pc = [15,-2730,360360,-32432400,1816214400,-54486432000,653837184000,-1307674368000]
    elif n==16: pc = [16,-3360,524160,-57657600,4151347200,-174356582400,3487131648000,-20922789888000]
    elif n==17: pc = [17,-4080,742560,-98017920,8821612800,-494010316800,14820309504000,-177843714048000,355687428096000]
    elif n==18: pc = [18,-4896,1028160,-160392960,17643225600,-1270312243200,53353114214400,-1067062284288000,6402373705728000]
    elif n==19: pc = [19,-5814,1395360,-253955520,33522128640,-3016991577600,168951528345600,-5068545850368000,60822550204416000,-121645100408832000]
    elif n==20: pc = [20,-6840,1860480,-390700800,60949324800,-6704425728000,482718652416000,-20274183401472000,405483668029440000,-2432902008176640000]
    elif n==21: pc = [21,-7980,2441880,-586051200,106661318400,-14079294028800,1267136462592000,-70959641905152000,2128789257154560000,-25545471085854720000,51090942171709440000]
    else:
        logger.info("intxncos: only defined for n upto 20 {}".format(n))
        z=[]

        
    if   n==0:  sc = [1]
    elif n==1:  sc = [1]
    elif n==2:  sc = [1,-2]
    elif n==3:  sc = [1,-6]
    elif n==4:  sc = [1,-12,24]
    elif n==5:  sc = [1,-20,120]
    elif n==6:  sc = [1,-30,360,-720]
    elif n==7:  sc = [1,-42,840,-5040]
    elif n==8:  sc = [1,-56,1680,-20160,40320]
    elif n==9:  sc = [1,-72,3024,-60480,362880]
    elif n==10: sc = [1,-90,5040,-151200,1814400,-3628800]
    elif n==11: sc = [1,-110,7920,-332640,6652800,-39916800]
    elif n==12: sc = [1,-132,11880,-665280,19958400,-239500800,479001600]
    elif n==13: sc = [1,-156,17160,-1235520,51891840,-1037836800,6227020800]
    elif n==14: sc = [1,-182,24024,-2162160,121080960,-3632428800,43589145600,-87178291200]
    elif n==15: sc = [1,-210,32760,-3603600,259459200,-10897286400,217945728000,-1307674368000]
    elif n==16: sc = [1,-240,43680,-5765760,518918400,-29059430400,871782912000,-10461394944000,20922789888000]
    elif n==17: sc = [1,-272,57120,-8910720,980179200,-70572902400,2964061900800,-59281238016000,355687428096000]
    elif n==18: sc = [1,-306,73440,-13366080,1764322560,-158789030400,8892185702400,-266765571072000,3201186852864000,-6402373705728000]
    elif n==19: sc = [1,-342,93024,-19535040,3047466240,-335221286400,24135932620800,-1013709170073600,20274183401472000,-121645100408832000]
    elif n==20: sc = [1,-380,116280,-27907200,5079110400,-670442572800,60339831552000,-3379030566912000,101370917007360000,-1216451004088320000,2432902008176640000]
    elif n==21: sc = [1,-420,143640,-39070080,8204716800,-1279935820800,140792940288000,-10137091700736000,425757851430912000,-8515157028618240000,51090942171709440000]

    y=x*x
    if np.mod(n,2)==0: z = x*np.cos(x)*np.polyval(pc,y) +   np.sin(x)*np.polyval(sc,y)
    if np.mod(n,2)==1: z =   np.cos(x)*np.polyval(pc,y) + x*np.sin(x)*np.polyval(sc,y)
    return(z)


# Differentiation formula accurate to 8th order
# Dont know how to handle anonymous functions and variable arguments so one function
# is defined for each number of arguments.
# Ugly but simple
def ndif0(f,x,dx):
    f4 = f(x+4*dx)-f(x-4*dx)
    f3 = f(x+3*dx)-f(x-3*dx)
    f2 = f(x+2*dx)-f(x-2*dx)
def theta_F2(r,a,L):
    s=a*L
    z=a*r
    f=r/np.maximum(1e-20,a)*theta_F1(r,a,L)
    for j in range(1,theta_N+1): f=f-2*L*j/sqrtpi*np.exp(-(z**2+(j*s)**2))*np.sinh(np.maximum(-700,np.minimum(700,2*s*z*j)))
    return(f)
    
def theta_periodic_gaussian(a,da,y,z,h,U,W,H):
    Fz=theta_F1(z,a,H)
    Fy=theta_F1(y,a,W)
    u = U*Fz*Fy
    v =-U*da*Fz*theta_F2(y,a,W)
    w =-U*da*Fy*theta_F2(z,a,H)
    return(u,v,w)

def flush_stats():
    global fpstats
    if mpirank==0: fpstats.flush()
    
def open_stats(fn,p,reset,problem):
    global basedir
    global abort
    global nstats
    global fpstats
    global fnsnm
    dtstats = dp.get_param(p,'dtstats')
    if dtstats is None: return

    fpstats=None
    nstats=None
    fnsnm=('n','it','t','T','dt','RSSMin','VMSMin','RSSMax','VMSMax','RSSMean','VMSMean',)
    var=problem.variables
    par=problem.parameters
    
    PIDX   = dp.get_param(p,'PIDX')
    PIDG   = dp.get_bool_param(p,'PIDG')
    dd     = dp.get_bool_param(p,'dd')
    db     = dp.get_bool_param(p,'db')
    ds     = dp.get_bool_param(p,'ds')
    ddserf = dp.get_bool_param(p,'ddserf')
    dbserf = dp.get_bool_param(p,'dbserf')
    dsserf = dp.get_bool_param(p,'dsserf')
    if 'b'  in var: fnsnm += ('minb','maxb',)#'Tb',
    if 's'  in var: fnsnm += ('mins','maxs',)#'Ts',
    if 'u'  in var: fnsnm += ('minu','maxu',)#'Tu',
    if 'v'  in var: fnsnm += ('minv','maxv',)#'Tv',
    if 'w'  in var: fnsnm += ('minw','maxw',)#'Tw',
    if 'EFB' in problem.parameters: fnsnm += ('FB',)
    if 'EFS' in problem.parameters: fnsnm += ('FS',)
    if 'EFU' in problem.parameters: fnsnm += ('FU',)
    if 'EFV' in problem.parameters: fnsnm += ('FV',)
    if 'EFW' in problem.parameters: fnsnm += ('FW',)
    if 'EWN' in problem.parameters: fnsnm += ('NS',)
    if 'Fx'  in problem.parameters: fnsnm += ('Fx',)
    if 'Fy'  in problem.parameters: fnsnm += ('Fy',)
    if 'Fz'  in problem.parameters: fnsnm += ('Fz',)
    if p['sType']=='gc': fnsnm += ('X',)
    if PIDX is not None:
        if PIDG: fnsnm += ('g',)
        else:    fnsnm += ('U',)
    if dd and db : fnsnm += ('dbu',)
    if dd and ds : fnsnm += ('dsu',)
    if db: fnsnm += ('dbmin','dbmax','dbh','dbw','db0','db1')
    if ds: fnsnm += ('dsmin','dsmax','dsh','dsw','ds0','ds1')
    if dd: fnsnm += ('ddmin','ddmax','ddh','ddw','dd0','dd1')
    if dp.get_bool_param(p,'divyz'): fnsnm += ('tops','topd')
    logger.info('open_stats: datasets {}'.format(fnsnm))
    if mpirank==0: 
        try:
            fn = basedir / 'stats.hdf5'
            if fn.is_file() and not reset:
                fpstats = h5py.File(str(fn), mode='a')
                nstats=len(fpstats['n'])
                for n in fnsnm:
                    if not n in fpstats: fpstats.create_dataset(n, maxshape=(None,), shape=(nstats,), dtype=np.float,  chunks=True)
            else:
                fpstats = h5py.File(str(fn), mode='w')
                nstats=0
                for n in fnsnm: fpstats.create_dataset(n, maxshape=(None,), shape=(0,), dtype=np.float, chunks=True)
                logger.info('Opened {}, n:{}'.format(fn,nstats))
            for m in fpstats:
                if not m in fnsnm:
                    logger.info('open_stats: deleting {} in {}'.format(m,fn))
                    del(fpstats[m])
        except:
            set_error('Failed to open stats file "{:}"'.format(fn))
    check_error()


def write_stats(problem,solver,domain,dt,X,U,g,se,gIxx,gIyy,gIzz,tops,topd):
    global fpstats
    global nstats
    global fnsnm
    global mem
    
    s=''

    MemTrack()
    if mpirank==0:
        if fpstats is not None:
            for n in fnsnm: fpstats[n].resize(nstats+1,axis=0)  
            fpstats['n'][nstats]    = nstats+1 
            fpstats['it'][nstats]   = solver.iteration
            fpstats['t'][nstats]    = solver.sim_time
            fpstats['T'][nstats]    = np.float(time.time())
            fpstats['dt'][nstats]   = dt
            fpstats['RSSMin'][nstats]  = mem['RSSMin']  
            fpstats['VMSMin'][nstats]  = mem['VMSMin']  
            fpstats['RSSMax'][nstats]  = mem['RSSMax']  
            fpstats['VMSMax'][nstats]  = mem['VMSMax']  
            fpstats['RSSMean'][nstats] = mem['RSSMean'] 
            fpstats['VMSMean'][nstats] = mem['VMSMean']  
            if 'bu'    in fnsnm: fpstats['dbu'][nstats]      = (gIzz*se['db']*se['dd']).sum()
            if 'su'    in fnsnm: fpstats['dsu'][nstats]      = (gIzz*se['ds']*se['dd']).sum()
            if 'dbh'    in fnsnm:
                B=problem.parameters['B']
                fpstats['dbh'][nstats]     = (gIzz*se['db']).sum()/B
                fpstats['dbw'][nstats]     = (gIzz*se['db']*(B-se['db'])).sum()/B**2
                fpstats['db0'][nstats]     = se['db'][0]
                fpstats['db1'][nstats]     = se['db'][-1]
                fpstats['dbmin'][nstats]   = se['db'].min()
                fpstats['dbmax'][nstats]   = se['db'].max()
            if 'dsh'    in fnsnm:
                S=problem.parameters['S']
                fpstats['dsh'][nstats]     = (gIzz*se['ds']).sum()/S
                fpstats['dsw'][nstats]     = (gIzz*se['ds']*(S-se['ds'])).sum()/S**2
                fpstats['ds0'][nstats]     = se['ds'][0]
                fpstats['ds1'][nstats]     = se['ds'][-1]
                fpstats['dsmin'][nstats]   = se['ds'].min()
                fpstats['dsmax'][nstats]   = se['ds'].max()
            if 'ddh'    in fnsnm:
                mindd=se['dd'].min()
                maxdd=se['dd'].max()
                sdd=(se['dd']-mindd)/(maxdd-mindd)
                fpstats['ddh'][nstats]     = (gIzz*sdd).sum()
                fpstats['ddw'][nstats]     = (gIzz*sdd*(1-sdd)).sum()
                fpstats['dd0'][nstats]     = se['dd'][0]
                fpstats['dd1'][nstats]     = se['dd'][-1]
                fpstats['ddmin'][nstats]   = mindd
                fpstats['ddmax'][nstats]   = maxdd
            if 'tops'  in fnsnm: fpstats['tops'][nstats]    = tops
            if 'topd'  in fnsnm: fpstats['topd'][nstats]    = topd
            if 'X'     in fnsnm: fpstats['X'][nstats]       = X
            if 'U'     in fnsnm: fpstats['U'][nstats]       = U 
            if 'g'     in fnsnm: fpstats['g'][nstats]       = g
            if 'Fx'    in fnsnm: fpstats['Fx'][nstats]      = problem.parameters['Fx']
            if 'Fy'    in fnsnm: fpstats['Fy'][nstats]      = problem.parameters['Fy']
            if 'Fz'    in fnsnm: fpstats['Fz'][nstats]      = problem.parameters['Fz']
            if 'X'     in fnsnm: s += ', X:{:6.3f}'.format(X)
        if 'U'    in fnsnm: s += ', U:{:6.4f}'.format(U) 
        if 'g'    in fnsnm: s += ', g:{:7.4f}'.format(g)

    for x in ('b','s','u','v','w'):
        if x in problem.variables:
            sx,minx,maxx= prange(x,solver.state[x])
            #Tx=jint(solver.state[x])
            if mpirank==0:
                s += sx
                if fpstats is not None:
                    fpstats['min'+x][nstats] = minx
                    fpstats['max'+x][nstats] = maxx
                    #fpstats['T'+x][nstats]   = Tx
    #if 'EWN' in problem.parameters:
    #x=wnoise*(noiseu*noiseu + noisev*noisev + noisew*noisew)
    #if noise>0: flow.add_property("wnoise*(noiseu*noiseu + noisev*noisev + noisew*noisew)",   name='NS')
    #y=np.sqrt(flow.volume_average('NS')/problem.parameters['EWN'])
    #   if mpirank==0:
    #       s += ', {}:{:6.4f}'.format('NS',y)
    #       if fpstats is not None: fpstats['NS'][nstats] = y
   
    #    for x in ('B','S','U','V','W'):
    #        wx='W'+x
    #        fx='F'+x
    #        if wq[wx] is not None:
    #            y=np.sqrt(flow.volume_average(fx)/wq[wx])
    #            if mpirank==0:
    #                s += ', {}:{:6.4f}'.format(fx,y)
    #                if fpstats is not None: fpstats[fx][nstats] = y
    fes,be,se,ue,ve,we=calc_uerror(problem,solver,domain)
    if mpirank==0:
        s=s+fes
        if fpstats is not None:
            if be is not None: fpstats['FB'][nstats] = be
            if se is not None: fpstats['FS'][nstats] = se
            if ue is not None: fpstats['FU'][nstats] = ue
            if ve is not None: fpstats['FV'][nstats] = ve
            if we is not None: fpstats['FW'][nstats] = we
        nstats += 1
    
    return(s)
    
def print_vars(a):
    for v in vars(a): logger.info('{}: {}'.format(v,eval('a.'+v)))
  
def close_stats():
    global fpstats
    if mpirank==0: fpstats.close()

def sizeof_fmt(num, suffix='B'):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

def print_item_sizes(a):
    for name, size in sorted(((name, sys.getsizeof(value)) for name,value in a),
                        key= lambda x: -x[1])[:10]:
        print("{:>30}: {:>8}".format(name,sizeof_fmt(size)))
    Barrier()
    
def log_all_sizes(a):
    global fnlog
    #locals().items()
    #globals().items()
    #print(a)
    #print(dir(pympler))
    for name, size in sorted(((name, asizeof.asizeof(value)) for name,value in a),
                            key= lambda x: -x[1])[:10]:
        append_file(fnlog ,"{:>30}: {:>8}\n".format(name,sizeof_fmt(size)))
        
# Set the forcing functions
# Make sure we are prepared for them to be updated as U changes
# it is important that fd,fu,fv and fw are copies and not references 
def set_ffun(problem,w,s,U,Tx,Ty,Tz,Uflg,AA,Px,JA,AAJ):
    f=None
    global Nx
    global Ny
    global Nz
    if w is not None:
        Py=1
        Pz=1

        if Nx>1 and 'x' not in w.meta: w.meta['x']={}
        if Ny>1 and 'y' not in w.meta: w.meta['y']={}
        if Nz>1 and 'z' not in w.meta: w.meta['z']={}
        
        if Nx>1:
            if 'parity' in w.meta['x']: Px = w.meta['x']['parity']
        if Ny>1:
            if 'parity' in w.meta['y']: Py = w.meta['y']['parity']
        if Nz>1:
            if 'parity' in w.meta['z']: Pz = w.meta['z']['parity']

        if Px is None: Px=1
        if Py is None: Py=1
        if Pz is None: Pz=1
        
        if Nx>1 and Tx=="SinCos": w.meta['x']['parity'] = Px
        if Ny>1 and Ty=="SinCos": w.meta['y']['parity'] = Py
        if Nz>1 and Tz=="SinCos": w.meta['z']['parity'] = Pz

        if Nx>1 and Tx=="Fourier": w.meta['x']['parity'] = 0
        if Ny>1 and Ty=="Fourier": w.meta['y']['parity'] = 0
        if Nz>1 and Tz=="Fourier": w.meta['z']['parity'] = 0

        if AA>1: aa=(1,AA)
        else:    aa=(1,)
        for a in aa:
            w.set_scales(a)     
            minx,maxx=jrange(w)
            logger.info('ffun ' + s+' {:3.1f} range: [{:9.6f}, {:9.6f}]'.format(a,minx,maxx))
            if a==AAJ: JA.savexyz('force/{:s}.hdf5'.format(s),s,w['g'],AAJ)
        if Uflg: f=copy.copy(w['g'])
        w['g']=U*w['g']
        problem.parameters[s] = w
        problem.parameters[s]['g']
    return(f)

# Set the weight functions
def set_wfun(problem,w,s,V,Tx,Ty,Tz,AA,JA,AAJ):
    if w is not None:
        if 'x' in w.meta:
            if   Tx=="SinCos":  w.meta['x']['parity'] = 1 
            elif Tx=="Fourier": w.meta['x']['parity'] = 0
        if 'y' in w.meta:
            if   Ty=="SinCos":  w.meta['y']['parity'] = 1 
            elif Ty=="Fourier": w.meta['y']['parity'] = 0
        if 'z' in w.meta:
            if   Tz=="SinCos":  w.meta['z']['parity'] = 1 
            elif Tz=="Fourier": w.meta['z']['parity'] = 0
        #        w['g']=np.maximum(0,np.minimum(1,w['g'])) 
        if AA>1: aa=(1,AA)
        else:    aa=(1,)
        for a in aa:
            w.set_scales(a)     
            minx,maxx=jrange(w)
            logger.info('wfun ' + s+' {:3.1f} range: [{:9.6f}, {:9.6f}]'.format(a,minx,maxx))
            if a==AAJ: JA.savexyz('force/{:s}.hdf5'.format(s),s,w['g'],AAJ)
 
        problem.parameters[s] = w
        problem.parameters[s]['g'] # This avoids some obscure bug
    
def mesh_info(domain,AA):
    debug=False
    r = comm.gather(list(domain.dist.grid_layout.start(scales=AA)),root=0)
    if mpirank != 0: r=np.ndarray((domain.dim,mpisize),dtype=np.int64)
    r = comm.bcast(r, root=0)
    r=np.asarray(r).transpose()
    if debug: logger.info('type(r) {}, r.shape {}'.format(type(r),r.shape))
    if debug: logger.info('r {} '.format(r))
    return(r)
 
def mesh_info2(domain,AA):
    debug=True
    r=np.ndarray((domain.dim,mpisize),dtype=np.int64)
    if mpirank == 0:
        r[:,0]=list(domain.dist.grid_layout.start(scales=AA))
        if debug: logger.info('{:4d}: {}'.format(0,r[:,0]))
    else:
        for j in range(1,mpisize):
            comm.send(domain.dist.grid_layout.start(scales=AA), dest=0)

    if debug: logger.info('mesh_info: send complete')
    if mpirank == 0:
        for j in range(1,mpisize):
            r[:,j]=comm.recv(source=j)
            if debug: logger.info('{:4d}: {}'.format(j,r[:,j]))

    r = comm.bcast(r, root=0)
    return(r)


def get_bases_type(b):
    tb=str(type(b))
    if tb=="<class 'dedalus.core.basis.Fourier'>":
        return('Fourier')
    elif tb=="<class 'dedalus.core.basis.Chebyshev'>":
        return('Chebyshev')
    elif tb=="<class 'dedalus.core.basis.SinCos'>":
        return('SinCos')

def print_bases_info(b):
    print('bases.type: {}'.format(get_bases_type(b)))
    for a in vars(b): print('bases.{}: {}'.format(a,eval('b.' + a)))   

def print_domain_info(domain,AA):
    if mpirank==0:
        dim=domain.dim
        #print("dir(domain):       {}".format(dir(domain)))
        #print("vars(distributor): {}".format(vars(domain.distributor)))
        #print("vars(dist):        {}".format(vars(domain.dist)))
        for a in vars(domain): print('domain.{}: {}'.format(a,eval('domain.' + a)))   
        for a in ('dim','dealias','hypervolume','global_coeff_shape','local_coeff_shape','dealias_buffer_size'):
            print('domain.{}: {}'.format(a,eval('domain.' + a)))   
        for b in domain.bases: print_bases_info(b)
        #for j in range(len(domain.bases)):
            #for a in vars(domain.bases[j]): print('domain.bases[{}].{}: {}'.format(j,a,eval('domain.bases[{}].{}'.format(j,a))))  
            #for a in vars(domain.bases[j].grid): print('domain.bases[{}].grid.{}: {}'.format(j,a,eval('domain.bases[{}].grid.{}'.format(j,a))))  
        print('domain.grid(1).shape:  {}'.format(domain.grid(1).shape))
        print('type(domain.grids(1)):  {}'.format(type(domain.grids(1))))
        for j in range(dim): print('domain.bases[{}].name: {}'.format(j,domain.bases[j]))
        print('domain.grids(1)[0][:,0,0] {} '.format(domain.grids(1)[0][:,0,0]))
        for a in domain.grids(1):  print('domain.grids({:3.1f}) type: {} shape: [{:3d} {:3d} {:3d}], rng [{:7.4f},{:7.4f}]'.format(1,type(a),a.shape[0],a.shape[1],a.shape[2],a.min(),a.max()))
        for a in domain.grids(AA): print('domain.grids({:3.1f}) type: {} shape: [{:3d} {:3d} {:3d}], rng [{:7.4f},{:7.4f}]'.format(AA,type(a),a.shape[0],a.shape[1],a.shape[2],a.min(),a.max()))
        print('domain.dist.rank:     {}'.format(domain.dist.rank))   
        print('domain.dist.size:     {}'.format(domain.dist.size))
        print('domain.dist.coords:   {}'.format(domain.dist.coords)) 
        print('domain.dist.mesh:     {}'.format(domain.dist.mesh))  
        for j in range(dim): print('domain.elements({}).shape:   {}'.format(j,domain.elements(j).shape))  
        print('domain.elements(0)[:,0,0]:   {}'.format(domain.elements(0)[:,0,0]))  
        print('domain.elements(1)[0,:,0]:   {}'.format(domain.elements(1)[0,:,0]))  
        print('domain.elements(2)[0,0,:]:   {}'.format(domain.elements(2)[0,0,:]))  
        for j in range(dim):
            for a in ('name','element_name','kx','ky','kz','base_grid_size','interval','dealias','grid','grid_dtype','coeff_dtype','elements','wavenumbers','coeff_size'):
                if a in vars(domain.bases[j]): print("domain.bases[{}].{}:   {}".format(j,a,eval('domain.bases[{}].{}'.format(j,a))))
        #print('domain.bases[{}]:   {}'.format(j,domain.bases[j]))
        #print('vars(domain.bases[{}]):   {}'.format(j,vars(domain.bases[j])))
        exit(1)
        #domain.bases: [<Fourier 140530315115152>, <Fourier 140530314755096>, <Chebyshev 140530314753864>]
        #domain.distributor: <dedalus.core.distributor.Distributor object at 0x7fcfc37db0b8>
        #domain.dist: <dedalus.core.distributor.Distributor object at 0x7fcfc37db0b8>
     
    print('                   ')
    print('domain.dist.grid_layout')
    print('  ext_coords     {}'.format(domain.dist.grid_layout.ext_coords))  
    print('  ext_mesh       {}'.format(domain.dist.grid_layout.ext_mesh))    
    print('  index          {}'.format(domain.dist.grid_layout.index))       
    print('  local          {}'.format(domain.dist.grid_layout.local))       
    print('  grid_space     {}'.format(domain.dist.grid_layout.grid_space))  
    print('                   ')
    print('  scales         {}'.format(1))
    print('    global_shape {}'.format(domain.dist.grid_layout.global_shape(scales=1))) 
    print('    local_shape  {}'.format(domain.dist.grid_layout.local_shape(scales=1))) 
    print('    slices       {}'.format(domain.dist.grid_layout.slices(scales=1)))         
    print('    start        {}'.format(domain.dist.grid_layout.start(scales=1)))           
    print('                   ')
    print('  scales         {}'.format(AA))
    print('    global_shape {}'.format(domain.dist.grid_layout.global_shape(scales=AA))) 
    print('    local_shape  {}'.format(domain.dist.grid_layout.local_shape(scales=AA))) 
    print('    slices       {}'.format(domain.dist.grid_layout.slices(scales=AA)))         
    print('    start        {}'.format(domain.dist.grid_layout.start(scales=AA)))           
    print('                   ')
    print('vars(domain.grid): {}'.format(vars(domain.grid)))
    print("domain.grid['axis']: {}".format(vars(domain.grid['axis'])))
    print('vars(domain.grids):{}'.format(vars(domain.grids)))
    print('domain.grids.axis:       {}'.format(domain.grids.axis))
    print('vars(domain):      {}'.format(vars(domain)))
    print('vars(distributor): {}'.format(vars(domain.distributor)))
    print('vars(dist):        {}'.format(vars(domain.dist)))
    
    print('vars(domain.dist.layout_references: {}'.format(dir(domain.dist.layout_reference)))
    print("vars(domain):      {}".format(vars(domain)))
    print("vars(distributor): {}".format(vars(domain.distributor)))
    print('vars(dist):        {}'.format(vars(domain.dist)))
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
    

def exit(a):
    #comm.MPI_Finalize()
    sys.exit(a)

# Doesn't work properly
def info_all(s):
    s=copy.copy(s)
    Barrier()
    if mpirank == 0: logger.info('{:4d}: {}'.format(mpirank,s))
    for j in range(1,mpisize):
        if mpirank == 0:
            r=comm.recv(source=j,tag=j)
            logger.info('{:4d}: {}'.format(j,r))
        else:  comm.send(s, dest=0,tag=mpirank)
    Barrier()


def remove_q_dep(y,z,f,n):
    Ny=y.size
    Nz=z.size
    y=np.reshape(y,(Ny,1))
    z=np.reshape(z,(1,Nz))
    q=np.arctan2(y,z)
    
    nn=2*n
    A=np.zeros((nn,  ),dtype=np.float)
    N=np.zeros((nn,nn),dtype=np.float)
    for j in range(n):
            A[j  ] = np.sum(np.sin((j+1)*q)*f)
            A[j+n] = np.sum(np.cos((j+1)*q)*f)
    for j in range(n):
        for k in range(n):
            N[j  ,k  ]=np.sum(np.sin((j+1)*q)*np.sin((j+1)*q))
            N[j+n,k  ]=np.sum(np.cos((j+1)*q)*np.sin((j+1)*q))
            N[j  ,k+n]=np.sum(np.sin((j+1)*q)*np.cos((j+1)*q))
            N[j+n,k  ]=np.sum(np.cos((j+1)*q)*np.cos((j+1)*q))
    
            
    A=comm.reduce(A, op=MPI.SUM,root=0)
    N=comm.reduce(N, op=MPI.SUM,root=0)

    if mpirank==0: c=np.reshape(np.linalg.solve(N,A),(nn,))
    else:          c=np.zeros((nn,),dtype=np.float)
    c = comm.bcast(c, root=0)
    
    for j in range(n):
            f=f-c[j]*np.sin((j+1)*q)-c[j+n]*np.cos((j+1)*q)
    return(f)


# project the function f onto the space spanned by the function w*x**j where j=0..n-1
def poly_proj(f,w,n):
    N=f.size
    x=np.linspace(-1,1, N, endpoint=True, dtype=np.float)
    A=np.zeros((n, ),dtype=np.float)
    N=np.zeros((n,n),dtype=np.float)
    for j in range(n): A[j] = np.sum(w*f*x**j)
 
    for j in range(n):
        for k in range(n):
            N[j,k]=np.sum(w**2*x**(j+k))
    c=np.linalg.solve(N,A)
    g=0*f
    for j in range(n): g += c[j]*w*x**j
    return(g)
    

    #mpiexec -n 4 ded_gc.py --preset ccle --Nx 256 --Ny 4 -W None -T 5 --Re 50 --dtjavrg 1 --dtjb 1 --dtju 1 --dtjv 1 --dtjw 1 --dtjy 0.1 --dtforce 1 --dtstats 1 --dtxyz 0 --dtmomb 0 --dtavrg 0 --dtgc 0 --dtleft 1 --dtright 1 --dtslice 0 --slicex 0 --slicey 0 --slicez 0 --dtb 0 --dtu 0 --dtv 0 --dtw 0 --dtx 0 --dty 0 --dtz 0 --dtxy 0 --dtyz 0 --dtxz 0 --lbc slip --ubc slip gc/test/05

# Calculate noise stepping parameters
# np.expm1(x)=exp(x)-1
# x' = a*x+b*W  where W is Gaussian noise std 1
# then x will approach gaussian with std sigma and autocorrelation time 1/lmbd
def noisedt(dt,lmbd,std):
    a=np.exp(-dt*lmbd);
    b=std*np.sqrt(1+a)*np.sqrt(-np.expm1(-dt*lmbd));
    return(a,b)

# Integrate a scalar or array for a noise equation dx/dt + lambda x = N
# With the noise chosen so that the standard deviation of x is std
# autocorrelation function decays as exp(-lambda*t)
def noiseint(x,dt,lmbd,std):
    (a,b)=noisedt(dt,lmbd,std)
    if np.isscalar(x):y=a*x+b*np.random.normal()
    else:             y=a*x+b*np.random.normal(size=x.shape())
    return(y)

def pinlet(y,z,inlet,R,FFTA):
    r  = np.sqrt(y*y+z*z)
    if inlet == "gaussian":
        w = 2*np.exp(-2*(r/R)**2)             # Make maximum and area pi*R^2 and int(f^2)=int(f)
        m = 2
    elif inlet == "circle":
        m = 1
        #        w  = smoothwf(r,10,R+W/2,R-W/2)
        #w  = ffta_heaviside(R-r,maxk,nfr)[0]
        w=FFTA.heaviside(R-r)[0]
        #write_hdf5('w.hdf','w',FFTA.heaviside(-r,np.sqrt(2))[0])
    else:
        set_error('pinlet: unknown inlet type {}'.format(inlet))
    w = m*clip01(w/m)

    return(w,m,r)

# Choose radius and scale factor so that
# int(f   dy dz) = A
# and
# int(f^2 dy dz) = A
def solve_inlet(y,z,SdA,inlet,R,FFTA):
    Ny=y.shape[1]
    Nz=z.shape[2]
    y=y.reshape([Ny,1])
    z=z.reshape([1,Nz])
    R0=copy.copy(R)
    dA=SdA/(math.pi*R**2)
    f = lambda RR: pinlet(y,z,inlet,RR,FFTA)[0].sum()**2*dA-(pinlet(y,z,inlet,RR,FFTA)[0]**2).sum()
    #logger.info('solve_inlet(f(R0[1/2 1 2]) {} {} {}'.format(f(R0/2),f(R0),f(2*R0)))
    R = opt.brentq(f,R0/2,2*R0)
    I1=( pinlet(y,z,inlet,R,FFTA)[0]    ).sum()*dA
    I2=((pinlet(y,z,inlet,R,FFTA)[0])**2).sum()*dA
    A = I1/I2
    return(R,A)
def parity_string(v,C=('x','y','z')):
        s='['
        for c in C:
            if c in v.meta:
                if 'parity' in v.meta[c]:
                    if v.meta[c]['parity'] is None: s = s + '   '
                    else:                           s = s + '{:3d}'.format(v.meta[c]['parity'])
                else:                               s = s + '  0'
            else:                                   s = s + '  0'
        return(s+']')

def print_parity(p,C=('x','y','z')):
    for v in p.variables:
        s='{:6s}'.format(v)
        for c in C:
            if c in p.meta[v]:
                if 'parity' in p.meta[v][c]:
                    if p.meta[v][c]['parity'] is None: s = s + ', ' + c + '   '
                    else:                              s = s + ', ' + c + '{:3d}'.format(p.meta[v][c]['parity'])
                else:                                  s = s + ', ' + c + '  0'
            else:                                      s = s + ', ' + c + '  0'
        logger.info('print_parity var: {}'.format(s))
    for v in p.parameters:
        x=p.parameters[v]
        if 'meta' in dir(x):
            s='{:6s}'.format(v)
            for c in ('x','y','z'):
                try:
                    s = s + ', ' + c + '{:3d}'.format(x.meta[c]['parity'])
                except:
                    s = s + ', ' + c + '  -'
            logger.info('print_parity par: {}'.format(s))
 

def grid_type(g):
    if   g.__doc__.find('Fourier')    !=-1: r='Fourier'
    elif g.__doc__.find('cos')        !=-1: r='SinCos'
    elif g.__doc__.find('Chebyshev')  !=-1: r='Chebyshev'
    else:
        logger.info('grid_type: {}'.format(g.__doc__))
        r='Unknown'
    return(r)

def nns(a,x):
    if x is not None: s=a.format(x)
    else: s=''
    return(s)

def nnsd(a,x,y):
    if x is not None and y is not None: s=a.format(x/y)
    else: s=''
    return(s)


def slice_odd(x):
    if x.start % 2 == 0 : s=slice(1+x.start,x.stop,2)
    else:                 s=slice(  x.start,x.stop,2) 
    return(s)

def slice_even(x):
    if x.start % 2 == 0 : s=slice(  x.start,x.stop,2)
    else:                 s=slice(1+x.start,x.stop,2) 
    return(s)

# Set the coefficients in the field xb to have a single hump at y of height 1
# I is the integral of the field
def ffta_single(f,domain,d,y):
    if 'parity' in f.meta[d]: parity=f.meta[d]['parity']
    else: parity=0
    dim     = domain.dim
    for k in range(dim):
        if domain.bases[k].name==d: break
    basis=domain.bases[k]
    if basis.name !=d:
        logger.info('Could not find the write direction')
    cslice  = domain.dist.coeff_layout.slices(scales=1)[k]
    lcshape = domain.dist.coeff_layout.local_shape(scales=1)[k]
    gcshape = domain.dist.coeff_layout.global_shape(scales=1)[k]
    N=gcshape-1
    w=basis.wavenumbers[cslice]
    i=np.arange(cslice.start,cslice.stop)
    sz=np.ones((dim,),dtype=np.int64)
    sz[k]=lcshape
    L=basis.interval[1]-basis.interval[0]
    #y=(y-basis.interval[0])/L
    #print(w/math.pi)
    A=  np.exp(sp.gammaln(2*N+1)-sp.gammaln(1+N+i)-sp.gammaln(1+N-i)-N*np.log(4))
    if parity==0: # periodic function
        f['c'] = np.reshape(np.exp(-1j*w*y)*A,sz)
        I=A[0]
    elif parity==1: # even function
        f['c'] = np.reshape(2*np.cos(w*y)*A/(1+np.heaviside(-i,1)),sz)
        I=A[0]
    elif parity==-1: # odd function
        f['c'] = np.reshape(4*np.sin(w*y)*A,sz)
        I=2*np.sum(f['c'][1:N+1:2]/w[1:N+1:2])
    return(f,I)

def ffta_single(f,domain,d,y):
    if 'parity' in f.meta[d]: parity=f.meta[d]['parity']
    else: parity=0
    dim     = domain.dim
    for j in range(dim):
        if domain.bases[j].name==d: k=j
        else: f.meta[domain.bases[j].name]['Constant']=True
    basis=domain.bases[k]
    if basis.name !=d:
        logger.info('Could not find the write direction')
    cslice  = domain.dist.coeff_layout.slices(scales=1)[k]
    lcshape = domain.dist.coeff_layout.local_shape(scales=1)[k]
    gcshape = domain.dist.coeff_layout.global_shape(scales=1)[k]
    N=gcshape-1
    w=basis.wavenumbers[cslice]
    i=np.arange(cslice.start,cslice.stop)
    sz=np.ones((dim,),dtype=np.int64)
    sz[k]=lcshape
    L=basis.interval[1]-basis.interval[0]
    A=  np.exp(sp.gammaln(2*N+1)-sp.gammaln(1+N+i)-sp.gammaln(1+N-i)-N*np.log(4))
    if parity==0: # periodic function
        f['c'] = np.reshape(np.exp(-1j*w*y)*A,sz)
        I=A[0]
    elif parity==1: # even function
        f['c'] = np.reshape(2*np.cos(w*y)*A/(1+np.heaviside(-i,1)),sz)
        I=A[0]
    elif parity==-1: # odd function
        f['c'] = np.reshape(4*np.sin(w*y)*A,sz)
        s=slice_odd(cslice)
        I=2*jsum(f['c'][s]/w[s])
    return(f,I)

def ffta_odd_constant(f,domain,d):

    if f.meta[d]['parity'] !=-1: set_error('ffta_odd_constant: Parity must be odd')
    dim     = domain.dim
    sz=np.ones((dim,),dtype=np.int64)
    for j in range(dim):
        if domain.bases[j].name==d: k=j
        else: f.meta[domain.bases[j].name]['Constant']=True
    #print('f[c].shape {}'.format(f['c'].shape ))   
    if domain.bases[k].name !=d: set_error('Could not find the write direction')
    check_error()
    s     = slice_odd(domain.dist.coeff_layout.slices(scales=1)[k])
    N     = domain.dist.coeff_layout.global_shape(scales=1)[k]
    i=np.arange(s.start,s.stop,2)
    f['c']=np.zeros((N,))
    sz[k] = len(i)
    
    f['c'][s]=np.reshape(4/(math.pi*i)*np.exp(2*sp.gammaln(N)-sp.gammaln(N+i)-sp.gammaln(N-i)),sz)
    return(f)


def ffta_single_NP(N,P,y,AA):

    if mpirank==0:
        sc=comm.Split(0,0)
        if P==0: g=de.Fourier('x', N, interval = (0,1), dealias=1)
        else:    g= de.SinCos('x', N, interval = (0,1), dealias=1)
        d=de.Domain([g], np.float64, comm=sc)
        f=d.new_field()
        f.meta['x']['parity']=P
        
        w=d.bases[0].wavenumbers
        N=len(w)
        
        i=np.arange(0,N)
        A=  np.exp(sp.gammaln(2*N+1)-sp.gammaln(1+N+i)-sp.gammaln(1+N-i)-N*np.log(4))
        if P==0: # periodic function
            f['c'] = np.exp(-1j*w*y)*A
            I=A[0]
        elif P==1: # even function
            f['c'] = 2*np.cos(w*y)*A/(1+np.heaviside(-i,1))
            I=A[0]
        elif P==-1: # odd function
            f['c'] = 4*np.sin(w*y)*A
            s=slice(1,N,2)
            I=2*jsum(f['c'][s]/w[s])
        f.set_scales(AA)
        z=np.array(copy.copy(f['g']))
        del(f)
    else:
        sc=comm.Split(1, mpirank)
        I=np.float64(0)
        z=np.array((np.int(AA*N),),dtype=np.float64)
    I   = comm.bcast(I,root=0)    
    z   = comm.bcast(z,root=0)    
    return(z,I)


#Function is 1 between y1 and y2 and zero outside this range
#Can be periodic, even or odd
def ffta_interval3(N,P,y1,y2,AA):

    if mpirank==0:
        sc=comm.Split(0,0)
        if P==0: g=de.Fourier('x', N, interval = (0,1), dealias=1)
        else:    g= de.SinCos('x', N, interval = (0,1), dealias=1)
        d=de.Domain([g], np.float64, comm=sc)
        f=d.new_field()
        if P !=0: f.meta['x']['parity']=P
        w=d.bases[0].wavenumbers
        N=len(w)
        w=w[slice(1,N)]
        n=np.arange(1,N)
        A=np.exp(2*sp.gammaln(2*N+1)-sp.gammaln(2*N+n+1)-sp.gammaln(2*N+1-n))/w
        if P==0: # periodic function But f does not endup periodic ???!!!???
            f['c'][1:N]=1j*A*(np.exp(-1j*w*y2)-np.exp(-1j*w*y1))
            f['c'][0]   = y2-y1
        elif P==1: # even function
            f['c'][1:N] = 2*A*(np.sin(w*y2)-np.sin(w*y1))
            f['c'][0]   = y2-y1
        elif P==-1: # odd function
            f['c'][1:N] = 4*A*(np.cos(w*y1)-np.cos(w*y2))
            s=slice(1,N,2)
        f.set_scales(AA)
        z=np.array(copy.copy(f['g']))
        del(f)
    else:
        sc=comm.Split(1, mpirank)
        z =np.array((np.int(AA*N),),dtype=np.float64)
    z = comm.bcast(z,root=0)    
    return(z)



# This returns a function that is 1 over as much of the domain as possible
# but is an odd function. i.e. it goes to zero at the end points
def ffta_odd_constant_N(N,A):
    K=np.int(N*A)
    if mpirank==0:
        sc=comm.Split(0,0)
        g=de.SinCos('x', N, interval = (0,1), dealias=1)
        d=de.Domain([g], np.float64, comm=sc)
        f=d.new_field()
        f.meta['x']['parity']=-1
        i=np.arange(1,N,2)
        s=slice(1,N,2)
        f['c']=np.zeros((N,))
        f['c'][s]=4/(math.pi*i)*np.exp(2*sp.gammaln(N)-sp.gammaln(N+i)-sp.gammaln(N-i))
        f.set_scales(A)
        z=np.array(copy.copy(f['g']))
        del(f)
    else:
        sc=comm.Split(1, mpirank)
        z=np.array((K,),dtype=np.float64)
    z   = comm.bcast(z,root=0)    
    return(z)

def ffta_odd_half_N(M,A):
    K=np.int(M*A)
    if mpirank==0:
        N=2*M
        sc=comm.Split(0,0)
        g=de.SinCos('x', N, interval = (0,1), dealias=A)
        d=de.Domain([g], np.float64, comm=sc)
        f=d.new_field()
        f.meta['x']['parity']=-1
        i=np.arange(1,N,2)
        s=slice(1,N,2)
        f['c']=np.zeros((N,))
        f['c'][s]=4/(math.pi*i)*np.exp(2*sp.gammaln(N)-sp.gammaln(N+i)-sp.gammaln(N-i))
        f.set_scales(A)
        z=np.array(copy.copy(f['g'][0:K]))
        del(f)
    else:
        sc=comm.Split(1, mpirank)
        z=np.array((K,),dtype=np.float64)
    z   = comm.bcast(z,root=0)    
    return(z)

# Single hump at zero width of hump falls off as 1/sqrt(N) but accuracy improves
def ffta_single_zero(x,P,N,M):
    I=np.float64(0)
    z=np.zeros_like(x)
    if mpirank==0:
        if P==1:
            s=np.sin(math.pi*x/2)**2
            c=np.cos(math.pi*x/2)**2
            a=1
            for n in range(M+1):
                z += a
                a *= (n+N-M)/(n+1)*s  # sp.pochhammer(N-M, n)/sp.factorial(n)
            z=z*c**(N-M)
        elif P==0:
            s=np.sin(math.pi*x)**2
            c=np.cos(math.pi*x)**2
            a=1
            N=np.int(N/2)
            for n in range(M+1):
                z += a
                a *= (n+N-M)/(n+1)*s  # sp.pochhammer(N-M, n)/sp.factorial(n)
            z=z*c**(N-M)
        elif P==-1: # Unclear what M would mean here
            s=np.sin(math.pi*x/2)
            c=np.cos(math.pi*x/2)
            z=c**(2*N-1)*s*np.sqrt(2*N-1)/(1-1/(2*N))**N
            
        I=z.sum()*(x[1]-x[0])
    I   = comm.bcast(I,root=0)    
    z   = comm.bcast(z,root=0)    
    return(z,I)
 

# Narrowest delta function so that z**2 is exactly represented on the grid 
# This is useful as a source term for divergence in u**2 terms are still exact
# Note that the width will decrease linearly with N for constant M
def ffta_delta2(x,y,P,N,M):
    I=np.float64(0)
    z=np.zeros((N,),dtype=np.float64)
    
    if mpirank==0:
        if   P==0: k = (N-2)*math.pi/(2*M)
        elif P==1: k = (N-1)*math.pi/(2*M)
        else:       set_error('ffta_delta: Parity must be 0 or 1')
        if   P== 1: x=np.minimum(math.pi/2,k*np.abs(np.mod(x-y+1,2)-1))
        elif P== 0: x=np.minimum(math.pi/2,k*np.abs(np.mod(x-y+1/2,1)-1/2))
        z=np.cos(x)**M
        I=z.sum()/N

    test_error()
    I   = comm.bcast(I,root=0)    
    z   = comm.bcast(z,root=0)    
    return(z,I)
   
    
# f(x) ~ delta
# This is primarily useful as a source term of flux
# It is not good as a forcing or limiting term since it doesn't exactly equal 1 at the grid points
# coefficient calculated by minimising int(f^2,x=1/2,1) over a fourier series with M+1 coefficients
# The function values are then calculated and recorded below
# need to take account that odd and even grids at at cell midpoints not boundaries
# Area of humps integrates to 1 
def ffta_delta3(P,N,M):
    if M>9:
        set_error('ffta_hump: M must be less than 10 P: {}, M: {}, N:{}'.format(P,M,N))
    if P!=-1 and P!=0 and P!=1:
        set_error('ffta_hump: P must be -1, 0, or 1 10 P: {}, M: {}, N:{}'.format(P,M,N)) 
    test_error()

    CE=([+4.098527206297015e-01],
        [+5.960444128326279e-01,-5.978786093400966e-02],
        [+7.356134969417889e-01,+0.000000000000000e+00,+7.233264033600564e-03],
        [+8.003752740192121e-01,+9.045458763111041e-02,+7.490901627989468e-04,-1.045187662449601e-03],
        [+8.387681994793682e-01,+1.744289039673268e-01,+0.000000000000000e+00,-1.593005091685743e-04,+1.617740801765207e-04],
        [+8.644712409215425e-01,+2.473337498551971e-01,+9.205229605441995e-03,+5.878488993339575e-05,+2.856800236299133e-05,-2.590555585444840e-05],
        [+8.829812517493421e-01,+3.098253541301349e-01,+2.594766363754438e-02,+0.000000000000000e+00,-8.070601584576941e-06,-4.959842786193255e-06,
         +4.229685213908158e-06],
        [+8.969841861405347e-01,+3.634673900281573e-01,+4.756211604497414e-02,+7.622230225249363e-04,+1.455397265069412e-06,+1.102470741243741e-06,
         +8.530149216403471e-07,-6.990505460705355e-07],
        [+9.079636078788050e-01,+4.097850761501867e-01,+7.202336813606877e-02,+3.011994184305247e-03,+0.000000000000000e+00,-4.256617274682213e-08,
         -1.558521760858841e-07,-1.462332109581450e-07,+1.164793046059847e-07],
        [+9.168110450614118e-01,+4.500677592597160e-01,+9.792682378416061e-02,+7.001609893896628e-03,+5.538451271435793e-05,+8.690073234837558e-09,
         -9.668485456810650e-09,+2.271786484524628e-08,+2.504073057437416e-08,-1.951996835609740e-08],
        [+9.240965191388775e-01,+4.853604737561212e-01,+1.243335824547687e-01,+1.272239030994610e-02,+2.976242921697164e-04,+0.000000000000000e+00,
         +5.869560862863153e-09,+3.095086555623624e-09,-3.393124387057006e-09,-4.286524906692770e-09,+3.284955271866932e-09],
        [+9.302023253549351e-01,+5.165002982049944e-01,+1.506301751576533e-01,+2.002701458703150e-02,+8.625091347585929e-04,+3.670762850645800e-06,
         -8.647486840781930e-10,-1.180847509240107e-09,-6.611800599864965e-10,+5.166200114751062e-10,+7.337695212956215e-10,-5.545634878745418e-10],
        [+9.353948162031303e-01,+5.441575967423303e-01,+1.764248316337823e-01,+2.870716695574022e-02,+1.862600999954534e-03,+2.620245413078129e-05,
         +0.000000000000000e+00,+2.279938082650908e-10,+1.850998747556491e-10,+1.258777207964574e-10,-7.987184534147236e-11,-1.256218051858424e-10,
         +9.384966083521305e-11],
        [+9.398654444541540e-01,+5.688717499337616e-01,+2.014753200523025e-01,+3.853826487188075e-02,+3.377418348283654e-03,+9.331325621978128e-05,
         +2.270429288786120e-07,-4.669253051399142e-11,-2.060138723200930e-11,-2.717757551706306e-11,-2.282896113106529e-11,+1.250265324424802e-11,
         +2.151002336807699e-11,-1.591292798006716e-11],
        [+9.437554584544370e-01,+5.910797075142037e-01,+2.256398095301727e-01,+4.930333255133965e-02,+5.453795220573187e-03,+2.373095986531051e-04,
         +2.112549021673434e-06,+0.000000000000000e+00,+3.011521933753973e-12,+7.123499998001703e-13,+3.905030021754204e-12,+4.041196748133103e-12,
         -1.977140670898991e-12,-3.683757877727392e-12,+2.702332344699793e-12])


    CP=([1],
        [1,+4.605219646724851e-03],
        [1,+7.503996593705215e-02,-3.369804330879776e-05],
        [1,+1.594330632310608e-01,+5.455926851328285e-06,-1.017987924916427e-06],
        [1,+2.376684762866712e-01,+1.043222863021468e-03,+6.973747542233899e-09,-2.887168941583899e-08],
        [1,+3.059383921172673e-01,+5.133947126034137e-03,+5.718347691268169e-09,+8.628372024299278e-10,-8.228341358114456e-10],
        [1,+3.646114133425967e-01,+1.300980487858995e-02,+8.140632050155445e-06,+3.221851528731741e-11,+2.698068570982809e-11,-2.362959895465135e-11],
        [1,+4.150053071989489e-01,+2.442496464955731e-02,+8.431045640169181e-05,+5.684779307309607e-12,-3.263419632097286e-13,+7.550043918275798e-13,
         -6.825523484873172e-13],
        [1,+4.584910800980817e-01,+3.873521575044544e-02,+3.530058481233631e-04,+4.605180684640788e-08,-2.177438868754538e-14,-2.203228206544676e-14,
         +2.078022660794820e-14,-1.979867263528280e-14],
        [1,+4.962620322504207e-01,+5.522305994260744e-02,+9.656377281040835e-04,+9.301269698573483e-07,+5.485773699383311e-15,-4.741922299696983e-16,
         -7.001201405464982e-16,+5.734434107923623e-16,-5.760345737120072e-16],
        [1,+5.293003169818323e-01,+7.323713382827429e-02,+2.061933773275513e-03,+6.204548905461306e-06,+2.106595458316921e-10,-2.568523198895626e-17,
         +9.635334995033390e-18,-1.940235490200283e-17,+1.593366298253801e-17,-1.679674229599595e-17],
        [1,+5.583992702728704e-01,+9.223660464436467e-02,+3.748437306713379e-03,+2.430187011817154e-05,+7.828175967824858e-09,+5.193488453417068e-18,
         +7.827116909900827e-19,+5.622786424411479e-19,-5.166866177936789e-19,+4.459506063411541e-19,-4.905970941733402e-19],
        [1,+5.841970784305069e-01,+1.117935849711689e-01,+6.092187829547772e-03,+6.939236850307271e-05,+8.055548797254791e-08,+8.263137602805995e-13,
         +1.318615068739500e-20,+7.946990755245851e-21,+1.803226504318407e-20,-1.360856427587078e-20,+1.256341428209495e-20,-1.434758873111417e-20],
        [1,+6.072083415695648e-01,+1.315796720142855e-01,+9.123374653370483e-03,+1.607748282199251e-04,+4.432579855910418e-07,+5.396838692968800e-11,
         +4.851645163437771e-21,-1.239807201687139e-22,-2.626494464998897e-22,+4.958335811098015e-22,-3.582245521433100e-22,+3.559491423340867e-22,
         -4.200153498576038e-22],
        [1,+6.278502399234935e-01,+1.513483739829103e-01,+1.284193989158117e-02,+3.213093580289857e-04,+1.673810906101129e-06,+8.336389960397660e-10,
         +2.880717906299260e-15,+2.546211631347504e-23,-1.857750099014604e-23,-1.434688188871249e-23,+1.281674852291496e-23,-9.462859896436801e-24,
         +1.013339192509904e-23,-1.230536185250575e-23],
        [1,+6.464632509695172e-01,+1.709182586574902e-01,+1.722510694751733e-02,+5.755108435311175e-04,+4.912747709141593e-06,+6.331435005116124e-09,
         +3.186701446467602e-13,+4.487942630905718e-24,-6.013183510031929e-25,-1.434800025807743e-25,-4.628323714321082e-25,+3.214071168491499e-25,
         -2.512299663926047e-25,+2.896569534308048e-25,-3.607445858574069e-25])

    CO=([+1.000000000000000e+00],
        [+9.623570038978341e-01,-8.765468670318903e-02],
        [+9.101226931099646e-01,+1.824805676469367e-01,+3.995533125060073e-02],
        [+8.506101279274844e-01,+4.260364614523306e-01,-1.032797973020977e-02,-7.651121295159519e-03],
        [+7.972436814837589e-01,+6.006662832516199e-01,+9.515283399373635e-03,+1.434198495143835e-03,+1.542193116869075e-03],
        [+7.512955511375477e-01,+7.234594998213116e-01,+6.637839682364347e-02,+1.961139275334172e-04,-2.211582783428859e-04,-2.947693861213940e-04],
        [+7.117871395236319e-01,+8.102941661208285e-01,+1.399604513983924e-01,+4.005754210229562e-04,-6.554202838284078e-05,+3.666711324056001e-05,
         +5.542891906686485e-05],
        [+6.775636081175063e-01,+8.720025534723474e-01,+2.188429064185567e-01,+7.127736926151096e-03,+2.220494668649159e-05,+1.327307724175259e-05,
         -6.243555633379748e-06,-1.027063194822467e-05],
        [+6.476421721533173e-01,+9.158680833333335e-01,+2.968651147789767e-01,+2.154641125291114e-02,+1.532958502199612e-05,-3.189555733815597e-06,
         -2.481898615020634e-06,+1.080322189741302e-06,+1.883863046563359e-06],
        [+6.212410039282404e-01,+9.468537051584389e-01,+3.708800879335577e-01,+4.302265405128065e-02,+6.228077763904015e-04,+6.667811134655501e-07,
         +4.127849028415000e-07,+4.524147501189836e-07,-1.884944216882671e-07,-3.428282107115735e-07],
        [+5.977472638192826e-01,+9.684019569318768e-01,+4.394388409641376e-01,+7.025448036639782e-02,+2.608972301400094e-03,+5.552989991589061e-07,
         -4.196062624508686e-09,-5.385100269523688e-08,-8.153275576805928e-08,+3.304200721827606e-08,+6.200255109382937e-08],
        [+5.766796716480526e-01,+9.829442237604772e-01,+5.020357046307879e-01,+1.018179697173380e-01,+6.540345416710769e-03,+4.756952866121665e-05,
         +5.089169343259292e-09,-8.872848204382457e-09,+7.128642451939746e-09,+1.460034175832762e-08,-5.806292514515928e-09,-1.115732090309261e-08],
        [+5.576572656998297e-01,+9.922271896007677e-01,+5.586776137115197e-01,+1.364042232734678e-01,+1.273827068703519e-02,+2.688760930684044e-04,
         +1.940752294511806e-08,+3.445753326482956e-09,+2.339211996258690e-09,-9.527702066379118e-10,-2.603601127542738e-09,+1.021501804838346e-09,
         +1.999437402895323e-09],
        [+5.403754612885416e-01,+9.975253420484818e-01,+6.096384134139301e-01,+1.729035349004602e-01,+2.131182377683984e-02,+8.337184440722163e-04,
         +3.299979819721097e-06,-4.647080307976002e-10,-6.643590722317443e-10,-4.763787357029944e-10,+1.277051805015050e-10,+4.628542401914760e-10,
         -1.797836319412636e-10,-3.570604982407735e-10],
        [+5.245881732125043e-01,+9.997817888994765e-01,+6.553190105584938e-01,+2.104181834166170e-01,+3.221246838477771e-02,+1.919009659995399e-03,
         +2.463546362233918e-05,+6.613864820110626e-10,+1.282043947217749e-10,+9.776989499930870e-11,+8.920559856480144e-11,-1.703871970591045e-11,
         -8.208176899117022e-11,+3.163933635559848e-11,+6.357500710461819e-11])


    I=np.float64(0)
    z=np.zeros((N,),dtype=np.float64)
    if mpirank==0:
        
        #if P != 0: x=(np.arange(0,M+1)+0.5)/(M+2)
        #else:      x=(np.arange(0,M+1))/(M+2)
        #if   P ==  1: C=np.cos(math.pi*x/2)**(2*(M+1))

        if   P ==  1: C=np.asarray(CE[M])
        elif P ==  0: C=np.asarray(CP[M])
        elif P == -1: C=np.asarray(CO[M])

        
        #if   P ==  0: I=(2*C.sum()-C[0])/N
        #else: I=C.sum()/N
        #C=C/I
        z[0:M+1]=C
        if   P== 0: z[N-M  :N] += np.flip(C[1:M+1])
        I=z.sum()/N

    I   = comm.bcast(I,root=0)    
    z   = comm.bcast(z,root=0)    
    return(z,I)

# odd constant
# returns a function that is 1 over most of the domain but has odd symmetry 
# Coefficients are calculating by minimising int((1-f)^2,x=1/2..1) with f(1)=1
# Number of coefficients in sin series is M+1 with M+1 interior points
def ffta_constant3(N,M):
    CO=([+7.071067811865475e-01],
        [+5.069911752244420e-01,+1.048187275370639e+00],
        [+4.513737218293190e-01,+9.713190591669934e-01,+9.933320266762983e-01],
        [+3.960420355148196e-01,+9.098328345958588e-01,+1.003015691540891e+00,+1.000836356098086e+00],
        [+3.604575554072329e-01,+8.604157803841613e-01,+9.974768501181859e-01,+9.992655503860745e-01,+1.000048756900463e+00],
        [+3.317055294062269e-01,+8.163062581910430e-01,+9.851402504716147e-01,+1.000116138744975e+00,+1.000167642025735e+00,+9.999149610632119e-01],
        [+3.092861337408186e-01,+7.786530908537965e-01,+9.696235628970070e-01,+9.997794691239811e-01,+9.999596115118618e-01,+9.999740620797646e-01,
         +1.000038678487665e+00],
        [+2.907049876077845e-01,+7.453510139143181e-01,+9.525038455608006e-01,+9.977086606984226e-01,+9.999972288764881e-01,+1.000014521664447e+00,
         +1.000000072984290e+00,+9.999861515257748e-01],
        [+2.751647700465852e-01,+7.159791282323538e-01,+9.348753412741138e-01,+9.940337240329432e-01,+9.999806620731560e-01,+1.000000060009068e+00,
         +9.999959462924172e-01,+1.000001893061538e+00,+1.000004293393153e+00],
        [+2.618508202125438e-01,+6.897174307208807e-01,+9.172087232659815e-01,+9.890208228834313e-01,+9.996638710461835e-01,+9.999988891632734e-01,
         +1.000000500936609e+00,+1.000000820595772e+00,+9.999990786158837e-01,+9.999987892704644e-01],
        [+2.503052886982728e-01,+6.661188233728717e-01,+8.998487377736846e-01,+9.829547909292746e-01,+9.989109397212529e-01,+9.999982990037833e-01,
         +1.000000277547367e+00,+9.999997256897890e-01,+9.999999246904569e-01,+1.000000304845309e+00,+1.000000311580587e+00],
        [+2.401586029187565e-01,+6.447491651651891e-01,+8.829535729124889e-01,+9.760764997707723e-01,+9.976649053596686e-01,+9.999525616847619e-01,
         +9.999998599724389e-01,+9.999999584254093e-01,+1.000000094771199e+00,+9.999999719484480e-01,+9.999999207874259e-01,+9.999999276521558e-01],
        [+2.311537681786137e-01,+6.252924574701847e-01,+8.666241808508276e-01,+9.685915673248549e-01,+9.959185063936782e-01,+9.998120045314867e-01,
         +9.999998500627408e-01,+1.000000035799624e+00,+9.999999977076431e-01,+9.999999764339713e-01,+1.000000019251636e+00,+1.000000015631827e+00,
         +1.000000014398830e+00],
        [+2.230900835179195e-01,+6.074794637961929e-01,+8.509017873328270e-01,+9.606634806122486e-01,+9.936923975076562e-01,+9.995345900519423e-01,
         +9.999935087593218e-01,+9.999999869286569e-01,+9.999999914176127e-01,+1.000000004648285e+00,+1.000000003644772e+00,+9.999999929744306e-01,
         +9.999999984052742e-01,+9.999999979641842e-01],
        [+2.158150464352705e-01,+5.910950308330348e-01,+8.358020989892738e-01,+9.524245035431750e-01,+9.910226479063029e-01,+9.990891562187399e-01,
         +9.999689594801942e-01,+9.999999867654165e-01,+1.000000002922487e+00,+1.000000001689564e+00,+9.999999979125738e-01,+1.000000000153246e+00,
         +1.000000001868211e+00,+9.999999995529008e-01,+9.999999999465629e-01])

    

    z=np.ones((N,),dtype=np.float64)
    C=np.asarray(CO[M])
    z[0:M+1]   = C
    z[N-M-1:N] = np.flip(C)
    return(z)
       


    
    
# f(x) ~ H(x+epsilon)
# Suitable for forcing and limiting since it evaluates very close to 1
def ffta_heaviside3(N,M):
    CC=([],
        [+1.286071037125326e-02],
        [+6.271796581543043e-02,+1.412260274952829e-04],
        [+1.051393254177300e-01,+3.772649213755421e-03,+8.386548696467935e-07],
        [+1.382599408610645e-01,+1.208441978752237e-02,+1.351478736273088e-04,+3.106720514419145e-09],
        [+1.645754406200722e-01,+2.296004565037766e-02,+8.843638899712638e-04,+3.227951385299871e-06,+7.853949447206996e-12],
        [+1.860132062293014e-01,+3.486594876342259e-02,+2.554693334809252e-03,+4.509148669744680e-05,+5.519169034475735e-08,+1.440606003001717e-14],
        [+2.038714332883750e-01,+4.695905167095065e-02,+5.154064956947938e-03,+2.047701343516511e-04,+1.699500838061753e-06,+7.096840508830224e-10,
         +2.004200411621180e-17],
        [+2.190284007008206e-01,+5.880487222333651e-02,+8.550833490439235e-03,+5.638973167245264e-04,+1.243171251523158e-05,+4.937880081919817e-08,
         +7.116977660993100e-12,+2.187085326890424e-20],
        [+2.320940912026030e-01,+7.019129726901110e-02,+1.257725800271547e-02,+1.178788435911466e-03,+4.766175858547887e-05,+5.929189233103336e-07,
         +1.141280313246707e-09,+5.724324349587762e-14,+1.921975271202317e-23],
        [+2.435043886881590e-01,+8.102574939591448e-02,+1.707421367470202e-02,+2.076425544617973e-03,+1.276530245214844e-04,+3.213773106482325e-06,
         +2.284140271406407e-08,+2.150267147717017e-11,+3.775842681739444e-16,+1.389252981206190e-26],
        [+2.535790923883089e-01,+9.127974360672412e-02,+2.190616979028266e-02,+3.260389886616910e-03,+2.732040512604488e-04,+1.117254106406786e-05,
         +1.772287621128457e-07,+7.264460133731700e-10,+3.367850475479932e-13,+2.079869076996595e-18,+8.404703590225977e-30],
        [+2.625583936484906e-01,+1.009587853834970e-01,+2.696346188277533e-02,+4.718514164759101e-03,+5.025244136081099e-04,+2.938287771657365e-05,
         +8.082283824416923e-07,+8.153061456092488e-09,+1.941363232628282e-11,+4.456153479253586e-15,+9.712489895479419e-21,+4.318082655597452e-33],
        [+2.706264428753874e-01,+1.100858253780487e-01,+3.215980575325566e-02,+6.429322992514143e-03,+8.296886134057359e-04,+6.394796728219044e-05,
         +2.636068995541027e-06,+4.920612116753494e-08,+3.179630629399270e-10,+4.423517351538438e-13,+5.048449484959188e-17,+3.894217534257294e-23,
         +1.907441568596978e-36],
        [+2.779269974308732e-01,+1.186920798971257e-01,+3.742846882685894e-02,+8.366687089818590e-03,+1.264263080307367e-03,+1.215351573311590e-04,
         +6.843761462523374e-06,+2.005696909162306e-07,+2.558884509166287e-09,+1.065418513728226e-11,+8.699847329519643e-15,+4.953476155401250e-19,
         +1.355288772136417e-25,+7.321724580392536e-40],
        [+2.845740803488343e-01,+1.268119266885863e-01,+4.271846117637584e-02,+1.050291410681536e-02,+1.811615572443511e-03,+2.087906933956335e-04,
         +1.508252423864208e-05,+6.254300635971838e-07,+1.312141811525791e-08,+1.150844331689248e-10,+3.102173512152950e-13,+1.492389002246289e-16,
         +4.250964068247833e-21,+4.132931433930561e-28,+2.464801565351292e-43])
    
    C=np.asarray(CC[M]) 
    z=np.zeros((N,),dtype=np.float64)
    z[  0:  M+0] = 1-np.flip(C)
    z[M+0:  M+1] = 0.5
    z[M+1:2*M+1] = C
    
    return(z)

# Band limited signals everything
# is built around the odd heaviside function and its derivative
# k sets the highest wavenumber of the function
# M sets the order of the highest vanishing derivative
# Heaviside function and derivative
# k is maximum wavenumber, M is order
# x are evaluation points
# H(x) = 1 x>>0
#      = 0 x<<0

def ffta_heaviside(x,k,M):
    xx=x;
    xc = (2*M+1)/k/2
    x  = np.minimum(math.pi/2,np.maximum(-math.pi/2,k*x/(2*M+1)))
    If = np.maximum(0,(x+math.pi)*xc)+np.maximum(0,xx-xc*math.pi)
    f  = 1/2+np.zeros_like(x)
    If0=0
    df = np.zeros_like(x)
    for n in range(M+1):
        sf = k*(2*n+1)/(2*M+1)
        c=2/math.pi*np.exp(2*sp.gammaln(M+3/2)-sp.gammaln(M+2+n)-sp.gammaln(M+1-n))/(2*n+1)
        If  += -c/sf*np.cos((2*n+1)*x)
        f   +=  c   *np.sin((2*n+1)*x)
        df  +=  c*sf*np.cos((2*n+1)*x)
    df0 = k*np.exp(sp.gammaln(2*M+1)-M*np.log(4)-2*sp.gammaln(M+1))
    return(f,df,df0,If)



def ffta_pheaviside(x,y,P,k,M,L):
    if   P ==  0:  f,df,df0=ffta_heaviside(np.mod(x-y+L/2,L)-L/2,k,M)
    else:
        af,adf,adf0=ffta_pheaviside( x,y,0,k,M,2*L)
        bf,bdf,bdf0=ffta_pheaviside(-x,y,0,k,M,2*L)
        if   P ==  1:
            f   = af+bf
            df  = adf+bdf
            df0 = adf0+bdf0
        else:
            f   = af-bf
            df  = adf-bdf
            df0 = adf0-bdf0
    return(f,df,df0)



def ffta_delta(x,y,P,k,M,L):
    if   P ==  0:  df=ffta_heaviside(np.mod(x-y+L/2,L)-L/2,k,M)[1]
    elif P ==  1:  df=ffta_delta(x,y,0,k,M,2*L)+ffta_delta(-x,y,0,k,M,2*L)
    elif P == -1:  df=ffta_delta(x,y,0,k,M,2*L)-ffta_delta(-x,y,0,k,M,2*L)
    return(df)

def ffta_interval(x,y1,y2,P,k,M,L):
    if   P ==  0:
        if y1>y2: f=1-ffta_interval(x,y2,y1,P,k,M,L)
        else:
            x=np.mod(x-y1+L/2,L)-L/2
            f=ffta_heaviside(x,k,M)[0]-ffta_heaviside(x+y1-y2,k,M)[0]
    elif P ==  1: f=ffta_interval(x,y1,y2,0,k,M,2*L)+ffta_interval(-x,y1,y2,0,k,M,2*L)
    elif P == -1: f=ffta_interval(x,y1,y2,0,k,M,2*L)-ffta_interval(-x,y1,y2,0,k,M,2*L)
    return(f)

# return the necessary radius to give the correct radius for an offset circle
def get_radius(r,Y,Z,Py,Pz):
    if Py==0 and Pz==0: return(r)
    
    A = math.pi*r**2
    if Py !=0: A=A/2
    if Pz !=0: A=A/2
    f = lambda R: circle_area(R,-Y,-Z,Py,Pz)-A
    
    #R=secant(f,A,r/2,r+np.sqrt(Y**2+Z**2))
    r1=r/2
    r2=r+np.sqrt(Y**2+Z**2)
    if False:
        f1=f(r1)
        f2=f(r2)
        print('get_radius: r1={:7.4f}, r2={:7.4f}, f1={:7.4f}, f2={:7.4f}'.format(r1,r2,f1,f2))   
    R=opt.brentq(f,r1,r2)
    return(R)

def circle_area(r,Y,Z,Py,Pz):
    #print('circle_area: r={:6.3f},  Y={:6.3f}, Z={:6.3f}, Py={:3d}, Pz={:3d}'.format(r,Y,Z,Py,Pz)) 
    if   Py == 0 and Pz == 0: A=r**2*math.pi
    elif Py == 0 and Pz != 0: A=r**2*jasin(Z/r)
    elif Py != 0 and Pz == 0: A=r**2*jasin(Y/r)
    else:                   A=r**2*jacircle(Y/r,Z/r)
    return(A)


def jacircle(p,q):
    if np.isscalar(p) and np.isscalar(q):
        if p**2+q**2>1:
            if   p>=0 and q>=0: A=0
            elif p< 0 and q>=0: A=jasin(q)
            elif p>=0 and q< 0: A=jasin(p)
            else:               A= math.pi-jasin(-p)-jasin(-q)
            
        else:
            A=-math.pi/4+p*q+jasin(p)/2+jasin(q)/2
    else:
        p,q=np.broadcast_arrays(p, q)
        A=np.zeros_like(p)
        for j in range(p.size) : A.flat[j]=jacircle(p.flat[j],q.flat[j])
    return(A)
    
# Return area of a circle with radius to the right of the line x=s  
def jasin(s):
    if np.isscalar(s):
        if   s>=1: A=0
        elif s>-1: A=math.pi/2-s*np.sqrt(1-s**2)-np.arcsin(s)
        else:      A=math.pi
    else:
        A=np.zeros_like(s)
        for j in range(s.size) : A.flat[j]=jasin(s.flat[j])
    return(A)
    

# solve f(x)=a using secant method (x1) and f(x2) must have different signs
# the function most be monotonic and continuous
# this latter restriction could be removed with a switch to a midpoint strategy if there are
# are multiple solutions
# usually takes 6 iteration
def secant(f,a,x1,x2):
    f1=f(x1)-a
    f2=f(x2)-a
    tol=1e-12
    if f1==0: return(x1)
    if f2==0: return(x2)
    x=x1
    f3=f1
    k=0
    while x1 != x2 and k < 100:
        k=k+1
        if f1*f2>=0:
            print('secant: f1*f2>=0 k={:3d}, x={:7.4f}, x1={:7.4f}, x2={:7.4f}, f={:7.1e}, f1={:7.1e}, f2={:7.1e}'.format(k,x,x1,x2,f3,f1,f2))
            quit(1)
            
        x = (x1*f2-x2*f1)/(f2-f1)
        f3=f(x)-a

        if abs(f3)<tol: break
        if f3*f1>0:
            f1=f3
            dx=np.abs(x1-x)
            x1=x
            #if dx<tol: x2=(x1+x2)/2
        else:
            f2=f3
            x2=x
            dx=np.abs(x2-x)
            #if dx<tol: x1=(x1+x2)/2
        print('k={:3d}, x={:7.4f}, x1={:7.4f}, x2={:7.4f}, dx={:7.1e}, f={:7.1e}, f1={:7.1e}, f2={:7.1e}'.format(k,x,x1,x2,dx,f3,f1,f2))
    return(x)

def add_zm_noise(f,inoise,s,V,rand):
    if f is None: return()
    mf1=jint(f)/V
    f.set_scales(1)
    f['g']
    srange(s + ' range before [{:6.3f},{:6.3f}]',f)
    f['g']+=+inoise*rand.standard_normal(f['g'].shape)
    f.set_scales(1)
    srange(s+ ' range after [{:6.3f},{:6.3f}]',f)
    mf2=jint(f)/V
    f['g']+=mf1-mf2
    mf3=jint(f)/V
    srange(s+ ' int(f) {:8.1e}',mf3-mf1)
    f.set_scales(AA)
    
def set_parity(f,Tx,Ty,Tz,Px,Py,Pz):
    if Tx is not None:
        if Tx == "SinCos":  f.meta['x']['parity'] = Px
        else:               f.meta['x']['parity'] = 0

    if Ty is not None:
        if Ty == "SinCos":  f.meta['y']['parity'] = Py
        else:               f.meta['y']['parity'] = 0

    if Tz is not None:
        if Tz == "SinCos":  f.meta['z']['parity'] = Pz
        else:               f.meta['z']['parity'] = 0


def cheb_coord(N,I=np.array([-1,1],dtype=np.float64)):
    x=I[0]+1+(I[1]-I[0])/2*np.cos((1+2*np.arange(N))*math.pi/(2*N))
    return(x)

def cheb_resample(f,N):
    minf=f.min()
    maxf=f.max()
    #logger.info('minf {} maxf {}'.format(minf,maxf))
    M=f.shape[0]
    x=cheb_coord(M)
    c=np.polynomial.chebyshev.chebfit(x,f,M-1)
    x=cheb_coord(N)
    f=np.polynomial.chebyshev.chebval(x, c)
    f=np.maximum(minf,np.minimum(maxf,f))
    #logger.info('ming {} maxg {}'.format(f.min(),g.max()))
    return(f)

def pchipi_resample(f,N,axis=0): # pchip resampling at interior points to size N along 
    minf=f.min()
    maxf=f.max()
    M=f.shape[axis]
    x=(0.5+np.arange(M))/M
    y=(0.5+np.arange(N))/N
    f=it.pchip_interpolate(x,f,y,axis=axis)
    f=np.maximum(minf,np.minimum(maxf,f))
    return(f)

#scipy.fftpack.dst
def fft_resample(f,N,parity=0,axis=0): # reall fft resampling with clipping:
    shape=f.shape
    M=shape[axis]
    if parity == -1:
        f = fft.dst( f, 2, M, axis=axis)
        f = fft.idst(f, 2, N, axis=axis)/(2*M)
    if parity ==  0:
        f = fft.rfft( f, M, axis=axis)
        f = fft.irfft(f, N, axis=axis)*N/M
    if parity ==  1:
        f = fft.dct( f, 2, M, axis=axis)
        f = fft.idct(f, 2, N, axis=axis)/(2*M)
    logger.info('jutil.fft_resample: M={}, N={}, Parity={}'.format(M,N,parity))
    return(f)

def fft_diff(f,dx,parity=0,axis=0): # fft derivative result should be multiplied by dx
    
    shape=f.shape
    M=shape[axis]
    N=np.int(M/2)
    if parity ==  0:
        w = 2j*math.pi*(np.mod(np.arange(M)+N,M)-N)/(M*dx)
        if M==2*N: w[N]=0
        #w=2j*math.pi*fft.fftfreq(M)/dx
        f = fft.fft( f, axis=axis)
        f = fft.ifft(w*f, axis=axis).real
    else:
        if parity == -1: f = np.concatenate((f,-np.flip(f,axis)),axis=axis)
        else:            f = np.concatenate((f,np.flip(f,axis)),axis=axis)
        f = jfft_diff(f,dx,0,axis)
        if   axis==0: f=f[0:M,...]
        elif axis==1: f=f[:,0:M,...]
        elif axis==2: f=f[:,:,0:M,...]
    return(f)

#ju.cheb_test(6)
def cheb_test(N):
    f=np.array([0, 1, 2, 3, 2, 1, 0])
    g=cheb_resample(f,N)
    print(g)



def log_param(p,nm,f='{:8.5f}'):
    for a in nm:
        if a in p:
            if p[a] is not None:
                logger.info(('{:8s} '+f).format(a+':',p[a]))


def gatheryz(gf,lf,ys,zs):
    if mpirank==0:
        gf[ys[0],zs[0]]=lf
        for k in range(1,mpisize): gf[ys[k],zs[k]]=comm.recv(source=k)
    else: comm.send(lf,dest=0)
    comm.barrier()
    #comm.Gatherv(sendbuf=lrur,recvbuf=(rur,yyzzcounts,yyzzstart,MPI.DOUBLE),root=0)
    #(buffer, counts, displacements, type)

def jtypes():
    return(('scalar','ndarray','x','fx','yz','fyz','xyz','fxyz'))

def jtype(x):
    if np.isscalar(x):           return('scalar')
    if isinstance(x,np.ndarray): return('ndarray')
    if isinstance(x,co.field.Field):
        #print('dir(x.domain) {}'.format(dir(x.domain)))
        #print('vars(x.domain) {}'.format(vars(x.domain)))
        gsize=x.domain.global_grid_shape(scales=1)
        if len(gsize)==3:
            if np.all(gsize==[Nx,  Ny,  Nz]):   return('fxyz')
            if np.all(gsize==[Nxx, Nyy, Nzz]):  return('fxyz')
            if np.all(gsize==[1,   Ny,  Nz]):   return('fyz')
            if np.all(gsize==[1,   Nyy, Nzz]):  return('fyz')
            if np.all(gsize==[Nx,  1,   1]):    return('fx')
            if np.all(gsize==[Nxx, 1,   1]):    return('fx')
        elif len(gsize)==2:
            if np.all(gsize==[Nx,  Nz]):   return('fxyz')
            if np.all(gsize==[Nxx, Nzz]):  return('fxyz')
            if np.all(gsize==[1,   Nz]):   return('fz')
            if np.all(gsize==[1,   Nzz]):  return('fz')
            if np.all(gsize==[Nx,  1]):    return('fx')
            if np.all(gsize==[Nxx, 1]):    return('fx')
        elif len(gsize)==1:
            if np.all(gsize==[Nx]):     return('fx')
            if np.all(gsize==[Nxx]):    return('fx')
 
        #logger.info('vars(x) {}'.format(vars(x)))
    return('Unknown')


def barrier_debug(s=' '):
    logger.error('barrier_debug: {}rank {} barrier 1'.format(s,mpirank))
    Barrier()
    logger.error('barrier_debug: {}rank {} barrier 2'.format(s,mpirank))


def print_v(a,x,s=''):
    if mpirank==0: typ=jtype(x)
    else:          typ=''
    typ = comm.bcast(typ, root=0)
    s= '{} {:6s} {:7s} '.format(s,a,typ)
    if   typ=='fx':      gx.print_v(x,s)
    elif typ=='fyz':     srange(     s + '[{:11.5f},{:11.5f}]' + ' {} {}'.format(x.domain.dist.grid_layout.global_shape(1)[1:],parity_string(x,('y','z'))),x)
    elif typ=='fxyz':    srange(     s + '[{:11.5f},{:11.5f}]' + ' {} {}'.format(x.domain.dist.grid_layout.global_shape(1),    parity_string(x)),x)
    elif typ=='scalar' and mpirank==0:  logger.info(s + ' {:11.5f}'.format(x))
    elif typ=='ndarray': logger.info(s + '[{:11.5f},{:11.5f}] {}'.format(x.min(),x.max(),x.shape))


# find zero in array y using linear interpolation and the distance away from that     
def find_zero(y,z=None):
    n=y.shape[0]
    x=-1
    for j in range(n-1):
        if y[j]*y[j+1]<0:
            x = y[j]/(y[j]-y[j+1]) #(1-x)*y[j]+x*y[j+1]=0   
            break
    #j=np.int64(np.floor(x))
    #r=x-j
    #logger.info('find_zero r {} {} {}'.format(y[j],(1-r)*y[j]+r*y[j+1],y[j+1]))
    if z is None: return(x+j)
    else: return(z.flatten()[j]*(1-x)+z.flatten()[j+1]*x)

    # Forcing function for density is to choose an error function that is positive wherever the inlet velocity
    # but with an amount of smoothing that depends on the grid size
    # Typically this produce inlet profiles that are too sharp
def db_forcing(u):
    n=u.shape[0]
    r=find_zero(u)
    #(1+sp.erf(Nz/5*state_ev['divz']/(state_ev['divz'].max()-state_ev['divz'].min())))/2
    s=(1+sp.erf((r-np.arange(n))/3))/2
    #i=np.int64(np.floor(r))
    #f=slice(i,i+2)
    #logger.info('f {}'.format(f))
    #logger.info('u {}'.format(u[f]))
    #logger.info('s {}'.format(s[f]-0.5))
    #quit(1)
    return(s)


# Advance f using logistic equation so that
# int (f-fmin)/(fmax-fmin) dz  = h
# f is field
# h is target h
# w is integration weights
# H is  height w.sum()
# f0 minimum field value
# f1 maximum field value
# scale so that int f = -1 and maximum value is unaltered

# equation for f rescaled to [0,1] is
#
# df    h-g         
# -- =  ---- f (1-f)
# dt     T
#
# g = int(f)
#
# h is target height
#
def logistic_h(f,h,w,H,f0=None,f1=None,scale=False,a=0):
    if f0 is None:f0 = f.min()
    if f1 is None:f1 = f.max()
    g   = (f-f0)/(f1-f0)
    hg  = (w*g).sum()                   # calculate hu
    hw  = (w*g*(1-g)).sum()             # calculate hu
    if a==0:  e  = (hg-h)/hw            # Linear solve to exactly match h
    else:     e  = np.exp((hg-h)/a)-1   # step forward a = T/dt
    e   = max(-1/2,min(1/2,e))
    g   = g/(1+e*(1-g))                 # Logistic equation solution
      
    if scale:
        hg2 = (w*g).sum()                   # calculate hu
        f = f1-(g-1)*(f1+1)/(hg2/H-1)       # Set int(g)/H=-1 and max(g)=maxg
    else:
        f = f0 + (f1-f0)*g
    if DED_DEBUG_LOGH:
        hg2 = (w*g).sum()                   
        ming2=f.min()
        maxg2=f.max()
        Ig=(w*f).sum()/H
        chQ=(h-hg2)/(h-hg)
        logger.info('logh: h {:4.2f} hQ {:8.6f} {:8.6f} chQ {:8.6f} min {:5.2f} {:5.2f} max {:5.2f} {:5.2f} I {:5.2f} a {:5.1e} e {:5.1e}'.format(h,hg,hg2,chQ,f0,ming2,f1,maxg2,Ig,a,e))
    return(f)

# Perform largest stable diffusion step according if function is non-monotonic
#
# D = int(|df|)/(max(f)-min(f))-1 diffusion rate
#
# D=0 if function is monotonic
#
# w are integration weights so that int(f dw) is conserved max(dw)<=1
def tvd_diffusion(f,w=1):
    df = np.diff(f)
    T1 = np.abs(df).sum()
    D  = T1/(f.max()-f.min())-1     # D zero if f is monotonic
    D  = min(1/4,max(0,10*D))
    f += D*w*np.diff(df,prepend=0,append=0)
    if DED_DEBUG_TVDD:
        T2  = np.abs(np.diff(f)).sum()    # D zero if f is monotonic
        logger.info('tvd_diffusion D:{:} T2/T1:{:}'.format(D,T2/T1))
    return(f)
