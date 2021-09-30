import numpy as np 
import copy as copy
from pathlib import Path
from dedalus import public as de
import jutil as ju
import h5py
import time as time
import ded_types as dp
import scipy.interpolate as it
from fractions import Fraction
debug_single=False
debug_flow=False
debug_avrg=False
debug_var=False
debug_slices=False


def jita(jrec,javrg,avrg,i0,i1,i2,dtnm,fnm,p):
    if not dtnm in p: return()
    dt=p[dtnm]
    if dt is None: return()
    if dt>0: javrg.ji.append(jit(jrec,javrg, avrg,i0,i1,i2, dt,fnm, p))

#self.scomm=new_comm(ju.comm,self.color,ju.mpirank)
def new_comm(ocomm,color,rank):
    color=np.int64(color)
    mx=ocomm.allreduce(np.int(color), op=ju.MPI.MAX)
    mn=ocomm.allreduce(np.int(color), op=ju.MPI.MIN)
    if mn==mx: ncomm=ocomm
    else: ncomm=ju.comm.Split(color,rank)
    #print('jrec.new_comm rank {} color {}'.format(rank,color))
    return(ncomm)
  
def get_wn_fn(d,nm,reset):
    if ju.mpirank==0:
        if not d.is_dir(): d.mkdir(parents=True)
    ju.Barrier()
    #ju.logger.info('get_wn: {} {}'.format(nm,d))
    wn = -1
    if not reset and ju.mpirank==0:
        p = sorted(d.glob('{}-*.hdf5'.format(nm)),reverse=True)
        for a in p:
            try:
                f=h5py.File(str(a),'r')
                #for b in f: print('   b: {}'.format(f[b]))
                if 'wn' in f:
                    wn = f['wn'][()]
                    #ju.logger.info('get_wn: {}, file: {}, wn: {}'.format(nm,a,wn))
                    break
                elif 'scales' in f:
                    if 'write_number' in f['scales']:
                        wn = f['scales']['write_number'][()]
                        #ju.logger.info('get_wn: {}, file: {}, write_number: {}'.format(nm,a,wn))
                        break
                else: ju.logger.info('get_wn: {}, file: {} no wn'.format(nm,a))
            except: ju.logger.info('get_wn: {}, file: {} failed to open'.format(nm,a))
    wn +=1
    wn = ju.comm.bcast(wn, root=0)
    #ju.logger.info('get_wn: {}, wn: {}'.format(nm,wn))
    return(d,wn)


class jsolve:
    def __init__(self,v,p):
        self.iteration=np.int(0)
        self.write_number=np.int(0)
        self.sim_time=np.float(0)
        self.state={}
        for a in v: self.state[a]=None
        self.state_ev={}
        if 'noised' in p:
            if 'x' in p['noised']: self.state_ev['dxnoise']=None
            if 'y' in p['noised']: self.state_ev['dynoise']=None
            if 'z' in p['noised']: self.state_ev['dznoise']=None
        if 'pmms'  in p:
            if p['pmss'] >0:
                self.state_ev['SY']=None
                self.state_ev['SZ']=None

        if p['sType']=='pm':
            if p['Forcing']==7:
                self.state_ev['vf']=None
            if p['Forcing']==6:
                self.state_ev['qvx']=None
                self.state_ev['Fx']=None


        return(s)

def check_rec_dt(solver,wn,dt):
    if dt>0: return(solver.sim_time  >= wn*dt)
    else:    return(solver.iteration >= wn)

def get_wn(solver,dt):
    
    if np.isscalar(dt):
        if dt>0: wn = np.int64(np.ceil(solver.sim_time/dt))
        else:    wn = np.int64(solver.iteration)
    else:
        wn=np.array(np.where(dt==0,solver.iteration,np.ceil(solver.sim_time/dt)),dtype=np.int64)
    return(wn)

class sdict:
    def __init__(self,aa,reset,ddict,nm,dt):
        self.debug=False
        self.var=ddict
        self.vnm=dict()
        if ju.mpirank==0:
            for a in ju.jtypes():  self.vnm[a] = list()
            for a in ddict.keys(): self.vnm[ju.jtype(ddict[a])].append(a)
        self.vnm=ju.comm.bcast(self.vnm,root=0)
        n=0
        for a in self.vnm.keys():
            if len(self.vnm[a])>0: n=n+1
        #for a in self.vnm: ju.logger.info('self.vnm[{}] {}'.format(a,self.vnm[a]))
        
        self.empty=n==0
        if self.empty:
            ju.logger.info('sdict: {} is empty'.format(nm))
            return
        if self.debug:
            ju.logger.info('sdict: {}'.format(nm))
            for a in self.vnm.keys():
                k=self.vnm[a]
                if len(k)>0: ju.logger.info('           {:7s}: {}'.format(a,k))
            
        self.dt       = dt
        self.solver   = aa.solver
        self.wn       = 0
        self.nm       = nm
        self.ddir     = Path(aa.dir)
        self.sim_time = Path(aa.dir)

        if ju.mpirank==0:
            self.wn=get_wn(self.solver,self.dt)
            d = aa.dir / self.nm
            if not d.is_dir(): d.mkdir(parents=True)
        self.wn = ju.comm.bcast(self.wn, root=0)
        ju.logger.info('jrec: nm: {:5s} dt: {:6.3f} num: {:5d}'.format(self.nm,self.dt,self.wn))

    def update(self,dt):
        if self.empty: return
        if not check_rec_dt(self.solver,self.wn,self.dt): return
        fn = str(self.ddir / self.nm / '{}-{:05d}.hdf5'.format(self.nm,self.wn))
        self.wn+=1
        if ju.mpirank==0:
            with h5py.File(fn, 'w') as f:  
                f.create_dataset('wn',data=self.wn)
                f.create_dataset('t', data=self.solver.sim_time)
                f.create_dataset('dt',data=dt)
                f.create_dataset('i', data=self.solver.iteration)
                for a in self.vnm['x']:       f.create_dataset(a,data=self.var[a])
                for a in self.vnm['scalar']:  f.create_dataset(a,data=self.var[a])
                for a in self.vnm['ndarray']: f.create_dataset(a,data=self.var[a])
                
        #if self.nm=='avar':
        #for a in self.ddict.keys():
        #    for b in self.ddict[a].keys(): ju.logger.info('jrec update nm=avar {} {} {}'.format(a,b,self.ddict[a][b]))
        #ju.srange('avar[dwdz] [{:6.3f},{:6.3f}]',self.var['dwdz'])
        #ju.srange('avar[bp]   [{:6.3f},{:6.3f}]',self.var['bp'])

        ju.gyz.write(fn,self.var,'a',self.vnm['fyz'])
        ju.gx.write( fn,self.var,'a',self.vnm['fx'])
        
class single_field:
    def __init__(self,a,reset,nm,dt):
        self.n         = len(nm)
        if self.n==0: return
        
        self.vnm = copy.copy(nm)
        self.nm  = copy.copy(nm)
            
        self.dt        = np.array(dt)
        self.dim       = a.szJ.dim
        self.solver    = a.solver
        self.nvar      = a.szJ.lsize.prod()  # How many elements  
        self.szA       = a.szJ
        self.AA        = a.AA
        self.AAJ       = a.AAJ
        self.reduce    = a.reduce
        
        self.color= self.nvar>0
        self.scomm=new_comm(ju.comm,self.color, ju.mpirank)
        
        self.wn=np.zeros((self.n),dtype=np.int64)

        self.ddir=Path(a.dir)

        if ju.mpirank==0:
            self.wn=get_wn(self.solver,self.dt)
            for nm in self.nm:
                d = a.dir / nm
                if not d.is_dir(): d.mkdir(parents=True)
                if debug_single: ju.logger.info('Dir: {}'.format(d))

        # Make sure these exactly the same on every node
        self.wn = ju.comm.bcast(self.wn, root=0)
        self.dt = ju.comm.bcast(self.dt, root=0)
        
        for j in range(0,self.n): ju.logger.info('jrec: nm: {:5s} dt: {:6.3f} num: {:5d}'.format(self.nm[j],self.dt[j],self.wn[j]))
        
        #if debug_single:
        #    for j in range(0,self.n):
        #        print('{:4d} n{:2d}: nm:{:2s} dt:{:5.2f} num:{:5d}'.format(ju.mpirank,j,self.nm[j],self.dt[j],self.wn[j]))

    def getvar(self):
        #for a in var: ju.logger.error('single_field.update:{:4d} {} {}'.format(ju.mpirank,a,var[a].shape))
        v=set()
        for j in range(self.n):
            if self.solver.sim_time >= self.wn[j]*self.dt[j]: v.add(self.vnm[j])
        
            #if len(v)>0: ju.logger.info('single_field.getvar: {}'.format(v))
        return(v)

    def update(self,dt,var):
        for j in range(self.n):
            if check_rec_dt(self.solver,self.wn[j],self.dt[j]): self.save(j,dt,self.reduce(var[self.vnm[j]],self.vnm[j]))

    def save(self,j,dt,xx):
        if self.color:
            if debug_flow: ju.logger.info('jrec.sf.save: j: {:1d}, dt {:5.3f}, nm: {:s}'.format(j,dt,self.nm[j]))
            nm=self.nm[j]
            wn=self.wn[j]
            fn = str(self.ddir / nm / '{}-{:05d}.hdf5'.format(nm,wn))
            if debug_single: ju.logger.info('single_field.save: {:4d} {} {} {} {} {}'.format(ju.mpirank,fn,xx.shape,self.szA.slices,self.szA.lsize,self.szA.gsize))
            with h5py.File(fn, 'w', driver='mpio', comm=self.scomm) as f:  #libver='latest'
                f.create_dataset('wn',data=wn)
                f.create_dataset('t', data=self.solver.sim_time)
                f.create_dataset('dt',data=dt)
                f.create_dataset('i', data=self.solver.iteration)
                d  = f.create_dataset(nm,self.szA.gsize, dtype=np.float) # , chunks=True) autochunking doesn't seem to help
                with d.collective: d[self.szA.slices]=xx
        self.wn[j]+=1
        if debug_flow: ju.logger.info('jrec.sf.save: Leaving, wn: {:4d}'.format(self.wn[j]))
 
class jit:
    def __init__(self,jrec,javrg,avrg,i0,i1,i2,dt,fnm,p):

        self.jdt=dt
        self.avrg=avrg
        self.fnm=fnm
        self.nm=javrg.nm
        self.var=dict()
        self.ii0=i0
        self.ii1=i1
        self.ii2=i2
        self.reduce=jrec.reduce
        self.jlsize=jrec.lsize
        self.slices=(slice(0,0,None),)
        if javrg.szA.dim==2:
            if i0 and i2:
                self.axis=(0,1)
                self.comm=javrg.szA.comm01
                self.d=jrec.Iw[0]*jrec.Iw[1]
                if javrg.szA.start[0]==0 and javrg.szA.start[1]==0: self.slices=tuple()
            elif i0:
                self.axis=(0)
                self.comm=javrg.szA.comm0
                self.d=jrec.Iw[0]
                if javrg.szA.start[0]==0: self.slices=(javrg.szA.slices[1],)
            elif i2:
                self.axis=(1)
                self.comm=javrg.szA.comm1
                self.d=jrec.Iw[1]
                if javrg.szA.start[1]==0: self.slices=(javrg.szA.slices[0],)
            else:
                self.axis=()
                self.comm=None
                self.d=None
                self.slices=(javrg.szA.slices[0],javrg.szA.slices[1])
        else:
            if i0 and i1 and i2:
                self.axis=(0,1,2)
                self.comm=javrg.szA.comm012
                if 'wyz' in jrec.parameters: self.d=jrec.Iw[0]*jrec.Iw[1]*jrec.Iw[2]*jrec.parameters['wyz']['g']
                else:                        self.d=jrec.Iw[0]*jrec.Iw[1]*jrec.Iw[2]
                if javrg.szA.start[0]==0 and javrg.szA.start[1]==0 and javrg.szA.start[2]==0: self.slices=tuple()
            elif i0 and i1:
                self.axis=(0,1)
                self.comm=javrg.szA.comm01
                self.d=jrec.Iw[0]*jrec.Iw[1]
                if javrg.szA.start[0]==0 and javrg.szA.start[1]==0: self.slices=(javrg.szA.slices[2],)
            elif i0 and i2:
                self.axis=(0,2)
                self.comm=javrg.szA.comm02
                self.d=jrec.Iw[0]*jrec.Iw[2]
                if javrg.szA.start[0]==0 and javrg.szA.start[2]==0: self.slices=(javrg.szA.slices[1],)
            elif i1 and i2:
                self.axis=(1,2)
                self.comm=javrg.szA.comm12
                if 'wyz' in jrec.parameters: self.d=jrec.Iw[1]*jrec.Iw[2]*jrec.parameters['wyz']['g']
                else:                        self.d=jrec.Iw[1]*jrec.Iw[2]
                if javrg.szA.start[1]==0 and javrg.szA.start[2]==0: self.slices=(javrg.szA.slices[0],)
            else:
                self.slices=(slice(0,0,None),slice(0,0,None))
                if i0:
                    self.axis=(0)
                    self.comm=javrg.szA.comm0
                    self.d=jrec.Iw[0]
                    if javrg.szA.start[0]==0 : self.slices=(javrg.szA.slices[1],javrg.szA.slices[2])
                elif i1:
                    self.axis=(1)
                    self.comm=javrg.szA.comm1
                    self.d=jrec.Iw[1]
                    if javrg.szA.start[1]==0 : self.slices=(javrg.szA.slices[0],javrg.szA.slices[2])
                elif i2:
                    self.axis=(2)
                    self.comm=javrg.szA.comm2
                    self.d=jrec.Iw[2]
                    if javrg.szA.start[2]==0 : self.slices=(javrg.szA.slices[0],javrg.szA.slices[1])
                else:
                    self.axis=()
                    self.comm=None
                    self.d=None
                    self.slices=(javrg.szA.slices[0],javrg.szA.slices[1],javrg.szA.slices[2])

        n=len(self.slices)
        self.size=np.ones((n,),dtype=np.int)
        for j in range(n): self.size[j]=self.slices[j].stop-self.slices[j].start
        if n==0: self.tsize=1
        else:    self.tsize=self.size.prod()

        self.color = self.tsize>0
        self.scomm=new_comm(ju.comm,self.color,ju.mpirank)
        f=list()
        if not i0: f.append(0)
        if not i1 and javrg.szA.dim==3: f.append(1)
        if not i2:  f.append(-1)
        self.gsize=javrg.szA.gsize[f]
        self.lsize=javrg.szA.lsize[f]
        #ju.logger.info('jrec.jit.init: javrg.gsize {}'.format(javrg.szA.gsize))
        #ju.logger.info('jrec.jit.init: javrg.lsize {}'.format(javrg.szA.lsize))
        #ju.logger.info('jrec.jit.init: self.gsize  {}'.format(self.gsize))
        #ju.logger.info('jrec.jit.init: self.lsize {}'.format(self.lsize))
        #quit(1)
        if self.fnm is not None:
            reset=False
            self.wn=0
            self.dir=self.fnm
            self.dir,self.wn=get_wn_fn(javrg.dir / self.fnm,self.fnm,reset)
            #print('jit: rank: {},  dir: {}, fnm: {}, dt: {}, wn: {}, nm:{}'.format(ju.mpirank,self.dir,self.fnm,self.jdt,self.wn,self.nm))
            #quit(1)
            self.wn=max(1,self.wn)
            self.wn  = ju.comm.bcast(self.wn,  root=0)    
            self.dir = ju.comm.bcast(self.dir,  root=0)    
            ju.logger.info('jit:  nm: {:5s} dt: {:6.3f} num: {:5d} nm: {}'.format(self.fnm,self.jdt,self.wn,self.nm))
            #if debug_slices: ju.logger.info('{} {} {} {}'.format(i0,i1,i2,avrg))
            if debug_slices: ju.logger.info('jit: slices: {}, size: {}, lsize: {}, self.gsize: {}'.format(self.slices,self.size,self.lsize,self.gsize))
 
    def print(self):
        ju.logger.info('jit: gsize  {}'.format(self.gsize))
        ju.logger.info('jit: lsize  {}'.format(self.lsize))
        ju.logger.info('jit: color  {}'.format(self.color))
        ju.logger.info('jit: slices {}'.format(self.slices))
        ju.logger.info('jit: axis   {}'.format(self.axis))
        ju.logger.info('jit: avrg   {}'.format(self.avrg))
        ju.logger.info('jit: jdt    {}'.format(self.jdt))
        ju.logger.info('jit: fnm    {}'.format(self.fnm))
        ju.logger.info('jit: nm     {}'.format(self.nm))
        ju.logger.info('jit: var    {}'.format(self.var))
        ju.logger.info('jit: d      {}'.format(self.d))

    def clear(self):
        if self.avrg:
            for a in self.nm:
                if not (a in ('ru','rb','rs','Au','Ab','As','Su','Sb','Ss') and len(self.lsize)>1): self.var[a]=np.zeros(self.lsize)
            #ju.logger.info('jit.clear self.var[{}].shape'.format(a,self.var[a].shape))
        else:
            for a in self.nm: self.var[a]=np.float(0)
        #ju.logger.info('jit.clear self.var.keys{}'.format(self.var.keys()))
        #quit()
        #ju.logger.info('jit.clear ' + self.fnm)
        #for a in self.nm: ju.logger.info('self.var[{:8s}] type {}'.format(a,type(self.var[a])));


        #for a in self.nm: ju.logger.info('clear: self.var[{:6s}].shape: {}'.format(a,self.var[a].shape))
        #ju.logger.info('clear: i0: {!r:5}, i1: {!r:5}, i2: {!r:5}, fnm: {:3s}, lsize: {}'.format(self.ii0,self.ii1,self.ii2,self.fnm,self.lsize))
        self.n  = 0
        self.t1 =  np.inf
        self.t2 = -np.inf
        self.i1 =  np.iinfo(np.int64).max
        self.i2 =  np.iinfo(np.int64).min
        self.dt = 0

    def int(self,x,p):
        if x is None: ju.logger.error(' jit.int x is None')
        if False: ju.logger.info('jrec.jit.int: Entering x.shape: {}, x.range: [{:6.3f} {:6.3f}], lsize: {}'.format(x.shape,ju.jmin(x),ju.jmax(x),self.jlsize))
        if not all(x.shape==self.jlsize): x = self.reduce(x,p) # Scale down if necessary
        if self.comm is None: return(x)
        x=self.comm.reduce(np.sum(x*self.d,axis=self.axis), root=0, op=ju.MPI.SUM)
        if x is not None:
            if debug_avrg: ju.logger.info('jrec.jit.int: Leaving  x.shape: {}'.format(x.shape))
            return(x)
        else:
            if debug_avrg: ju.logger.info('jrec.jit.int: Leaving size: {}, x.shape: None'.format(self.size))
            return(np.zeros(()))
              
    def intnr(self,x,sf=1):
        if x is None: ju.logger.error(' jit.int x is None')
        if False: ju.logger.info('jrec.jit.int: Entering x.shape: {}, x.range: [{:6.3f} {:6.3f}], lsize: {}'.format(x.shape,ju.jmin(x),ju.jmax(x),self.jlsize))
        if self.comm is None: return(x)
        x=self.comm.reduce(np.sum(x*self.d*sf,axis=self.axis), root=0, op=ju.MPI.SUM)
        if x is not None:
            if debug_avrg: ju.logger.info('jrec.jit.int: Leaving  x.shape: {}'.format(x.shape))
            return(x)
        else:
            if debug_avrg: ju.logger.info('jrec.jit.int: Leaving size: {}, x.shape: None'.format(self.size))
            return(np.zeros(()))
              
class javrg:
    def __init__(self,a,p,reset):
        self.szA=a.szJ
        self.dim    = a.szJ.dim
        self.solver = a.solver

        nm  = list()
        self.vr  = list()
        self.vsp = list()
        self.v1a = list()
        self.v1b = list()
        self.v2a = list()
        self.v2b = list()
        self.v2c = list()
        if ju.mpirank==0:
            pv=set(a.pv)
            self.v1a=list(set(
                ('d','p','s','b','u','v','w','SR','OM','SS','SB',#'rfu','rfv','rfw',
                 'dudx','dvdx','dwdx','dudy','dvdy','dwdy','dudz','dvdz','dwdz',
                 'ududx','ududy','ududz','udvdx','udvdy','udvdz','udwdx','udwdy','udwdz',
                 'vdudx','vdudy','vdudz','vdvdx','vdvdy','vdvdz','vdwdx','vdwdy','vdwdz',
                 'wdudx','wdudy','wdudz','wdvdx','wdvdy','wdvdz','wdwdx','wdwdy','wdwdz')
            ) & pv)
            for aa in self.v1a: self.v1b.append(aa)
            if 'b'  in pv: self.vsp.append('Eb')
            if 's'  in pv: self.vsp.append('Es')
            if 'b'  in pv: self.vr.append('rb')
            if 's'  in pv: self.vr.append('rs')
            if 'u'  in pv: self.vr.append('ru')
            if 'b'  in pv: self.vr.append('Ab')
            if 's'  in pv: self.vr.append('As')
            if 'u'  in pv: self.vr.append('Au')
            if 'b'  in pv: self.vr.append('Sb')
            if 's'  in pv: self.vr.append('Ss')
            if 'u'  in pv: self.vr.append('Su')
            
            v2a=list(('p','s','b','u','v','w','b','b','b','b','s','s','s','u','u','v','b','s','u','v','w'))
            v2b=list(('p','s','b','u','v','w','s','u','v','w','u','v','w','v','w','w','d','d','d','d','d'))
            nm = copy.copy(self.v1b)

            for j in self.vsp: nm.append(j)
            for j in self.vr:  nm.append(j)

            for j in range(len(v2a)):
                if v2a[j] in pv and v2b[j] in pv:
                    self.v2a.append(v2a[j])
                    self.v2b.append(v2b[j])
                    self.v2c.append(v2a[j]+v2b[j])
                    nm.append(v2a[j]+v2b[j])
            nm.sort()
        self.nm  = ju.comm.bcast(nm,  root=0)
        #ju.logger.info('dir(comm) {}'.format(dir(ju.comm)))
        #ju.logger.info('comm.info {}'.format(ju.comm.info))
        #print('rank: {} {}/{}, nm: {}'.format(ju.mpirank,ju.comm.rank,ju.comm.size,self.nm))
        #quit(1)
        self.vr  = ju.comm.bcast(self.vr,  root=0)
        self.vsp = ju.comm.bcast(self.vsp, root=0)
        self.v1a = ju.comm.bcast(self.v1a, root=0)
        self.v1b = ju.comm.bcast(self.v1b, root=0)
        self.v2a = ju.comm.bcast(self.v2a, root=0)
        self.v2b = ju.comm.bcast(self.v2b, root=0)
        self.v2c = ju.comm.bcast(self.v2c, root=0)

        self.dir = a.dir
        self.ji=list()

        jita(a,self, False, False, False, False, 'dtj',    'state', p)
        jita(a,self, False,  True, False, False, 'dtjx',   'x',     p)
        jita(a,self, False, False,  True, False, 'dtjy',   'y',     p)
        jita(a,self, False, False, False,  True, 'dtjz',   'z',     p)
        jita(a,self, False,  True,  True, False, 'dtjxy',  'xy',    p)
        jita(a,self, False, False,  True,  True, 'dtjyz',  'yz',    p)
        jita(a,self, False,  True, False,  True, 'dtjxz',  'xz',    p)
        jita(a,self, False,  True,  True,  True, 'dtjxyz', 'xyz',   p)
        jita(a,self,  True, False, False, False, 'dtja',   'a',     p)
        jita(a,self,  True,  True, False, False, 'dtjax',  'ax',    p)
        jita(a,self,  True, False,  True, False, 'dtjay',  'ay',    p)
        jita(a,self,  True, False, False,  True, 'dtjaz',  'az',    p)
        jita(a,self,  True,  True,  True, False, 'dtjaxy', 'axy',   p)
        jita(a,self,  True, False,  True,  True, 'dtjayz', 'ayz',   p)
        jita(a,self,  True,  True, False,  True, 'dtjaxz', 'axz',   p)
        jita(a,self,  True,  True,  True,  True, 'dtjaxyz','axyz',  p)
    
        self.n=len(self.ji)
        if self.n==0: return
        self.avrg = False
        for a in self.ji: self.avrg = self.avrg | a.avrg
        for a in self.ji: a.clear()
        
    def print(self):
        print('javrg rank:{}/{}'.format(ju.mpirank,ju.mpisize))
        print('dir(javrg) {}'.format(dir(self)))
        print('vars(javrg) {}'.format(vars(self)))
        #print(' gsize   {}'.format(self.gsize))
        #print(' lsize   {}'.format(self.lsize))
        #print(' xslice  {}'.format(self.xslice))
        #print(' zslice  {}'.format(self.zslice))
        #print(' nx      {}'.format(self.nx))
        #print(' nz      {}'.format(self.nz))
        ju.Barrier()
        
    def getvar(self):
        if self.n==0: return(set())
        if self.avrg: return(set(self.v1a))
        for a in self.ji:
            if  self.solver.sim_time >= a.wn*a.jdt: return(set(self.v1a))
        return(set())
          
    def update(self, dt, v):
        if debug_flow: ju.logger.info('javrg.update: Entering')
        for a in self.ji: 
            if a.avrg:
                #print('a.var.keys {}'.format(a.var.keys()))
                #print('v.keys {}'.format(v.keys()))
                #print('v1a {}'.format(self.v1a))
                #print('v1b {}'.format(self.v1b))
                #ju.logger.info('self.v1a {}'.format(self.v1a))
                #ju.logger.info('self.vba {}'.format(self.v1b))
                #ju.logger.info('var.keys {}'.format(a.var.keys()))
                #ju.logger.info('v.keys   {}'.format(v.keys()))
                #for j in range(len(self.v1a)):
                #    if v[self.v1a[j]] is None: ju.logger.info('v[{}] is None v1a[{}]'.format(self.v1a[j],j))
                #    if v[self.v1b[j]] is None: ju.logger.info('v[{}] is None v1b[{}]'.format(self.v1b[j],j))
                #    x1=v[self.v1a[j]];
                #    x2=v[self.v1a[j]];
                #    ju.logger.info('{} {} {} {} {} {} {}'.format(j,self.v1a[j],self.v1b[j],type(x1),type(x2),x1.shape,x2.shape))
                #    x3=a.int(x1)
                for j in range(len(self.v1a)): a.var[self.v1b[j]] += dt*a.int(v[self.v1a[j]],self.v1a[j])
                for j in range(len(self.v2a)): a.var[self.v2c[j]] += dt*a.int(v[self.v2a[j]]*v[self.v2b[j]],self.v2a[j]+self.v2b[j])
                if 'Es' in self.vsp:           a.var['Es']        += dt*a.int(ju.mix_entropy(v['s']),'s')
                if 'Eb' in self.vsp:           a.var['Eb']        += dt*a.int(ju.mix_entropy(v['b']),'b')
                if a.var['u'].ndim<2:
                    if 'ru' in self.vr:            a.var['ru']    += dt*             a.var['u']    /np.sqrt(np.maximum(ju.realeps,a.var['uu']))
                    if 'Au' in self.vr:            a.var['Au']    += dt*             a.var['u'] **2/        np.maximum(ju.realeps,a.var['uu'])
                    if 'Su' in self.vr:            a.var['Su']    += dt*np.maximum(0,a.var['uu'])  /                              a.var['u']
                    if 'rb' in self.vr:            a.var['rb']    += dt*np.maximum(0,a.var['b'])   /np.sqrt(np.maximum(ju.realeps,a.var['bb']))
                    if 'Ab' in self.vr:            a.var['Ab']    += dt*np.maximum(0,a.var['b'])**2/        np.maximum(ju.realeps,a.var['bb'])
                    if 'Sb' in self.vr:            a.var['Sb']    += dt*np.maximum(0,a.var['bb'])  /        np.maximum(ju.realeps,a.var['b'])
                    if 'rs' in self.vr:            a.var['rs']    += dt*np.maximum(0,a.var['s'])   /np.sqrt(np.maximum(ju.realeps,a.var['ss']))
                    if 'As' in self.vr:            a.var['As']    += dt*np.maximum(0,a.var['s'])**2/        np.maximum(ju.realeps,a.var['ss'])
                    if 'Ss' in self.vr:            a.var['Ss']    += dt*np.maximum(0,a.var['ss'])  /        np.maximum(ju.realeps,a.var['s'])
   
                a.n  += 1
                a.dt += dt
                a.t1=min(a.t1,self.solver.sim_time)
                a.i1=min(a.i1,self.solver.iteration)
                a.t2=max(a.t1,self.solver.sim_time)
                a.i2=max(a.i1,self.solver.iteration)
            elif self.solver.sim_time >= a.wn*a.jdt:
                a.n  = 1
                a.dt = dt
                a.t1 = self.solver.sim_time
                a.i1 = self.solver.iteration
                a.t2 = self.solver.sim_time
                a.i2 = self.solver.iteration
                for j in range(len(self.v1a)): a.var[self.v1b[j]] = a.int(v[self.v1a[j]],self.v1a[j])
                for j in range(len(self.v2a)): a.var[self.v2c[j]] = a.int(v[self.v2a[j]]*v[self.v2b[j]],self.v2a[j]+self.v2b[j])
                if 'Es' in self.vsp:           a.var['Es']        = a.int(ju.mix_entropy(v['s']),'s')
                if 'Eb' in self.vsp:           a.var['Eb']        = a.int(ju.mix_entropy(v['b']),'b')
                if a.var['u'].ndim<2:
                    if 'ru' in self.vr:            a.var['ru']    =              a.var['u']    /np.sqrt(np.maximum(ju.realeps,a.var['uu']))
                    if 'Au' in self.vr:            a.var['Au']    =              a.var['u'] **2/        np.maximum(ju.realeps,a.var['uu'])
                    if 'Su' in self.vr:            a.var['Su']    = np.maximum(0,a.var['uu'])  /                              a.var['u']
                    if 'rb' in self.vr:            a.var['rb']    = np.maximum(0,a.var['b'])   /np.sqrt(np.maximum(ju.realeps,a.var['bb']))
                    if 'Ab' in self.vr:            a.var['Ab']    = np.maximum(0,a.var['b'])**2/        np.maximum(ju.realeps,a.var['bb'])
                    if 'Sb' in self.vr:            a.var['Sb']    = np.maximum(0,a.var['bb'])  /        np.maximum(ju.realeps,a.var['b'])
                    if 'rs' in self.vr:            a.var['rs']    = np.maximum(0,a.var['s'])   /np.sqrt(np.maximum(ju.realeps,a.var['ss']))
                    if 'As' in self.vr:            a.var['As']    = np.maximum(0,a.var['s'])**2/        np.maximum(ju.realeps,a.var['ss'])
                    if 'Ss' in self.vr:            a.var['Ss']    = np.maximum(0,a.var['ss'])  /        np.maximum(ju.realeps,a.var['s'])
  
                #ju.logger.error('self.solver.sim_time {},a.wn*a.jdt  {}, a.wn {}'.format(self.solver.sim_time,a.wn*a.jdt,a.wn))
            if self.solver.sim_time >= a.wn*a.jdt and a.n>0:
                self.save(a)
                a.clear()
        if debug_flow: ju.logger.info('javrg.update: Exiting')
        
    def save(self,a):
        if a.n==0: return
        if a.color:
            #if ju.mpirank==0:
            #    if not a.dir.is_dir(): a.dir.mkdir()
            #a.scomm.barrier()
            fn = str(a.dir) + '/{}-{:05d}.hdf5'.format(a.fnm,a.wn)
            if debug_avrg: ju.logger.info('javrg.save {} {:8.4f}'.format(fn,self.solver.sim_time))
            # ju.logger.info('javrg.save {}'.format(fn))
            with h5py.File(fn, 'w', driver='mpio', comm=a.scomm) as f:
                #self.szA.h5c(f)
                f.create_dataset('wn',data=a.wn)
                f.create_dataset('t1',data=a.t1)
                f.create_dataset('t2',data=a.t2)
                f.create_dataset('dt',data=a.dt)
                f.create_dataset('i1',data=a.i1)
                f.create_dataset('i2',data=a.i2)
                f.create_dataset('n', data=a.n)
                for b in a.var:
                    #ju.logger.info('  {}: size {} slices {}'.format(b,a.gsize,a.slices))
                    d = f.create_dataset(b, a.gsize, dtype=np.float)
                    with d.collective: d[a.slices]=a.var[b]
        a.wn+=1
        
class jsz:
    def __init__(self,domain,AA):       
        self.AA     = AA
        self.dim    = domain.dim
        self.gsize  = domain.dist.grid_layout.global_shape(scales=AA)
        self.lsize  = domain.dist.grid_layout.local_shape(scales=AA)
        self.slices = domain.dist.grid_layout.slices(scales=AA)
        self.start  = domain.dist.grid_layout.start(scales=AA)
        if debug_flow: ju.logger.info('jsz.init: Calling mesh_info')
        self.mi     = ju.mesh_info(domain,AA)

        # Setup communicators for integrating in different dimensions and combinations
        if debug_flow: ju.logger.info('jsz.init: Preparing communicators')

        self.color0   = -1  
        self.color1   = -1  
        self.color2   = -1  
        self.color01  = -1 
        self.color02  = -1 
        self.color12  = -1 
        self.color012 = -1
        
        if self.dim==3:
            for k in range(ju.mpisize):
                if self.mi[0,k]==0 and self.mi[1,ju.mpirank]==self.mi[1,k] and self.mi[2,ju.mpirank]==self.mi[2,k]: self.color0=k
                if self.mi[1,k]==0 and self.mi[0,ju.mpirank]==self.mi[0,k] and self.mi[2,ju.mpirank]==self.mi[2,k]: self.color1=k
                if self.mi[2,k]==0 and self.mi[0,ju.mpirank]==self.mi[0,k] and self.mi[1,ju.mpirank]==self.mi[1,k]: self.color2=k
                if self.mi[0,k]==0 and self.mi[1,k]==0                     and self.mi[2,ju.mpirank]==self.mi[2,k]: self.color01=k
                if self.mi[0,k]==0 and self.mi[2,k]==0                     and self.mi[1,ju.mpirank]==self.mi[1,k]: self.color02=k
                if self.mi[1,k]==0 and self.mi[2,k]==0                     and self.mi[0,ju.mpirank]==self.mi[0,k]: self.color12=k
                if self.mi[0,k]==0 and self.mi[1,k]==0                     and self.mi[2,k]==0:                     self.color012=k
        else:
            for k in range(ju.mpisize):
               if self.mi[0,k]==0 and self.mi[1,ju.mpirank]==self.mi[1,k]: self.color0=k
               if self.mi[1,k]==0 and self.mi[0,ju.mpirank]==self.mi[0,k]: self.color1=k
               if self.mi[0,k]==0 and self.mi[1,k]==0:                     self.color01=k
            self.color2   = k
            self.color02  = self.color0
            self.color12  = self.color1
            self.color012 = self.color01
               
        self.key0   = ju.mpirank
        self.key1   = ju.mpirank
        self.key2   = ju.mpirank
        self.key01  = ju.mpirank
        self.key02  = ju.mpirank
        self.key12  = ju.mpirank
        self.key012 = ju.mpirank
            
        if self.color0   == ju.mpirank: self.key0=0
        if self.color1   == ju.mpirank: self.key1=0
        if self.color2   == ju.mpirank: self.key2=0
        if self.color01  == ju.mpirank: self.key01=0
        if self.color02  == ju.mpirank: self.key02=0
        if self.color12  == ju.mpirank: self.key12=0
        if self.color012 == ju.mpirank: self.key012=0

        #ju.logger.error('{:3d} color1: {:3d} key1: {:3d}'.format(ju.mpirank, self.color1, self.key1))
        if debug_flow: ju.logger.info('jsz.init: Creating communicators')
        self.comm0   = new_comm(ju.comm,self.color0,   self.key0)
        self.comm1   = new_comm(ju.comm,self.color1,   self.key1)
        self.comm2   = new_comm(ju.comm,self.color2,   self.key2)
        self.comm01  = new_comm(ju.comm,self.color01,  self.key01)
        self.comm02  = new_comm(ju.comm,self.color02,  self.key02)
        self.comm12  = new_comm(ju.comm,self.color12,  self.key12)
        self.comm012 = new_comm(ju.comm,self.color012, self.key012)
        self.rank0   = self.comm0.Get_rank()
        self.rank1   = self.comm1.Get_rank()
        self.rank2   = self.comm2.Get_rank()
        self.rank01  = self.comm01.Get_rank()
        self.rank02  = self.comm02.Get_rank()
        self.rank12  = self.comm12.Get_rank()
        self.rank012 = self.comm012.Get_rank()

        if debug_flow: ju.logger.info('jsz.init: Creating slices')
         
        self.slicesc=np.ndarray((self.dim), dtype=type(slice))
        self.slicesc[0] = (slice(0,0,None),)
        self.slicesc[1] = (slice(0,0,None),)
        if self.dim==3:
           self.slicesc[2] = (slice(0,0,None),)
           if self.start[0]==0 and self.start[1]==0: self.slicesc[2] = (self.slices[2],)
           if self.start[0]==0 and self.start[2]==0: self.slicesc[1] = (self.slices[1],)
           if self.start[1]==0 and self.start[2]==0: self.slicesc[0] = (self.slices[0],)
        else:
            if self.start[0]==0:                     self.slicesc[1] = (self.slices[1],)
            if self.start[1]==0:                     self.slicesc[0] = (self.slices[0],)
            
        
        if debug_flow: ju.logger.info('jsz.init: Creating sizes')
    
        #print('vars(domain) {}'.format(vars(domain)))
        #print('domain.bases.name {}'.format(domain.bases[:].name))
        #print('dir(domain) {}'.format(dir(domain)))
        #print('dir(domain.bases) {}'.format(dir(domain.bases)))
        #print('vars(domain.grid) {}'.format(vars(domain.grid)))
        #print('vars(domain.grids) {}'.format(vars(domain.grids(1))))
        #print('domain.grid.axis {}'.format(domain.grid.axis))
        #quit()
        if debug_flow: ju.logger.info('jsz.init: Preparing coordinates')
        self.cn=list()
        for j in range(self.dim): self.cn.append(domain.bases[j].name)
        self.z=np.array([],dtype=np.float)
        if self.dim==3:
            self.c=list((np.zeros((0)),np.zeros((0)),np.zeros((0))))
            if self.start[1]==0 and self.start[2]==0: self.c[0]=np.reshape(copy.copy(domain.grids(AA)[0][:,0,0]),(self.lsize[0]))
            if self.start[0]==0 and self.start[2]==0: self.c[1]=np.reshape(copy.copy(domain.grids(AA)[1][0,:,0]),(self.lsize[1]))
            if self.start[0]==0 and self.start[1]==0: self.c[2]=np.reshape(copy.copy(domain.grids(AA)[2][0,0,:]),(self.lsize[2]))
        else:
            self.c=list((np.zeros((0)),np.zeros((0))))
            #print('domain.grids(AA)[0] {}'.format(domain.grids(AA)[0]))
            #print('self.start {}'.format(domain.grids(AA)[0]))
            #print('self.c   ] {}'.format(domain.grids(AA)[0]))
            if self.start[1]==0: self.c[0]=np.reshape(copy.copy(domain.grids(AA)[0][:,0]),(self.lsize[0]))
            if self.start[0]==0: self.c[1]=np.reshape(copy.copy(domain.grids(AA)[1][0,:]),(self.lsize[1]))

    def h5c(self,f):
        for k in range(self.dim):
            x   = f.create_dataset(self.cn[k], (self.gsize[k],), dtype=np.float)
            x[self.slicesc[k]]=self.c[k]
  
    def print(self):
        for j in range(ju.mpisize):
            ju.Barrier()
            if j==ju.mpirank:
                ju.logger.error('jsz: {:3d}, AA: {}, gsize: {}, lsize: {}, start: {}'.format(ju.mpirank,self.AA,self.gsize,self.lsize,self.start))
                ju.logger.error('        slices: {}'.format(self.slices))
                #ju.logger.error('         slices1.shape: {} slices2.shape: {}'.format(self.slices1.shape,self.slices2.shape))
                #ju.logger.error('         size1: {} size2: {} size3 {}:'.format(self.size1,self.size2,self.size3))
                #for j in range(self.dim): ju.logger.error('  slices1[{:1d},0]={}, slices1[{1:d},1]={}'.format(j,self.slices1[j,0],j,self.slices1[j,1]))
                #for j in range(self.dim): ju.logger.error('  slices1[{:1d},0]={}, slices1[{:1d},1]={}'.format(j,self.slices1[j,0],j,self.slices1[j,1]))
                #for j in range(len(self.slices2)): ju.logger.error('jsz: {:3d} slices2[{:1d}]'.format(ju.mpirank,self.slices2[j]))
                #ju.logger.error('        slices1: {}'.format(self.slices1))
                #ju.logger.error('        slices2: {}'.format(self.slices2))
                #ju.logger.error('        cn: {}'.format(self.cn))
                #ju.logger.error('         c: {}'.format(self.c))
                #ju.logger.error('     {}  {} {}'.format(j,len(self.cn),len(self.c)))
                for j in range(self.dim):
                    if len(self.c[j])>0: ju.logger.error('        {:1d}  {:1s} {:4d} [{:7.3f} {:7.3f}]'.format(j,self.cn[j],len(self.c[j]),self.c[j][0],self.c[j][-1]))
                    else:                ju.logger.error('        {:1d}  {:1s} {:4d}'.format(j,self.cn[j],len(self.c[j])))
               
class jrec:
    fld = None
    def __init__(self,domain,ddir,p,reset,dfd):
        if debug_flow: ju.logger.info('jrec.init: Init started')
        self.dfd = dfd # dynamic divergence
            
        self.Chebz  = p['Tz']=='Cheb'
        self.AA  = p['AA']  # Antaliasing level used in computations
        self.AAS = p['AAS'] # Antialising level for saving program state
        self.AAJ = p['AAJ'] # Antialising level for saving output files
        # Probably should add an extra level for calculating quadratic terms
        self.sz1 = jsz(domain,1)
        if   self.AA ==1:        self.szA = self.sz1
        else:                    self.szA = jsz(domain,self.AA)
        if   self.AAS==1:        self.szS = self.sz1
        elif self.AAS==self.AA:  self.szS = self.szA
        else:                    self.szS = jsz(domain,self.AAS)
        if   self.AAJ==1:        self.szJ = self.sz1
        elif self.AAJ==self.AA:  self.szJ = self.szA
        elif self.AAJ==self.AAS: self.szJ = self.szS
        else:                    self.szJ = jsz(domain,self.AAJ)
        self.domain = domain
        self.Iw=ju.int_weights(self.AAJ,domain)[0]
        ju.check_error()

        
        self.dtjadv  = dp.get_bool_param(p,'dtjadv')  # Save advective terms
        self.dtjsw   = dp.get_bool_param(p,'dtjsw')   # Save strain rate and rotation rate
        self.dtjd1   = dp.get_bool_param(p,'dtjd1')   # Save first order derivatives
        self.interval=np.ndarray((self.szA.dim,2),dtype=np.float)
        self.btype=np.ndarray((self.szA.dim,),dtype=np.dtype('S9'))   # Can't write variable length strings with parallel hdf5
        self.dim=self.szA.dim
        for j in range(self.szA.dim):
            self.interval[j,:]=domain.bases[j].interval
            self.btype[j]=ju.get_bases_type(domain.bases[j])
                        
        self.PC=dict()
        self.PC['x']=np.ones((self.dim,))
        self.PC['y']=np.ones((self.dim,))
        self.PC['z']=np.ones((self.dim,))
        if 'Tx' in p:
            self.Tx=p['Tx']
            if self.Tx=='SinCos': self.PC['x'][0]=-1
        if 'Ty' in p:
            self.Ty=p['Ty']
            if self.Ty=='SinCos': self.PC['y'][1]=-1
        if 'Tz' in p:
            self.Tz=p['Tz']
            if self.Tz=='SinCos': self.PC['z'][-1]=-1

        self.dir    = ddir
        self.var    = dict()
        self.fno    = Path('')
        jrec.fld=domain.new_field()  
        self.fld    = jrec.fld
        ju.logger.info('Initialising jrec: AA={:3.1f}'.format(self.AA))

        if debug_flow: ju.logger.info('jrec.init: Finding variables')
        
        if debug_flow: ju.logger.info('jrec.init: Init complete')
    
        self.gl=domain.dist.grid_layout
        self.grids=domain.grids
        self.gsize=self.gl.global_shape(scales=self.AAS)
        self.lsize=self.gl.local_shape(scales=self.AAS)
        self.slices=self.gl.slices(scales=self.AAS)
        self.color=np.all(self.lsize>0)
        self.scomm=new_comm(ju.comm,self.color, ju.mpirank)
        xx,self.write_number=get_wn_fn(self.dir/'final','state',reset)
        ju.logger.info('jrec.init: self.write_number {}'.format(self.write_number))

        with h5py.File(str(self.dir / 'coord.hdf5'), 'w', driver='mpio', comm=ju.comm) as f:
            for k in range(self.sz1.dim):
                x   = f.create_dataset(self.sz1.cn[k], (self.sz1.gsize[k],), dtype=np.float)
                x[self.sz1.slicesc[k]]=self.sz1.c[k]
            for k in range(self.szA.dim):
                x   = f.create_dataset('A'+self.szA.cn[k], (self.szA.gsize[k],), dtype=np.float)
                x[self.szA.slicesc[k]]=self.szA.c[k]
            for k in range(self.szJ.dim):
                x   = f.create_dataset('J'+self.szJ.cn[k], (self.szJ.gsize[k],), dtype=np.float)
                x[self.szJ.slicesc[k]]=self.szJ.c[k]
            for k in range(self.szS.dim):
                x   = f.create_dataset('S'+self.szS.cn[k], (self.szS.gsize[k],), dtype=np.float)
                x[self.szS.slicesc[k]]=self.szS.c[k]
            f.create_dataset('AA', data=self.AA)
            f.create_dataset('AAS',data=self.AAS)
            f.create_dataset('AAJ',data=self.AAJ)
        self.fcomm=None
        self.oy=None
        
    def set_solver(self,solver,problem,p,reset,state_ev,avar):
        self.parameters = problem.parameters

        self.variables  = problem.variables
        self.solver     = solver
        epv=list()
        if self.dtjadv:
            if self.dim==2: epv += list(('ududx','wdudz','udwdx','wdwdz'))           
            else:           epv += list(('ududx','vdudy','wdudz','udvdx','vdvdy','wdvdz','udwdx','vdwdy','wdwdz'))
        if self.dtjd1:
            if self.dim==2: epv += list(('dudx',       'dwdx',                     'dudz',       'dwdz'))
            else:           epv += list(('dudx','dvdx','dwdx','dudy','dvdy','dwdy','dudz','dvdz','dwdz'))
        if self.dtjsw:
            epv += list(('SR','OM'))
            if 's' in self.variables: epv += list(('SS',))
            if 'b' in self.variables: epv += list(('SB',))
        self.pv = list(set(('rfu','rfv','rfw','p','b','s','u','v','w')) & (set(self.variables) | set(self.parameters)))+epv
        if self.dfd and 'fd' in problem.parameters: self.pv.append('d')
            
        ju.logger.info('pv:  {}'.format(self.pv))
        #ju.logger.info('epv: {}'.format(epv))
        # It is very important that all the names are in the correct order otherwise h5py deadlocks can occur
        # So we do construct things on the root node and then broadcast them
        if debug_flow: ju.logger.info('jrec.init: Calling javrg')
        self.javrg = javrg(self,p,reset)
        sf=list([])
        dt=list([])
        if ju.mpirank==0:
            for a in self.pv:
                dtnm='dtj{}'.format(a)
                if dtnm in p:
                    if p[dtnm] > 0:
                        sf.append(a)
                        dt.append(p[dtnm])
                
        sf=ju.comm.bcast(sf, root=0)        
        dt=ju.comm.bcast(dt, root=0)
        if debug_flow: ju.logger.info('jrec.ssolver: Calling single_field')
        self.sf = single_field(self,reset,sf,dt)

        if 'dtjsev'  in p: self.sev=sdict(self,reset,state_ev,'sev',p['dtjsev'])
        else:              self.sev=None
        if 'dtjavar' in p: self.avar=sdict(self,reset,avar,'avar',p['dtjavar'])
        else:              self.avar=None

        pv=list(set(('p','b','s','u','v','w')) & (set(self.variables)))
        self.parity=dict()
        self.parity['d']=np.ones((self.dim,))
        self.parity['dd']=np.ones((self.dim,))
        if self.dim==3:  cc=('x','y','z')
        else:            cc=('x','z')
        for p1 in pv:
            self.parity[p1]=self.get_parity(solver.state[p1])
            self.parity[p1+'d']=self.get_parity(solver.state[p1])
            for c in cc: self.parity['d'+p1+'d'+c]=self.get_parity(solver.state[p1])*self.PC[c]
            for p2 in pv:
                self.parity[p1+p2]=self.get_parity(solver.state[p1])*self.get_parity(solver.state[p2])
                for c in cc: self.parity[p1+'d'+p2+'d'+c]=self.get_parity(solver.state[p1])*self.get_parity(solver.state[p2])*self.PC[c]
        ju.write_dict('parity.hdf5',self.parity)
        
    def reduce(self,a,p):
        if self.AA==self.AAJ or ju.jtype(a) != 'ndarray' : return(a)
        if np.all(a.shape==self.lsize): return
        self.fld.set_scales(self.AA)
        if type(p) is str: p=self.parity[p] 
        self.fld.meta['x']['parity']=p[0]
        if self.dim==3: self.fld.meta['y']['parity']=p[1]
        self.fld.meta['z']['parity']=p[-1]
        self.fld['g']=a
        self.fld.set_scales(self.AAJ)
        return(np.copy(self.fld['g'])) # This copy is critical

    def square(self,a):
        if self.AA==self.AAJ: return(a**2)
        if np.all(a.shape==self.lsize): return
        self.fld.set_scales(self.AA)
        self.fld.meta['x']['parity']=1
        if self.dim==3: self.fld.meta['y']['parity']=1
        self.fld.meta['x']['parity']=1
        self.fld['g']=a**2
        self.fld.set_scales(self.AAJ)
        return(np.copy(self.fld['g']))

    # Do all calculations with AA scale on the fields and then use the above to go back to AAJ
    # This should reduce aliasing
    def update(self, dt):
        if ju.mpirank==0: pv = self.javrg.getvar() | self.sf.getvar()
        else: pv =set()
        pv=list(pv)
        
        pv=ju.comm.bcast(pv, root=0)
        
        if len(pv)==0: return
        
        if debug_var or debug_flow: ju.logger.info('jrec.update: t:{:7.3f}, pv:{}'.format(self.solver.sim_time,pv))
        v=dict()
        for a in pv: # Take care to minimise the memory footprint
            if a in self.variables:
                self.solver.state[a].set_scales(self.AA) # Do not use keep_data=True This destroys performance
                self.solver.state[a]['g']         # Switch into grid space  
                v[a]=self.solver.state[a]['g']    # Get a reference to the grid space data but don't copy it
            elif a in self.parameters:
                ju.logger.info('{} in self.parameters'.format(a))
                self.parameters[a].set_scales(self.AA) # Do not use keep_data=True This destroys performance
                self.parameters[a]['g']         # Switch into grid space  
                v[a]=self.parameters[a]['g']    # Get a reference to the grid space data but don't copy it
        if 'd' in pv:
            self.parameters['fd'].set_scales(self.AA) 
            self.parameters['fd']['g']         
            v['d']=self.parameters['fd']['g']    
 
        if self.dtjadv or self.dtjsw or self.dtjd1:
            der=dict()
            
            if self.dim==3: aa=list(('u','v','w'))
            else:           aa=list(('u',    'w'))
            if self.dim==3: cc=('x','y','z')
            else:           cc=('x',    'z')
            for c in cc:    der[c]=dict()
            if self.dtjsw:
                if 's' in self.variables: aa += list(('s',)) 
                if 'b' in self.variables: aa += list(('b',)) 
            for a in aa:
                self.solver.state[a].set_scales(self.AA)
                aa=np.copy(self.solver.state[a]['g'])
                for c in cc:
                    if [a+c] in self.variables:
                        self.solver.state[a+c].set_scales(self.AA)
                        der[c][a]=np.copy(self.solver.state[a+c]['g'])
                    else:
                        self.solver.state[a].differentiate(c,out=self.fld)
                        self.fld['g']
                        der[c][a]=np.copy(self.fld['g'])
                self.solver.state[a]['g']=aa

        if self.dtjadv:
            #ju.logger.info('pv {}'.format(pv.keys()))
            #ju.logger.info('v.keys {}'.format(v.keys()))
            #ju.logger.info('der.keys {}'.format(der.keys()))
            #ju.logger.info('der[x].keys {}'.format(der['x'].keys()))
            #ju.logger.info('der[z].keys {}'.format(der['z'].keys()))
            if self.dim==2:
                v['ududx']=self.reduce(v['u']*der['x']['u'],self.parity['ududx'])
                v['udwdx']=self.reduce(v['u']*der['x']['w'],self.parity['udwdx'])
                v['wdudz']=self.reduce(v['w']*der['z']['u'],self.parity['wdudz'])
                v['wdwdz']=self.reduce(v['w']*der['z']['w'],self.parity['wdwdz'])
            else:                    
                v['ududx']=self.reduce(v['u']*der['x']['u'],self.parity['ududx'])
                v['udvdx']=self.reduce(v['u']*der['x']['v'],self.parity['udvdx'])
                v['udwdx']=self.reduce(v['u']*der['x']['w'],self.parity['udwdx'])
                v['vdudy']=self.reduce(v['v']*der['y']['u'],self.parity['vdudy'])
                v['vdvdy']=self.reduce(v['v']*der['y']['v'],self.parity['vdvdy'])
                v['vdwdy']=self.reduce(v['v']*der['y']['w'],self.parity['vdwdy'])
                v['wdudz']=self.reduce(v['w']*der['z']['u'],self.parity['wdudz'])
                v['wdvdz']=self.reduce(v['w']*der['z']['v'],self.parity['wdvdz'])
                v['wdwdz']=self.reduce(v['w']*der['z']['w'],self.parity['wdwdz'])
        if self.dtjsw:
            if self.dim==3:   #Strain rate tensor squared
                v['SR']   = self.square(der['x']['u'])+self.square(der['y']['v'])+self.square(der['z']['w']) \
                    + ( self.square(der['x']['v']+der['y']['u']) \
                    +   self.square(der['x']['w']+der['z']['u']) \
                    +   self.square(der['y']['w']+der['z']['v']))/2
                v['OM'] =  (self.square(der['z']['v']-der['y']['w'])+self.square(der['x']['w']-der['z']['u'])+self.square(der['y']['u']-der['x']['v']))/2
                if 's' in self.variables: v['SS'] =  self.square(der['x']['s'])+self.square(der['y']['s'])+self.square(der['z']['s'])
                if 'b' in self.variables: v['SB'] =  self.square(der['x']['b'])+self.square(der['y']['b'])+self.square(der['z']['b'])
                
            else:
                v['SR']   = self.square(der['x']['u'])+self.square(der['z']['w'])+self.square(der['x']['w']+der['z']['u'])/2   
                v['OM']   =                                                       self.square(der['x']['w']-der['z']['u'])/2
                if 's' in self.variables: v['SS'] =  self.square(der['x']['s'])+self.square(der['z']['s'])
                if 'b' in self.variables: v['SB'] =  self.square(der['x']['b'])+self.square(der['z']['b'])
                
        if self.dtjd1:
            if self.dim==2:
                v['dudx']=self.reduce(der['x']['u'],self.parity['dudx'])
                v['dwdx']=self.reduce(der['x']['w'],self.parity['dwdx'])
                v['dudz']=self.reduce(der['z']['u'],self.parity['dudz'])
                v['dwdz']=self.reduce(der['z']['w'],self.parity['dwdz'])
            else:
                v['dudx']=self.reduce(der['x']['u'],self.parity['dudx'])
                v['dvdx']=self.reduce(der['x']['v'],self.parity['dvdx'])
                v['dwdx']=self.reduce(der['x']['w'],self.parity['dwdx'])
                v['dudy']=self.reduce(der['y']['u'],self.parity['dudy'])
                v['dvdy']=self.reduce(der['y']['v'],self.parity['dvdy'])
                v['dwdy']=self.reduce(der['y']['w'],self.parity['dwdy'])
                v['dudz']=self.reduce(der['z']['u'],self.parity['dudz'])
                v['dvdz']=self.reduce(der['z']['v'],self.parity['dvdz'])
                v['dwdz']=self.reduce(der['z']['w'],self.parity['dwdz'])
 
        #ju.logger.info('pv {}'.format(pv))
        #ju.logger.info('v.keys {}'.format(v.keys()))
        
             
        if debug_flow: ju.logger.info('jrec.update: calling jrec.javrg.update')
        self.javrg.update(dt,v)
        
        if debug_flow: ju.logger.info('jrec.update: calling jrec.sf.update')
        self.sf.update(dt,v)

        if self.sev  is not None: self.sev.update(dt)
        if self.avar is not None: self.avar.update(dt)
        if debug_flow: ju.logger.info('jrec.update: Complete')
        #ju.logger.error('{:4d} {} {}'.format(ju.mpirank,a,self.var[a].shape))
        #self.savexyz(ju.basedir/'Jdudx6.hdf5','dudx',v['dudx'],self.AAJ)
        #self.savexyz(ju.basedir/   'Ju2.hdf5',   'u',self.reduce(v['u'],self.parity['u']),   self.AAJ)
        #quit()
        
    def save(self):
        for a in self.javrg.ji: 
            if a.avrg:
                self.javrg.save(a)
                a.clear()


    def saven(self,fn,nm,data,AA,d):
        start  = self.gl.start(scales=AA)
        dim=len(start)
        if   d=='x': n=0
        elif d=='y': n=1
        elif d=='z': n=dim-1
        gsize  = (self.gl.global_shape(scales=AA)[n],)
        lsize  = (self.gl.local_shape(scales=AA)[n],)
        if dim==3:
            if   n==0: gsl=(slice(0,lsize),0,0)
            elif n==1: gsl=(0,slice(0,lsize),0)
            elif n==2: gsl=(0,0,slice(0,lsize))
        else:
            if   n==0: gsl=(slice(0,lsize),0)
            elif n==1: gsl=(0,slice(0,lsize))
            
        start=np.delete(start,[n])
        slices = self.gl.slices(scales=AA)[n]
        color  = int(lsize[0]>0 and np.all(start==0))
        c=new_comm(ju.comm,color, ju.mpirank)
        #ju.logger.info('fn:{}, nm:{}, d:{}, n:{}, start:{} lsize:{} slices:{}'.format(fn,nm,d,n,start,lsize,slices))
        if color==1:
            with h5py.File(str(self.dir / fn), 'w', driver='mpio', comm=c) as f:
                #                x=f.create_dataset(self.szA.cn[n], gsize, dtype=np.float)
                #ju.logger.info('slices:{}, gsl:{}, lsize:{}, x.shape:{}'.format(slices,gsl,lsize,x.shape))
                #x[slices]=np.reshape(self.grids(AA)[n][gsl],lsize)
                if type(nm)==list or type(nm)==tuple:
                    for (a,b) in zip(nm, data):
                        d  = f.create_dataset(a,gsize, dtype=np.float)
                        with d.collective: d[slices]=np.reshape(b,lsize)
                else:
                    d  = f.create_dataset(nm,gsize, dtype=np.float)
                    with d.collective: d[slices]=np.reshape(data,lsize)

    def saved(self,fn,aa,d,v,AA,mode='w'):
        debug_saved=False
        if debug_saved:ju.logger.info('jrec.saved: fn: {}, nm: {}, d: {}, AA: {}'.format(fn,aa,d,AA))

        if d.find('y')>=0 and self.dim==2: return

        nm=list()
        if isinstance(v,tuple):
            vo=v
            v=dict()
            for a, vv in zip(aa, vo):
                nm.append(a)
                v[a]=vv
        elif ju.mpirank==0:
            for a in set(aa) & set([*v]):
                if v[a] is not None: nm.append(a)
        nm=ju.comm.bcast(nm, root=0)
        if debug_saved:ju.logger.info('jrec.saved: nm: {}'.format(nm))
 
        if len(nm)==0: return

        nd = len(d)
        mn=-np.ones((nd,),dtype=np.int64)
        for k in range(nd):
            if   d[k]=='x': mn[k]=0
            elif d[k]=='y': mn[k]=1
            elif d[k]=='z': mn[k]=self.dim-1
        if np.any(mn<0): 
            ju.logger.info('jrec.saved: fn:{} nm:{} dim:{} d:{} mn:{}'.format(fn,nm,self.dim,d,mn))
            return()
        k=np.array(list(set(range(self.dim))-set(mn)))
        if debug_saved:
            ju.logger.info('jrec.saved:  k:{}'.format(k))
            ju.logger.info('jrec.saved:  d:{}'.format(d))
            #ju.logger.info('jrec.saved:  k:{}  type(k):{}  type(k[0]):{} shape:{}'.format( k, type(k), type(k[0]), k.shape))
            ju.logger.info('jrec.saved:  k:{}'.format( k))
            ju.logger.info('jrec.saved: mn:{} type(mn):{} type(mn[0]):{} shape:{}'.format(mn,type(mn),type(mn[0]),mn.shape))
        gsize  = self.gl.global_shape(scales=AA)[mn]
        lsize  = self.gl.local_shape(scales=AA)[mn]
        if len(k)>0: start = self.gl.start(scales=AA)[k]
        else:        start = np.array([0])
        sslices = self.gl.slices(scales=AA)#mn
        slices=()
        for k in range(len(mn)): slices=slices+(sslices[mn[k]],)
        #ju.logger.info('slices {}'.format(slices))
        #ju.logger.info('type(slices) {}'.format(type(slices)))
        #quit(1)
        color  = np.all(start==0) and np.all(lsize>0)
        c=new_comm(ju.comm,np.int(color), ju.mpirank)

 
        p=str(self.check_dir(self.dir / fn))
        if color:
            with h5py.File(p, mode, driver='mpio', comm=c) as f:
                for a in nm:
                    if a in f: del f[a]
                    d  = f.create_dataset(a,gsize, dtype=np.float)
                    if ju.jtype(v[a])[0]=='f':
                        v[a].set_scales(AA)
                        with d.collective: d[slices]=v[a]['g'].reshape(lsize)
                    else:
                        with d.collective: d[slices]=v[a].reshape(lsize)

                    
    def savex(self,fn,nm,data,AA):
        gsize  = (self.gl.global_shape(scales=AA)[n],)
        lsize  = (self.gl.local_shape(scales=AA)[n],)
        start  = self.gl.start(scales=AA)
        dim=len(start)
        start=start[slice(1,dim)]
        slices = self.gl.slices(scales=AA)[n]
        color  = int(lsize[n]>0 and np.all(start==0))
        c=new_comm(ju.comm,color, ju.mpirank)
        if color==1:
            with h5py.File(str(self.dir / fn), 'w', driver='mpio', comm=c) as f:
                x=f.create_dataset(self.szA.cn[n], gsize, dtype=np.float)
                if dim==3: x[slices]=np.reshape(self.grids(AA)[n][:,0,0],lsize)
                else:      x[slices]=np.reshape(self.grids(AA)[n][:,0],lsize)
                if type(nm)==list or type(nm)==tuple:
                    for (a,b) in zip(nm, data):
                        d  = f.create_dataset(a,gsize, dtype=np.float)
                        with d.collective: d[slices]=np.reshape(b,lsize)
                else:
                    d  = f.create_dataset(nm,gsize, dtype=np.float)
                    with d.collective: d[slices]=np.reshape(data,lsize)

    
    def saveyz(self,fn,nm,data,AA):
        lsize  = self.gl.local_shape(scales=AA)
        dim    = len(lsize)
        if dim==2:
            self.saven(fn,nm,data,AA,'y')
            return
        gsize  = self.gl.global_shape(scales=AA)[[1,2]]
        slices = self.gl.slices(scales=AA)
        slices = (slices[1],slices[2])
        cn     = [self.szA.cn[1],self.szA.cn[2]]
        start  = self.gl.start(scales=AA)
        color  = int(start[0]==0 and lsize[1]>0 and lsize[2]>0)
        lsize  = [lsize[1],lsize[2]]
        c=new_comm(ju.comm,color, ju.mpirank)
        if color==1:
            p=self.check_dir(self.dir / fn)
            with h5py.File(str(p), 'w', driver='mpio', comm=c) as f:
                y=f.create_dataset(cn[0], (gsize[0],), dtype=np.float)
                z=f.create_dataset(cn[1], (gsize[1],), dtype=np.float)
                if  start[0]==0 and start[2]==0: y[slices[0]]=np.reshape(self.grids(AA)[1][0,:,0],(lsize[0]))
                if  start[0]==0 and start[1]==0: z[slices[1]]=np.reshape(self.grids(AA)[2][0,0,:],(lsize[1]))
                if type(nm)==list or type(nm)==tuple:
                    for (a,b) in zip(nm, data):
                        #ju.logger.info('gsize {}, lsize {}, slices {}, b.shape {}'.format(gsize,lsize,slices,b.shape))
                        d  = f.create_dataset(a,gsize, dtype=np.float)
                        with d.collective: d[slices]=b.reshape(lsize)
                else:
                    d  = f.create_dataset(nm,gsize, dtype=np.float)
                    with d.collective: d[slices]=data.reshape(lsize)

    def savexz(self,fn,nm,data,AA):
        lsize  = self.gl.local_shape(scales=AA)
        dim    = len(lsize)
        dx=0
        dz=dim-1
        gsize  = self.gl.global_shape(scales=AA)[[0,dz]]
        slices = self.gl.slices(scales=AA)
        slices = (slices[0],slices[dz])
        cn     = [self.szA.cn[0],self.szA.cn[dz]]
        
        start  = self.gl.start(scales=AA)
        
        color  = np.int(np.all(lsize[[0,-1]]>0))
        lsize  = [lsize[0],lsize[dz]]
        c=new_comm(ju.comm,color, ju.mpirank)
        if color==1:
            p=self.check_dir(self.dir / fn)
            f=h5py.File(str(p), 'w', driver='mpio', comm=c)
            x=f.create_dataset(cn[0], (gsize[0],), dtype=np.float)
            z=f.create_dataset(cn[1], (gsize[1],), dtype=np.float)
            if dim==3:
                if  start[1]==0 and start[2]==0: x[slices[0]]=np.reshape(self.grids(AA)[0][:,0,0],(lsize[0]))
                if  start[1]==0 and start[0]==0: z[slices[1]]=np.reshape(self.grids(AA)[2][0,0,:],(lsize[1]))
            else:
                if                  start[1]==0: x[slices[0]]=np.reshape(self.grids(AA)[0][:,0],(lsize[0]))
                if                  start[0]==0: z[slices[1]]=np.reshape(self.grids(AA)[1][0,:],(lsize[1]))
            if type(nm)==list or type(nm)==tuple:
                for (a,dd) in zip(nm, data):
                    d  = f.create_dataset(a,gsize, dtype=np.float)
                    with d.collective: d[slices]=dd.reshape(lsize)
            else:
                d  = f.create_dataset(nm,gsize, dtype=np.float)
                with d.collective: d[slices]=data.reshape(lsize)
            f.close()
                        
    def savez(self,fn,nm,data,AA):
        lsize  = self.gl.local_shape(scales=AA)
        dim    = len(lsize)
        dx=0
        dz=dim-1
        gsize  = self.gl.global_shape(scales=AA)[[dz]]
        slices = self.gl.slices(scales=AA)
        slices = (slices[dz],)
        cn     = [self.szA.cn[dz]]
        
        start  = self.gl.start(scales=AA)
        
        color  = np.int(np.all(lsize[[dz]]>0))
        lsize  = [lsize[dz]]
        c=new_comm(ju.comm,color, ju.mpirank)
        if color==1:
            p=self.check_dir(self.dir / fn)
            f=h5py.File(str(p), 'w', driver='mpio', comm=c)
            z=f.create_dataset(cn[0], (gsize[0],), dtype=np.float)
            if dim==3:
                if  start[1]==0 and start[0]==0: z[slices[0]]=np.reshape(self.grids(AA)[2][0,0,:],(lsize[0]))
            else:
                if                  start[0]==0: z[slices[0]]=np.reshape(self.grids(AA)[1][0,:],(lsize[0]))
            if type(nm)==list or type(nm)==tuple:
                for (a,dd) in zip(nm, data):
                    d  = f.create_dataset(a,gsize, dtype=np.float)
                    with d.collective: d[slices]=dd.reshape(lsize)
            else:
                d  = f.create_dataset(nm,gsize, dtype=np.float)
                with d.collective: d[slices]=data.reshape(lsize)
            f.close()
                        
    def check_dir(self,p):
        if ju.mpirank==0:
            if not p.parent.is_dir(): p.parent.mkdir()
        ju.Barrier()
        return(p)

    def savexyz(self,fn,nm,data,AA):
        debug_savexyz=False
        gsize  = self.gl.global_shape(scales=AA)
        lsize  = self.gl.local_shape(scales=AA)
        start  = self.gl.start(scales=AA)
        slices = self.gl.slices(scales=AA)
        color  = np.all(lsize>0)
        c=new_comm(ju.comm,np.int(color), ju.mpirank)
        fnn=self.check_dir(self.dir / fn)
        #if data is None:
        #    ju.logger_info('jrec.savexyz: data is None fn: {}, nm: {}, AA: {}'.format(fn,nm,AA))
        #    ju.set_error('jrec.savexyz: data is None fn: {}, nm: {}, AA: {}'.format(fn,nm,AA))
        #ju.check_error()
        if debug_savexyz: ju.logger.info('jrec.savexyz: fn: {}, nm: {}, AA: {}, data.shape: {}'.format(fn,nm,AA,data.shape))
        if color:
            with h5py.File(str(fnn), 'w', driver='mpio', comm=c) as f:
                if type(nm)==list or type(nm)==tuple:
                    for (a,b) in zip(nm, data):
                        d  = f.create_dataset(a,gsize, dtype=np.float)
                        with d.collective: d[slices]=b
                else:
                    d  = f.create_dataset(nm,gsize, dtype=np.float)
                    #ju.logger.info('type(data) {}, slices'.format(type(data),slices))
                    #ju.logger.info('data.shape() {}'.format(data.shape))
                    with d.collective: d[slices]=data
                   
    def print(self):
        print('jrec rank:{}/{}'.format(ju.mpirank,ju.mpisize))
        print(' start   {}'.format(self.start))
        print(' gsize   {}'.format(self.gsize))
        print(' lsize   {}'.format(self.lsize))
        print(' slices  {}'.format(self.slices))
        print(' xslice  {}'.format(self.xslice))
        print(' yslice  {}'.format(self.yslice))
        print(' zslice  {}'.format(self.zslice))
        print(' x0      {}'.format(self.x0))
        print(' y0      {}'.format(self.y0))
        print(' z0      {}'.format(self.z0))
        print(' nx      {}'.format(self.nx))
        print(' ny      {}'.format(self.ny))
        print(' nz      {}'.format(self.nz))
        print(' AA      {}'.format(self.AA))
        print(' AAS     {}'.format(self.AAS))
        print(' AAJ     {}'.format(self.AAJ))
        ju.Barrier()

        
        
  
        #        for axis, basis in enumerate(domain.bases):
        #            print('axis {}'.format(axis)) 0,2,3
        #            print('basis {}'.format(basis)) x,y,z
        #            print('basis.name {}'.format(basis.name)) x,y,z
        #            print('basis.element_label {}'.format(basis.element_label)) k
        #            print('basis.grid(AA) {}'.format(basis.grid(AA)))
        # domain.bases[j].name
 
        #self.var['x']=domain.bases[0].grid(AA)
        #self.var['z']=domain.bases[self.dim-1].grid(AA)


    # Need to set up a communicator with non-zero nodes for unaliased data
    # some systems hang when the slice is zero
    #gl=domain.dist.grid_layout.
    # impliment a jsave object for the program state make it the same as single field?
    
    # Need to be very careful not to load or save data from nodes that have no data
    # It works on some implimentations of parallel hdf5 but not others so we setup a separate communicator
    # Ideally all the communicators would be setup in the overarching jrec classe
    # with everything else subclasses of this so that each commuincator only needs to be setup once
    
    # add Nx Ny Nz
    # add in everything we need to create the coordinate system
    # Need to know interval and SinCos, Fourier or Cheb 
    #ps=list([])
    #px=list([])
    #if ju.mpirank==0:
    #    ps = list(set(('Fx','Fy','Fz')) & set(state_ev.keys()))   # Scalar parameters
    #px = list(set(('qvx',))         & set(state_ev.keys()))   # 'x' only parameters
    #ps=ju.comm.bcast(ps, root=0)        
    #px=ju.comm.bcast(px, root=0)
    #self.ps=dict()
    #self.px=dict()
    #for a in ps: self.ps[a]=state_ev[a]
    #for a in px: self.px[a]=state_ev[a]
    #ju.logger.info('jrec.ssolver: self.ps  {}'.format(self.ps.keys()))
    #ju.logger.info('jrec.ssolver: self.px  {}'.format(self.px.keys()))
    #ju.logger.info('jrec.ssolver: state_ev {}'.format(state_ev.keys()))

    def get_parity(self,y):
        jt=ju.jtype(y)
        
        if   jt=='fxyz':
            if self.dim==3: cn='xyz'
            else: cn='xz'
        elif jt=='fyz':  cn='yz'
        elif jt=='fx':   cn='x'
        n=len(cn)
        p=np.zeros((n,),dtype=np.int)
        for k in range(n):
            m=y.meta[cn[k]]
            if 'parity' in m: p[k]=m['parity']
        return(p)
    
    def save_field(self,tg,a,y,AA):
        debug_save_field=False
        if type(tg) is str:
            fn=self.check_dir(self.dir / tg)
            tg=h5py.File(fn, 'w', driver='mpio', comm=self.scomm)
        else: fn=None
        y.set_scales(AA)
        y['g']
        if debug_save_field: ju.logger.info('jrec.save_field: {}, AA:{}, shape:{}'.format(a,AA,x.shape))
        sl=self.gl.slices(scales=AA)
        if tg is not None:
            d = tg.create_dataset(a,self.gl.global_shape(scales=AA), dtype=np.float)
            with d.collective: d[sl]=y['g']
            d.attrs['parity']=self.get_parity(y) 
        if fn is not None: tg.close()
    def load_field_slow(self,tg,a,y,AA): # very slow but why ?
        #ju.logger.info('jrec.load_field: {:6s} color:{} AA:{}'.format(a,self.color,AA))
        y.set_scales(AA)
        y['g']
        #ju.logger.info('y[g].shape {}'.format(y['g'].shape))
        sl=self.gl.slices(scales=AA)
        #print('rank: {} sl {}, ls {}, AA {} ygs {} {}'.format(ju.mpirank,sl,ls,AA,y['g'].shape,self.color))
        if tg is not None:
            if a in tg:
                d = tg[a]
                if not np.all(y['g'].shape==d[sl].shape):
                    ju.logger.error('load_field: shape mismatch rank: {} yg.shape {}, d[sl].shape {} sl {}'.format(ju.mpirank,y['g'].shape,d[sl].shape,sl))
                with d.collective: y['g']=d[sl]
                py=self.get_parity(y)
                pd=d.attrs['parity']
                if not np.all(pd==py): ju.logger.info('jrec.load_field: {:6s} parity error {} {}'.format(a,pd,py))
                #else:                  ju.logger.info('jrec.load_field: {:6s} {} {}'.format(a,pd,py))
            else:
                ju.logger.info('jrec.load_field "{}" does not exist'.format(a))
                y['g']=0

    def load_field(self,tg,a,y,AA,param,sl1=None,sl2=None,rsx=False,rsy=False):
        y.set_scales(AA)
        y['g']
        jtp=ju.jtype(y)
        if sl2 is None: dim2=self.dim
        else:           dim2=len(sl2)
        if sl1 is None: dim1=dim2
        else:           dim1=len(sl1)
        #if jtp=='fyz':  ju.logger.info('load_field: type({}), {}.shape: {}, sl1: {}, sl2:, {}'.format(a,jtp,y['g'].shape,sl1,sl2)) 
        if   dim1==3 and dim2==2: py=self.get_parity(y)[[0,2]]
        elif dim1==2 and dim2==1: py=self.get_parity(y)[[0]]
        else:                     py=self.get_parity(y)
        if tg is not None:
            if a in tg:
                if sl2 is None: sl2=self.gl.slices(scales=AA)
                #if jtp=='fyz': sl2[0]=1
                #if jtp=='fyz' and sl1 is not none: sl1[0]=1
                d = tg[a]
                if rsx:
                    if self.dim==2:
                        with d.collective: y['g']=ju.fft_resample(d[:,sl2[1]],       self.gl.local_shape(scales=AA)[0],py[0])
                    else:
                        with d.collective: y['g']=ju.fft_resample(d[:,sl2[1],sl2[2]],self.gl.local_shape(scales=AA)[0],py[0])
                elif sl1 is None:#ju.logger.error('#:{}, {}.shape: {}, sl2:, {}, d[sl2].shape: {}, d.shape {}'.format(ju.mpirank,a,y['g'].shape,sl2,d[sl2].shape,d.shape)) 
                    with d.collective: y['g']=d[sl2]
                elif dim1==dim2:
                    with d.collective: y['g'][sl1]=d[sl2]
                elif dim1==3 and dim2==2:
                    with d.collective: y['g'][sl1]=np.expand_dims(d[sl2],1)
                else:
                    ju.logger.info('jrec.load_field: Incompatible slices {} {}'.format(sl1,sl2))
                    ju.set_error('Incompatible slices')
                    return
                if 'parity' in d.attrs:
                    pd=d.attrs['parity']
                    if not np.all(pd==py):
                        ju.logger.info('jrec.load_field "{}" parity error {} {}'.format(a,pd,py))
                        #ju.logger.info('jrec.load_field: {:6s} {} {}'.format(a,pd,py))
            else:
                ju.logger.info('jrec.load_field "{}" does not exist'.format(a))
                y['g']=0
        if sl1 is not None and param['exty']=='copy' and dim2==3 and dim1==3:
            if self.fcomm is None:
                slices = self.gl.slices(scales=AA)
                shape  = self.gl.local_shape(scales=AA) # find which node have the same x and z
                #print('rank: {}, slices: {}, shape: {}, exty; {}'.format(ju.mpirank,slices,shape,param['exty']))
                # make a communicator with all those that have the same z coordinates and same periodic y copies
                rshape = ju.comm.bcast(shape,root=0)
                if np.any(rshape != shape): # Check that they all have the same number of y grid points
                    print('jrec.load_field shapes on all nodes must match rank: {} shape: {}'.format(ju.mpirank,shape))
                    ju.set_error('jrec.load_field shapes do not match');
                ju.check_error()
                ny=np.int(slices[1].start/shape[1])
                nz=np.int(slices[2].start/shape[2])
                if tg is None: nny=0
                else:          nny=ny
                maxny=1+ju.comm.allreduce(nny, op=ju.MPI.MAX)
                maxnz=1+ju.comm.allreduce(nny, op=ju.MPI.MAX)
                nn=(ny%maxny)+maxny*nz
                self.oy = np.int(ny / maxny) % 2
                if np.int(ny/maxny)==0: rank=0
                else:                  rank=ju.mpirank
                self.fcomm = ju.comm.Split((ny%maxny)+maxny*nz,rank)
                #print('ju.rank: {:3d}, rank: {:3d}, ny: {:2d}, maxny: {:2d}, nz: {:2d}, maxnz: {:2d}, nn: {:3d}, oy: {:1d}'.format(ju.mpirank,rank,ny,maxny,nz,maxnz,nn,oy))

           
            y['g'] = self.fcomm.bcast(y['g'],root=0)
            if py[1] != 0 and self.oy == 1:
                if   py[1] == 1: y['g'] = np.flip( y['g'],axis=1)   # even
                elif py[1] ==-1: y['g'] = np.flip(-y['g'],axis=1)   # odd
                #if   py[1] == 1: print('jrec.load_field:rank: {} flipping even'.format(ju.mpirank))
                #elif py[1] ==-1: print('jrec.load_field:rank: {} flipping odd,'.format(ju.mpirank))
             
    # Can remove from here saving any variables that don't have time derivatives since these can be calculated, ie ux, vz, wz p
    def save_state(self,dt,state_ev,fn=None):
        debug_save_state=False
        if debug_save_state: ju.logger.info('jrec.save_state: dt: {}, fn:,{}'.format(dt,fn))
        self.print_range('jrec.save_state: ',self.solver.state,state_ev)

        tg=None
        wt=np.float(time.time())
        wt=ju.comm.bcast(wt, root=0)
        fnn = Path('')
        savenm = ju.comm.bcast(self.solver.state.field_names,root=0)
        for a in savenm:
            self.solver.state[a].set_scales(self.AAS)
            self.solver.state[a]['g']
            #self.savexyz('state-'+a+'.hdf5',a,self.solver.state[a]['g'],self.AAS)
        nm0    = list()
        nmf    = list()
        nmfx   = list()
        nmfyz  = list()
        nmfxyz = list()
        if ju.mpirank==0:
            for a in state_ev.keys():
                typ=ju.jtype(state_ev[a])
                if   typ=='fyz' or typ=='fxyz': nmf.append(a)
                if   typ=='fx':   nmfx.append(a)
                elif typ=='fyz':  nmfyz.append(a)
                elif typ=='fxyz': nmfxyz.append(a)
                else:             nm0.append(a)
        nm0    = ju.comm.bcast(nm0)
        nmf    = ju.comm.bcast(nmf)
        nmfx   = ju.comm.bcast(nmfx)
        nmfyz  = ju.comm.bcast(nmfyz)
        nmfxyz = ju.comm.bcast(nmfxyz)
        for a in nmf:
            state_ev[a].set_scales(self.AAS)
            state_ev[a]['g']
       
        if ju.mpirank==0:
            if fn is None: fn  = self.dir / 'final' / 'state-{:05d}.hdf5'.format(self.write_number)
            else: fn=self.dir / fn
            fno = Path(str(fn) + '.old')
            fnn = Path(str(fn) + '.new')
                        
            #ju.logger.info('jrec.save_state "{:s}"'.format(str(fn)))
            with h5py.File(str(fnn), 'w') as f:
                sg = f.create_group('scales')
                sg.create_dataset('interval',     data=self.interval)
                sg.create_dataset('gsize',        data=self.gsize)
                sg.create_dataset('iteration',    data=self.solver.iteration )
                sg.create_dataset('sim_time',     data=self.solver.sim_time  )
                sg.create_dataset('timestep',     data=dt)
                sg.create_dataset('wall_time',    data=wt)
                sg.create_dataset('write_number', data=self.write_number)
                sg.create_dataset('type',         data=self.btype)
                sg.create_dataset('AA',           data=self.AA)
                sg.create_dataset('AAS',          data=self.AAS)
                tg = f.create_group('tasks')
                for a in nm0:  tg.create_dataset( a, data=state_ev[a])
  
        if debug_save_state: ju.logger.info('jrec.save_state: saving fields')

        fnn=ju.comm.bcast(fnn, root=0)
        if self.color:
            with h5py.File(str(fnn), 'a', driver='mpio', comm=self.scomm) as f:
                tg = f['tasks']
                for a in nmfxyz+savenm:
                    d = tg.create_dataset(a, self.gsize, dtype=np.float)
                    with d.collective: d[self.slices]=self.solver.state[a]['g']
                    d.attrs['parity']=self.get_parity(self.solver.state[a]) 

        ju.Barrier()
        if ju.gyz.color:
            with h5py.File(str(fnn), 'a', driver='mpio', comm=ju.gyz.comm) as f:        
                tg = f['tasks']
                for a in nmfyz:
                    d  = tg.create_dataset(a,ju.gyz.gsize, dtype=np.float)
                    with d.collective: d[ju.gyz.slices]=state_ev[a]['g'].reshape(ju.gyz.lsize)
                    d.attrs['parity']=self.get_parity(state_ev[a]) 
           
        ju.Barrier()
        if ju.gx.color:
            with h5py.File(str(fnn), 'a') as f:        
                for a in nmfx:
                    d=f['tasks'].create_dataset(a,data=state_ev[a]['g'])
                    d.attrs['parity']=self.get_parity(state_ev[a]) 
           
        ju.Barrier()
        if debug_save_state: ju.logger.info('jrec.save_state: renaming files')
        if ju.mpirank==0:
            if fn.is_file(): fn.rename(fno)
            fnn.rename(fn)
            self.fno=fno
            ju.logger.info('jrec.save_state "{:s}" ({:8.3f} s) i={} dt={:8.6f}'.format(str(fn),np.float(time.time()-wt),self.solver.iteration,dt))
        ju.Barrier()      
        if debug_save_state: ju.logger.info('jrec.save_state: finished')

        
    def load_state(self,fn,solver,state_ev,param):
        wt=np.float(time.time())
        dt=np.float(0)      
        gsz=self.szA.gsize
        solver.iteration=np.int(0)
        solver.write_number=np.int(0)
        solver.sim_time=np.float(0)
        failed=0
        gsz=None
        tg=None
        ts=None
        
        sl1=None
        sl2=None
        dsl1=None
        dsl2=None
        ichange  = False
        szchange = False
        rsx=False
        rsy=False
        AA=1
        tgnm=None
        if ju.mpirank==0:
            with h5py.File(str(fn), 'r') as f:
                ts=f['scales']
                tg=f['tasks']
                if 'gsize' in ts: gsz = ts['gsize'][:]
                if 'timestep'     in ts: dt                  = ts['timestep'][()]
                if 'iteration'    in ts: solver.iteration    = ts['iteration'][()]
                if 'sim_time'     in ts: solver.sim_time     = ts['sim_time'][()]
                if 'write_number' in ts: solver.write_number = ts['write_number'][()]
                if 'type'         in ts: btype               = ts['type'][()]
                if 'interval'     in ts: interval            = ts['interval'][()]
                else:
                    ju.logger.info('jrec.load_state: {} does not contain interval data'.format(fn))
                    interval=self.interval
                if 'db' in state_ev and 'bp' in tg:
                    ju.logger.info('Renaming db to bp')
                    state_ev['bp'] = state_ev.pop('db')
                tgnm = list(set(tg) & set(state_ev.keys()))
                tgnmc=tgnm.copy()
                for a in tgnmc:
                    jt=ju.jtype(state_ev[a])
                    if jt=='scalar':
                        state_ev[a]=tg[a][()]
                        #ju.logger.info('jrec.load_state: sev: {:6s} {:7.4f} scalar'.format(a,state_ev[a]))
                        tgnm.remove(a)
                    elif jt=='ndarray':
                        f=np.array(tg[a])
                        sz1=state_ev[a].shape
                        sz2=f.shape
                        nsz1=len(sz1)
                        nsz2=len(sz2)
                        if np.prod(sz1)==np.prod(sz2): state_ev[a]=f.reshape(sz1)
                        elif sz1[0]==sz2[0]:
                            n=min(sz1[1],sz1[2])
                            state_ev[a][:,0:n]=np.array(tg[a])[:,0:n]
                        elif a in ('divz','bp','db','dd'):
                            state_ev[a]=ju.cheb_resample(np.array(tg[a]),sz1[0]).reshape(sz1)
                            ju.logger.info('jrec.load_state: cheb resample: {:6s} [{} to {}]'.format(a,sz2[0],sz1[0]))
                        elif a in ('divx'):
                            state_ev[a]=ju.pchipi_resample(np.array(tg[a]),sz1[0]).reshape(sz1)
                            ju.logger.info('jrec.load_state: pchip resample: {:6s} [{} to {}]'.format(a,sz2[0],sz1[0]))
                        elif a in ('topd'):
                            state_ev[a]=ju.fft_resample(ju.fft_resample(np.array(tg[a]),sz1[0],axis=0),sz1[1],axis=1).reshape(sz1)
                        else: ju.set_error('jrec.load_state: shapes for field "{}" differ {} {}'.format(a,sz1,sz1))
                        #ju.logger.info('jrec.load_state: sev: {:6s} [{:7.4f} {:7.4f}] {} {}'.format(a,state_ev[a].min(),state_ev[a].max(),sz1,sz2))
                        tgnm.remove(a)
                    elif jt=='fx':
                        f=np.array(tg[a])
                        sz1=state_ev[a]['g'].shape[0]
                        sz2=f.shape[0]
                        if sz1==sz2:
                            state_ev[a]['g']=f
                        else:
                            state_ev[a].set_scales(sz2/sz1)
                            state_ev[a]['g']=f
                            state_ev[a].set_scales(1)
                        tgnm.remove(a)

            if gsz is None:
                ju.logger.info('jrec.load_state: scales/gsize does not exist {}'.format(str(fn)))
                gsz=self.sz1.gsize
            
            dr1=(self.interval[:,1]-self.interval[:,0])/self.sz1.gsize
            dr2=(     interval[:,1]-     interval[:,0])/gsz
            dim1=len(dr1)
            dim2=len(dr2)
            dsl2=np.arange(0,dim2)
            if dim1==dim2:            dsl1=np.arange(0,dim1)
            elif dim1==3 and dim2==2: dsl1=np.array([0,2])
            else:
                ju.logger.info('jrec.load_state: Cannot change dimensions from {} to {}'.format(dim2,dim1))
                ju.set_error('dimensions changed')

            ichange  = not np.all(interval[dsl2] == self.interval[dsl1])
            szchange = not np.all(     gsz[dsl2] == self.sz1.gsize[dsl1])
            
            if dim1==dim2:
                if not np.all(btype==self.btype):
                    ju.logger.info('jrec.load_state: basis type changed {}'.format(fn))
                    ju.logger.info('           from: {}'.format(btype))
                    ju.logger.info('             to: {}'.format(self.btype))
                    ju.set_error('basis type changed')
            else:
                ichange=True
                ju.logger.info('jrec.load_state: Number of dimensions has changed from {} to {}'.format(dim2,dim1))
 
            if ichange: 
                AA=Fraction(np.int64(np.round(dr1[0]*self.sz1.gsize[0]/dr2[0])),np.int64(self.sz1.gsize[0]))
                ju.logger.info('jrec.load_state: interval changed {}'.format(fn))
                ju.logger.info('                 [{}] [{}] '.format(interval[:,0],self.interval[:,0]))
                ju.logger.info('                 [{}] [{}] '.format(interval[:,1],self.interval[:,1]))
                ju.logger.info('             AA: {}'.format(AA))
                ju.logger.info('            dr1: {}'.format(dr1))
                ju.logger.info('            dr2: {}'.format(dr2))
                ju.logger.info('            sz1: {}'.format(self.sz1.gsize))
                ju.logger.info('            gsz: {}'.format(gsz))
                sl1=list()  # slices into local copy of field
                sl2=list()  # slices into field stored in datafile
                if dim2==2 and dim1==3: # changing from 2 to 3 dimensions
                    i0=(max(self.interval[0,0],interval[0,0]),max(self.interval[2,0],interval[1,0]))
                    i1=(min(self.interval[0,1],interval[0,1]),min(self.interval[2,1],interval[1,1]))
                    # Should this be dr1 not dr2 ???
                    sl1.append(slice( np.int64(max(i0[0],self.interval[0,0])/dr2[0]),np.int64(min(self.interval[0,1],i1[0])/dr2[0]) ))
                    sl2.append(slice( np.int64(max(i0[0],     interval[0,0])/dr2[0]),np.int64(min(     interval[0,1],i1[0])/dr2[0]) ))
                    sl1.append(slice(0,np.int64(AA*self.szS.lsize[1])))
                    sl1.append(slice( np.int64(max(i0[1],self.interval[2,0])/dr2[1]),np.int64(min(self.interval[2,1],i1[1])/dr2[1]) ))
                    sl2.append(slice( np.int64(max(i0[1],     interval[1,0])/dr2[1]),np.int64(min(     interval[1,1],i1[1])/dr2[1]) ))
                elif dim1==dim2: 
                    for k in range(self.dim):
                        i0=max(self.interval[k,0],interval[k,0])
                        i1=min(self.interval[k,1],interval[k,1])
                        #ju.logger.info(' k {}       [{:6.4f} {:6.4f}]'.format(k,i0,i1))
                        #ju.logger.info('     state [{:6.4f} {:6.4f}]'.format(interval[k,0],interval[k,1]))
                        #ju.logger.info(' new state [{:6.4f} {:6.4f}]'.format(self.interval[k,0],self.interval[k,1]))
                        sl1.append(slice( np.int64(max(i0,self.interval[k,0])/dr2[k]),np.int64(min(self.interval[k,1],i1)/dr2[k]) ))
                        sl2.append(slice( np.int64(max(i0,     interval[k,0])/dr2[k]),np.int64(min(     interval[k,1],i1)/dr2[k]) ))
                else:
                    ju.logger.info('jrec.load_state: Cannot go from dim2:{} to dim1:{}'.format(dim2,dim1))
                    ju.set_error('Mismatched sizes')
   
                ju.logger.info('jrec.load_state: AA:{}, sz1: {}, sz2: {}'.format(AA,np.int64(gsz/AA),gsz))
                #ju.logger.info('            sl1: {}'.format(sl1))
                #ju.logger.info('            sl2: {}'.format(sl2))
                # These are not correct for extending chebbychev polynomials
            else:
                if np.all(gsz==self.sz1.gsize):
                    AA=1
                    ju.logger.info('jrec.load_state: saved is unaliased  AA:{}, sz: {}'.format(AA,gsz))
                elif np.all(gsz==self.szA.gsize):
                    AA=self.AA
                    ju.logger.info('jrec.load_state: saved is aliased  AA:{}, sz: {}'.format(AA,gsz))
                else:
                    AA=Fraction(np.int64(gsz[-1]),np.int64(self.sz1.gsize[-1]))
                    ju.logger.info('')
                    ju.logger.info('jrec.load_state: Resizing using set_scales(), AA: {}'.format(AA))
                    ju.logger.info('jrec.load_state: saved size: {}'.format(gsz))
                    ju.logger.info('jrec.load_state: sim   size: {} '.format(self.sz1.gsize))
                    ju.logger.info('jrec.load_state: load  size: {}'.format(np.array(AA*self.sz1.gsize,dtype=np.int64)))
                    rsx=np.int64(AA*self.sz1.gsize[0]) != gsz[0] # Resample x
                    if dim1==3 and dim2==3: rsy=np.int64(AA*self.sz1.gsize[1]) != gsz[1] # Resample y
                    else: rsy = False
                    if rsx: ju.logger.info('jrec.load_state: Resample x')
                    if rsy: ju.logger.info('jrec.load_state: Resample y')
                    if np.int64(AA*self.sz1.gsize[dsl1[1]])!=gsz[dsl2[1]]:
                        ju.logger.info('jrec.load_state: sizes for {} do not match AA {} s'.format(str(fn),AA))
                        ju.set_error('Mismatched sizes')
        ju.check_error()
 
        rsx                                               = ju.comm.bcast(rsx,  root=0)
        rsy                                               = ju.comm.bcast(rsy,  root=0)
        dsl1                                              = ju.comm.bcast(dsl1, root=0)
        dsl2                                              = ju.comm.bcast(dsl2, root=0)
        tgnm                                              = ju.comm.bcast(tgnm, root=0)
        for a in state_ev: # don't broadcast fields
            if  np.isscalar(state_ev[a]) or isinstance(state_ev[a], np.ndarray):
                state_ev[a]                               = ju.comm.bcast(state_ev[a], root=0)
        AA                                                = ju.comm.bcast(AA,root=0)   
        sl1                                               = ju.comm.bcast(sl1,root=0)   
        sl2                                               = ju.comm.bcast(sl2,root=0)   
        solver.iteration    = solver.initial_iteration    = ju.comm.bcast(solver.iteration,    root=0)  
        solver.sim_time     = solver.initial_sim_time     = ju.comm.bcast(solver.sim_time,     root=0)  
        solver.write_number = solver.initial_write_number = ju.comm.bcast(solver.write_number, root=0)  
        dt                                                = ju.comm.bcast(dt,                  root=0)  
        ichange                                           = ju.comm.bcast(ichange,                  root=0)  
        
        ls=copy.copy(self.gl.local_shape(scales=AA)) # The copy here is critical
        if ichange:
            #y.set_scales(AA)
            # Set the scale of each dimension
            sl=self.gl.slices(scales=AA)
            #ju.logger.info('jrec.load_state: sl ={}'.format(sl))
            #ju.logger.info('jrec.load_state: sl1={}'.format(sl1))
            #ju.logger.info('jrec.load_state: sl2={}'.format(sl2))
            for k in range(len(dsl1)):
                k1=dsl1[k]
                k2=dsl2[k]
                sl2[k2]=slice(max(sl2[k2].start,sl[k1].start),min(sl2[k2].stop,sl[k1].stop))
                ls[k2]=sl2[k2].stop-sl2[k2].start
                sl1[k1]=slice(max(sl1[k1].start,0),max(sl1[k1].start,0)+ls[k2])
            #ju.logger.info('jrec.load_state: sl1={}'.format(sl1))
            #ju.logger.info('jrec.load_state: sl2={}'.format(sl2))
            #quit(1)
            sl1=tuple(sl1)
            sl2=tuple(sl2)
        if AA==self.AAS and not ichange:
            gs=self.gsize
            color=self.color
            scomm=self.scomm
        else:
            gs=self.gl.global_shape(scales=AA)
            color=np.all(ls>0)
            scomm=new_comm(ju.comm,np.int(color), ju.mpirank)   

            
        # All nodes must be called so that transforms and transposes can occur even if they then
        if color:
            f=h5py.File(str(fn), 'r', driver='mpio', comm=scomm)
            tg=f['tasks']
        else:
            tg=None
        tgnmc=tgnm.copy()
        for a in tgnmc:
            if ju.jtype(state_ev[a])=='fxyz':
                ju.logger.info('jrec.load_state: state_ev[{}] fxyz'.format(a))
                self.load_field(tg, a,     state_ev[a], AA,param,sl1,sl2,rsx,rsy)
                tgnm.remove(a)
        for a in solver.state.field_names:
            self.load_field(tg, a, solver.state[a], AA,param,sl1,sl2,rsx,rsy)
        if color: f.close()
        
        if ju.gyz.color:
            f=h5py.File(str(fn), 'r', driver='mpio', comm=ju.gyz.comm)
            tg=f['tasks']
        else:
            tg=None
        eyz=None
        tgnmc=tgnm.copy()
        for a in tgnmc:
            if ju.jtype(state_ev[a])=='fyz':
                tgnm.remove(a)
                if AA==1: self.load_field(tg, a, state_ev[a], 1,param,sl1,ju.gyz.slices,rsx,rsy)
                else:
                    AA=Fraction(AA)
                    if eyz is None:
                        if AA>1: ed=ju.gridyz(ju.gyz.p,ju.gyz.A*AA,ju.gyz.m, Nx=AA._numerator)
                        else:    ed=ju.gridyz(ju.gyz.p,1, ju.gyz.m, Nx=AA._numerator)
                        #ed.print()
                        eyz=ed.domain.new_field()
                        eyz.meta['x']['parity']=0
                    eyz.meta['y']['parity']=state_ev[a].meta['y']['parity']
                    eyz.meta['z']['parity']=state_ev[a].meta['z']['parity']
                    
                        
                    ju.logger.info('jrec.load_state: state_ev[{}] fyz'.format(a))
                    #ju.logger.info('AA: {}/{}, sl1: {}, ed.slices: {} , rsx: {}, rsy: {} edsl: {}'.format(AA._numerator,AA._denominator,sl1,ju.gyz.slices,rsx,rsy,ed.domain.dist.grid_layout.slices(scales=1)))
                    eyz.set_scales(1)
                    self.load_field(tg, a, eyz, 1,param,sl1,ed.slices,rsx,rsy)
                    eyz.set_scales(1/AA)
                    state_ev[a]['g']=eyz['g'][0:1,...]
        if ju.gyz.color: f.close()

        if 'bp' in state_ev:  state_ev['db']   = state_ev.pop('bp')
  
        if ju.mpirank==0:
            for a in tgnm: ju.logger.info('Failed to load {}'.format(a))
 
        ju.check_error()
        
        # This doesn't work because the boundary conditions are not correctly imposed 
        #for a in self.diffnm: eval("solver.state['{}'].differentiate('z', out=solver.state['{}'])".format(a[0:-1],a))
        #ju.logger.info('max uz error {}'.format(ju.jmax(np.abs(solver.state['uz']['g']-solver.state['u'].differentiate('z')['g']))))
        #ju.logger.info('max wz error {}'.format(ju.jmax(np.abs(solver.state['wz']['g']-solver.state['w'].differentiate('z')['g']))))
        #ju.logger.info('max bz error {}'.format(ju.jmax(np.abs(solver.state['bz']['g']-solver.state['b'].differentiate('z')['g']))))
        #quit(1)
        ju.logger.info('jrec.load_state: {} ({:8.3f} s) '.format(fn,np.float(time.time()-wt)))
        ju.logger.info("jrec.load_state: iteration:    {:7d}".format(solver.iteration))
        ju.logger.info("jrec.load_state: sim time:     {:7.2f}".format(solver.sim_time))
        ju.logger.info("jrec.load_state: write_number: {:7d}".format(solver.write_number))
        ju.logger.info("jrec.load_state: timestep:     {:7.5f}".format(dt))

        self.print_range('jrec.load_state: ',solver.state,state_ev)
        #    if np.ndims(state_ev['divz'])==1: np.reshape(state_ev['divz'],(state_ev['divz'].shape[0],1))
          
        return(dt,state_ev)
 
    def print_range(self,s,state,sev):
        for a in state.field_names: ju.print_v(a,state[a],s+'var ')
        for a in sev:               ju.print_v(a,sev[a],  s+'sev ')

                
   
    def end(self):
       if ju.mpirank==0:
            if self.fno.is_file(): self.fno.unlink()
