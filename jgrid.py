import h5py
from mpi4py import MPI
import ded_types as dp
from dedalus import public as de
import numpy as np


class gridx: # Only on mpirank 0
    def __init__(self,p,A,comm):
        self.jtype = 'fx'
        self.mpirank=comm.Get_rank()
        self.color = self.mpirank==0
        self.comm  = comm.Split(int(self.color), self.mpirank)

        if not self.color: return

        self.L  = p['Length']
        self.Nx = np.int(A*p['Nx'])
        self.Tx = p['Tx']
        self.A  = A
        self.LL = (0,self.L)
        self.Rx = self.L

        if   self.Tx=="SinCos":  self.xb=de.SinCos(   'x', self.Nx, interval=self.LL, dealias=1)
        elif self.Tx=="Fourier": self.xb=de.Fourier(  'x', self.Nx, interval=self.LL, dealias=1)
        elif self.Tx=="Cheb":    self.xb=de.Chebyshev('x', self.Nx, interval=self.LL, dealias=1)
        self.x  = self.xb.grid(scale=1)
        self.dx = np.mean(np.diff(np.asarray(self.x)))
        
        self.domain = de.Domain((self.xb,), grid_dtype=np.float64, comm=self.comm)
        self.gsize  = (self.Nx,)
        self.lsize  = (self.Nx,)
        self.slices = (slice(0,self.Nx),)

    def new_field(self,Px):
        if self.color:
            f=self.domain.new_field()
            if self.Tx == "SinCos":  f.meta['x']['parity'] = Px
            else:                    f.meta['x']['parity'] = 0
            return(f)
        else:          return(None)
        
    def create_dataset(self,x,nm,f):
        if not self.color: return
        f.create_dataset(nm,data=x['g'])
        
    def write(self,fn,x,m,vnm=None):
        if not self.color: return
        if vnm is None: vnm=x.keys()
        if len(vnm)==0: return
        with h5py.File(fn, m) as f:
            for a in vnm:
                d = f.create_dataset(a,data=x[a]['g'])
                d.attrs['parity']=self.get_parity(x[a]) 

    def get_parity(self,f):
        if not self.color: return
        p=np.zeros((1,),dtype=np.int)
        if 'parity' in f.meta['x']: p[0]=f.meta['x']['parity']
        return(p)
    
    def intyz(self,f,w,o,sl=None): # o=int(f*w dy dz)
        if sl is None:
            y = comm.reduce(np.sum(w*f,         axis=(1,2)), op=MPI.SUM)
            if self.color: o['g'] = y 
        else:
            y = comm.reduce(np.sum(w*f[sl,...], axis=(1,2)), op=MPI.SUM)
            if self.color: o['g'][sl] = y

    def intyz2(self,f1,w1,f2,w2,o,sl=None): # o=int(f*w dy dz)
        if sl is None:
            y = comm.reduce(np.sum(w1*f1,axis=(1,2))+np.sum(w2*f2,axis=(1,2)), op=MPI.SUM)
            if self.color: o['g'] = y 
        else:
            y = comm.reduce(np.sum(w1*f1[sl,...], axis=(1,2)), op=MPI.SUM)+comm.reduce(np.sum(w2*f2[sl,...], axis=(1,2)), op=MPI.SUM)
            if self.color: o['g'][sl] = y

    def set_parity(self,f,Px):
        if not self.color: return
        if self.Tx == "SinCos":  f.meta['x']['parity'] = Px
        else:                    f.meta['x']['parity'] = 0

    def print_v(self,x,s=''):
        #logger.error('{}, gridx.print_v: rank: {}, self.color {}, type(x): {}'.format(s,self.mpirank,self.color,type(x)))
        #logger.info('{} [{:9.3f},{:9.3f}]'.format(s,jmin(x),jmax(x)))
        #logger.error('gridx.print_v: 1 {}'.format(self.mpirank))
        if self.color: logger.info('{}[{:9.3f},{:9.3f}]'.format(s,x['g'].min(),x['g'].max()))
        #logger.error('gridx.print_v: 2 {}'.format(self.mpirank))
        #Barrier()
        #quit()
        #if self.color: logger.info('{} [{:9.3f},{:9.3f}]'.format(s,jmin(x,self.comm),jmax(x,self.comm)))
        #logger.error('gridx.print_v: 3 {}'.format(self.mpirank))
        #Barrier()
        #quit()
        #logger.error('gridx.print_v: 3')
        #logger.info('{} [{:9.3f},{:9.3f}]'.format(s,jmin(x),jmax(x)))
        #logger.error('gridx.print_v: 4')

    def integrate(self,f):
        if self.color: return(f.integrate('x')['g'][0])
        else:          return(None)

    def differentiate(self,f,out):
        if self.color: out['g']=f.differentiate('x')['g']

class gridyz:
    def __init__(self,p,A,comm,m,Nx=1):
        self.mpirank=comm.Get_rank()
        self.jtype='fyz'
        self.p  = p
        self.A  = A
        self.m  = m
        self.W  = p['Width']
        self.H  = p['Height']
        self.Nx = Nx
        self.Ny = np.int(A*p['Ny'])
        self.Nz = np.int(A*p['Nz'])
        self.Ty = dp.get_param(p,'Ty')
        self.Tz = dp.get_param(p,'Tz')
        
        sType = p['sType']
        
        self.WW=(0,self.W)
        self.HH=(0,self.H)
        if sType=='pm':
            if self.Ty != "SinCos": self.WW=(-self.W/2,self.W/2)
            if self.Tz != "SinCos": self.HH=(-self.H/2,self.H/2)
            
            
        self.Ry=self.WW[1]-self.WW[0]
        self.Rz=self.HH[1]-self.HH[0]
        self.Area=self.Ry*self.Rz
        x_basis=de.Fourier(  'x', self.Nx, interval=(0,1), dealias=1)
        if   self.Ty=="SinCos":  self.yb=de.SinCos(   'y', self.Ny, interval=self.WW, dealias=1)
        elif self.Ty=="Fourier": self.yb=de.Fourier(  'y', self.Ny, interval=self.WW, dealias=1)
        elif self.Ty=="Cheb":    self.yb=de.Chebyshev('y', self.Ny, interval=self.WW, dealias=1)
        if   self.Tz=="SinCos":  self.zb=de.SinCos(   'z', self.Nz, interval=self.HH, dealias=1)
        elif self.Tz=="Fourier": self.zb=de.Fourier(  'z', self.Nz, interval=self.HH, dealias=1)
        elif self.Tz=="Cheb":    self.zb=de.Chebyshev('z', self.Nz, interval=self.HH, dealias=1)
        #print(vars(self))
        if self.Ty is not None:
            self.y = self.yb.grid(scale=1)
            self.dy=np.mean(np.diff(np.asarray(self.y)))
        if self.Tz is not None:
            self.z = self.zb.grid(scale=1)
            self.dz=np.mean(np.diff(np.asarray(self.z)))
        if self.Ty is not None:
            self.domain=de.Domain([x_basis, self.yb, self.zb], grid_dtype=np.float64, mesh=m)
        else:
            self.domain=de.Domain([x_basis, self.zb], grid_dtype=np.float64)
        gl=self.domain.dist.grid_layout
        self.sl1 = gl.slices(scales=1)
        if self.Ty is not None:
            self.gsize  = gl.global_shape(scales=1)[[1,2]]
            self.lsize  = gl.local_shape(scales=1)[[1,2]]
            self.slices = (self.sl1[1],self.sl1[2])
        else:
            self.gsize  = gl.global_shape(scales=1)[[1]]
            self.lsize  = gl.local_shape(scales=1)[[1]]
            self.slices = (self.sl1[1],)

        self.color  = np.all(self.lsize>0)
        self.comm=comm.Split(int(self.color), self.mpirank)

    def set_parity(self,f,Py,Pz):
        if self.Ty == "SinCos":  f.meta['y']['parity'] = Py
        else:                    f.meta['y']['parity'] = 0
        if self.Tz == "SinCos":  f.meta['z']['parity'] = Pz
        else:                    f.meta['z']['parity'] = 0

    def print(self):
        logger.info('A  {}'.format(self.A))
        logger.info('m  {}'.format(self.m)) 
        logger.info('WW {}'.format(self.WW)) 
        logger.info('HH {}'.format(self.HH)) 
        logger.info('Nx {}'.format(self.Nx))
        logger.info('Ny {}'.format(self.Ny))
        logger.info('Nz {}'.format(self.Nz))
        logger.info('Ty {}'.format(self.Ty))
        logger.info('Tz {}'.format(self.Tz))     
        logger.info('sl1 {}'.format(self.sl1))
        logger.info('gsize  {}'.format(self.gsize))
        logger.info('lsize  {}'.format(self.lsize))
        logger.info('slices {}'.format(self.slices))

    def create_dataset(self,x,nm,f):
        #logger.info('type(x): {}, x.shape: {}'.format(type(x),x.shape))
        if self.color:
            d  = f.create_dataset(nm,self.gsize, dtype=np.float)
            with d.collective: d[self.slices]=x.reshape(self.lsize)

    def get_parity(self,f):
        p=np.zeros((2,),dtype=np.int)
        if 'parity' in f.meta['y']: p[0]=f.meta['y']['parity']
        if 'parity' in f.meta['z']: p[1]=f.meta['z']['parity']
        return(p)

    def write(self,fn,x,m,vnm=None):
        #print('gridyz.write started: rank: {} fn: {}, m: {}, self.color: {}'.format(self.mpirank,fn,m,self.color))
        #traceback.print_stack(file=sys.stdout)
        if vnm is None: vnm=x.keys()
        if len(vnm)==0: return
        for a in vnm: x[a]['g']
        if self.color:
            with h5py.File(fn, m, driver='mpio', comm=self.comm) as f:
                for a in vnm:
                    d  = f.create_dataset(a,self.gsize, dtype=np.float)
                    with d.collective: d[self.slices]=x[a]['g'].reshape(self.lsize)
                    d.attrs['parity']=self.get_parity(x[a]) 
                    
    def intx(self,f,w,o,sl=None):
        if sl is None: o['g']=np.sum(w*f,        axis=(0,))
        else:          o['g']=np.sum(w*f[sl,...],axis=(0,)) 

    def new_field(self,Py,Pz):
        if self.color:
            f=self.domain.new_field()
            if self.Ty == "SinCos":  f.meta['y']['parity'] = Py
            else:                    f.meta['y']['parity'] = 0
            if self.Tz == "SinCos":  f.meta['z']['parity'] = Pz
            else:                    f.meta['z']['parity'] = 0
            return(f)
        else:          return(None)

class gridz:
    def __init__(self,p,A,comm,m,Nx=1,Ny=1):
        self.mpirank=comm.Get_rank()
        self.jtype='fz'
        self.p  = p
        self.A  = A
        self.m  = m
        self.H  = p['Height']
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = np.int(A*p['Nz'])
        Ty = dp.get_param(p,'Ty')
        self.Tz = dp.get_param(p,'Tz')
        
        sType = p['sType']
        
        self.WW=(0,self.W)
        self.HH=(0,self.H)
        if sType=='pm':
            if self.Tz != "SinCos": self.HH=(-self.H/2,self.H/2)
            
            
        self.Rz=self.HH[1]-self.HH[0]
        x_basis=de.Fourier(  'x', self.Nx, interval=(0,1), dealias=1)
        if   self.Tz=="SinCos":  self.zb=de.SinCos(   'z', self.Nz, interval=self.HH, dealias=1)
        elif self.Tz=="Fourier": self.zb=de.Fourier(  'z', self.Nz, interval=self.HH, dealias=1)
        elif self.Tz=="Cheb":    self.zb=de.Chebyshev('z', self.Nz, interval=self.HH, dealias=1)
        #print(vars(self))


        if self.Tz is not None:
            self.z = self.zb.grid(scale=1)
            self.dz=np.mean(np.diff(np.asarray(self.z)))
        if Ty is not None:
            y_basis=de.Fourier(  'y', self.Ny, interval=(0,1), dealias=1)
            self.domain=de.Domain([x_basis, y_basis, self.zb], grid_dtype=np.float64, mesh=m)
        else:
            self.domain=de.Domain([x_basis, self.zb], grid_dtype=np.float64)
        gl=self.domain.dist.grid_layout
        self.sl1 = gl.slices(scales=1)
        self.gsize  = gl.global_shape(scales=1)[[-1]]
        self.lsize  = gl.local_shape(scales=1)[[-1]]
        self.slices = (self.sl1[-1],)

        self.color  = np.all(self.lsize>0)
        self.comm=comm.Split(int(self.color), self.mpirank)

    def set_parity(self,f,Pz):
        if self.Tz == "SinCos":  f.meta['z']['parity'] = Pz
        else:                    f.meta['z']['parity'] = 0

    def print(self):
        logger.info('A  {}'.format(self.A))
        logger.info('m  {}'.format(self.m)) 
        logger.info('WW {}'.format(self.WW)) 
        logger.info('HH {}'.format(self.HH)) 
        logger.info('Nx {}'.format(self.Nx))
        logger.info('Ny {}'.format(self.Ny))
        logger.info('Nz {}'.format(self.Nz))
        logger.info('Tz {}'.format(self.Tz))     
        logger.info('sl1 {}'.format(self.sl1))
        logger.info('gsize  {}'.format(self.gsize))
        logger.info('lsize  {}'.format(self.lsize))
        logger.info('slices {}'.format(self.slices))

    def create_dataset(self,x,nm,f):
        #logger.info('type(x): {}, x.shape: {}'.format(type(x),x.shape))
        if self.color:
            d  = f.create_dataset(nm,self.gsize, dtype=np.float)
            with d.collective: d[self.slices]=x.reshape(self.lsize)

    def get_parity(self,f):
        p=np.zeros((2,),dtype=np.int)
        if 'parity' in f.meta['z']: p[0]=f.meta['z']['parity']
        return(p)

    def write(self,fn,x,m,vnm=None):
        #print('gridyz.write started: rank: {} fn: {}, m: {}, self.color: {}'.format(self.mpirank,fn,m,self.color))
        #traceback.print_stack(file=sys.stdout)
        if vnm is None: vnm=x.keys()
        if len(vnm)==0: return
        for a in vnm: x[a]['g']
        if self.color:
            with h5py.File(fn, m, driver='mpio', comm=self.comm) as f:
                for a in vnm:
                    d  = f.create_dataset(a,self.gsize, dtype=np.float)
                    with d.collective: d[self.slices]=x[a]['g'].reshape(self.lsize)
                    d.attrs['parity']=self.get_parity(x[a]) 
                    
    def intx(self,f,w,o,sl=None):
        if sl is None: o['g']=np.sum(w*f,        axis=(0,))
        else:          o['g']=np.sum(w*f[sl,...],axis=(0,)) 

    def new_field(self,Pz):
        if self.color:
            f=self.domain.new_field()
            if self.Tz == "SinCos":  f.meta['z']['parity'] = Pz
            else:                    f.meta['z']['parity'] = 0
            return(f)
        else:          return(None)

