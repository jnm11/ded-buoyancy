import numpy as np 
import copy as copy
from pathlib import Path
import jutil as ju
import h5py


    
class javrg:
    def __init__(self,solver,problem,domain,ddir,AA,reset,jdt):
        self.dim=len(domain.dist.grid_layout.global_shape(scales=1))
        #ju.logger.info('global_shape {}'.format(domain.dist.grid_layout.global_shape(scales=AA)))
        self.size   = domain.dist.grid_layout.global_shape(scales=AA)
        self.lsize  = domain.dist.grid_layout.local_shape(scales=AA)
        self.slices = domain.dist.grid_layout.slices(scales=AA)
        self.AA = AA
        self.start  = domain.dist.grid_layout.start(scales=AA)
        self.end    = domain.dist.grid_layout.start(scales=AA)
        self.jdt    = jdt
        if self.dim==3:
            self.nx=self.size[0]
            self.ny=self.size[1]
            self.nz=self.size[2]
            self.xslice = self.slices[0]
            self.yslice = self.slices[1]
            self.zslice = self.slices[2]
            self.x=copy.copy(domain.grids(AA)[0][:,0,0])
            self.y=copy.copy(domain.grids(AA)[1][0,:,0])
            self.z=copy.copy(domain.grids(AA)[2][0,0,:])

            self.sy     = self.start[1]==0
            if self.sy:
                self.lsize  = self.lsize[[0,2]]
                self.slices = (self.slices[0],self.slices[2])
                self.aslices = (self.slices[0],slice(0,1,None),self.slices[1])
                self.intyx  = self.x
                self.intyy  = self.y
                self.intyz  = self.z
                self.yint1slices = (self.xslice,self.zslice)
                self.yint2slices = (self.xslice,slice(0,1,None),self.zslice)
            else:
                self.yslice = slice(0,0,None)
                self.lsize  = [0,0]
                self.slices = (slice(0,0,None),slice(0,0,None))
                self.aslices = (slice(0,0,None),slice(0,0,None),slice(0,0,None))
                self.intyx  = np.array([],dtype=np.float)
                self.intyy  = np.array([],dtype=np.float)
                self.intyz  = np.array([],dtype=np.float)
                self.yint1slices = (slice(0,0,None),slice(0,0,None))
                self.yint2slices = (slice(0,0,None),slice(0,0,None),slice(0,0,None))
            self.nm = ('p','b','u','v','uu','vv','ww','bw','bb','Wx','Wy','Wz','S')  #'bu','uw',
        else:
            self.nx=self.size[0]
            self.ny=1
            self.nz=self.size[1]
            self.xslice = self.slices[0]
            self.yslice = slice(0,1,None)
            self.zslice = self.slices[1]
            self.x=copy.copy(domain.grids(AA)[0][:,0])
            self.y=np.array([0],dtype=np.float)
            self.z=copy.copy(domain.grids(AA)[1][0,:])
            self.sy = True
            self.nm = ('p','b','u','uu','ww','bw','bb','Wy','S')  #'bu','uw',
            self.yint1slices = (self.xslice,self.zslice)
            self.yint2slices = (self.xslice,self.zslice)
 

        self.intysize   = np.array([self.nx,self.nz],dtype=np.int)

                
          
        #print('self.start  {}'.format(self.start))
        #print('self.size   {}'.format(self.size))
        #print('self.lsize  {}'.format(self.lsize))
        #print('self.slices {}'.format(self.slices))
        #print('self.aslice {}'.format(self.aslices))
        
        self.solver=solver
        self.dir=Path(str(ddir)+'/javrg')
        if ju.mpirank==0:
            if not self.dir.is_dir(): self.dir.mkdir(parents=True)
        ju.comm.Barrier()
        if 'p' in problem.variables: self.P = solver.state['p']
        if 'b' in problem.variables: self.B = solver.state['b']
        if 's' in problem.variables: self.S = solver.state['s']
        if 'u' in problem.variables: self.U = solver.state['u']
        if 'v' in problem.variables: self.V = solver.state['v']
        if 'w' in problem.variables: self.W = solver.state['w']

        self.var=dict()
        #        for axis, basis in enumerate(domain.bases):
        #            print('axis {}'.format(axis)) 0,2,3
        #            print('basis {}'.format(basis)) x,y,z
        #            print('basis.name {}'.format(basis.name)) x,y,z
        #            print('basis.element_label {}'.format(basis.element_label)) k
        #            print('basis.grid(AA) {}'.format(basis.grid(AA)))
        # domain.bases[j].name
 
        #self.var['x']=domain.bases[0].grid(AA)
        #self.var['z']=domain.bases[self.dim-1].grid(AA)
        self.t2=copy.copy(solver.sim_time)
        self.i2=copy.copy(solver.iteration)
        self.clear()

        ju.logger.info('javrg: directory {}'.format(self.dir))
        self.write_num = -1
        if not reset and ju.mpirank==0:
            p = sorted(self.dir.glob('javrg-*.hdf5'),reverse=True)
            for a in p:
                with h5py.File(a,'r') as f:
                    #for b in f: print('   b: {}'.format(f[b]))
                    if 'write_num' in f:
                        self.write_num = f['write_num'][()]
                        ju.logger.info('javrg: file {}, write_num: {}'.format(a,self.write_num))
                        break
                    else:
                        ju.logger.info('javrg: file {} no write_num'.format(a))
        self.write_num +=1
        self.write_num = ju.comm.bcast(self.write_num, root=0)
        ju.logger.info('javrg: write_num {}'.format(self.write_num))

    def fint(self,x):
        if self.dim==3: x=x.integrate('y')
        return(np.reshape(x['g'][self.aslices],self.lsize))

    def fsqr(self,x):
        x=x**2
        x=x.evaluate()
        if self.dim==3: x=x.integrate('y')
        return(np.reshape(x['g'][self.aslices],self.lsize))

    def fprod(self,A,B):
        x=A*B
        x=x.evaluate()
        if self.dim==3: x=x.integrate('y')
        return(np.reshape(x['g'][self.aslices],self.lsize))

    def dsqr(self,A,B):
        x=(A-B)**2
        x=x.evaluate()
        if self.dim==3: x=x.integrate('y')
        return(np.reshape(x['g'][self.aslices],self.lsize))
            
    def clear(self):
        x=np.zeros(self.lsize, dtype=np.float64)
        for a in self.nm: self.var[a]  = copy.copy(x)
        self.n=0
        self.t1=copy.copy(self.t2)
        self.i1=copy.copy(self.i2)
        self.dt=0
         
    def update(self, dt):
        self.n  += 1
        self.dt += dt

        if 'p'  in self.var: self.var['p']  += dt*self.fint(self.P)
        if 'b'  in self.var: self.var['b']  += dt*self.fint(self.B)
        if 'u'  in self.var: self.var['u']  += dt*self.fint(self.U)
        if 'v'  in self.var: self.var['v']  += dt*self.fint(self.V)
        if 'uu' in self.var: self.var['uu'] += dt*self.fsqr(self.U)
        if 'vv' in self.var: self.var['vv'] += dt*self.fsqr(self.V)
        if 'ww' in self.var: self.var['ww'] += dt*self.fsqr(self.W)
        #if 'bu' in self.var: self.var['bu'] += dt*self.fprod(self.B,self.U)
        if 'bw' in self.var: self.var['bw'] += dt*self.fprod(self.B,self.W)
        if 'bb' in self.var: self.var['bb'] += dt*self.fprod(self.B,self.B)
        #if 'uw' in self.var: self.var['uw'] += dt*self.fprod(self.U,self.W)

      
        ux = self.U.differentiate('x')
        uz = self.U.differentiate('z')
        wx = self.W.differentiate('x')
        wz = self.W.differentiate('z')
        if self.dim==3:
            uy = self.U.differentiate('y')
            vx = self.V.differentiate('x')
            vy = self.V.differentiate('y')
            vz = self.V.differentiate('z')
            wy = self.W.differentiate('y')
        if 'Wx' in self.var: self.var['Wx'] += dt*self.dsqr(wy,vz)
        if 'Wy' in self.var: self.var['Wy'] += dt*self.dsqr(uz,wx)
        if 'Wz' in self.var: self.var['Wy'] += dt*self.dsqr(vx,uy)
        if self.dim==3: x=2*ux**2+2*uy**2+2*wz**2+(uy+vx)**2+(uz+wx)**2+(vz+wy)**2
        else:           x=2*ux**2+2*wz**2+(uz+wx)**2
        x=x.evaluate()
        if self.dim==3: x=x.integrate('y')
        x=np.reshape(x['g'][self.aslices],self.lsize)
        self.var['S']  += dt*x

        if self.solver.sim_time >= (1+self.write_num)*self.jdt:
            self.save()
            self.clear()
            

    def save(self):

        self.t2=copy.copy(self.solver.sim_time)
        self.i2=copy.copy(self.solver.iteration)

        fn = Path(str(self.dir) + '/javrg-{:05d}.hdf5'.format(self.write_num))
        ju.logger.info('Writing javrg to  {}'.format(fn))
        #ju.logger.info('self.nx  {}'.format(self.nx))
        #ju.logger.info('self.nz  {}'.format(self.nz))
        #ju.logger.info('self.x  {}'.format(self.x))
        #ju.logger.info('self.z  {}'.format(self.z))

        with h5py.File(fn, 'w', driver='mpio', comm=ju.comm) as f:
            dx  = f.create_dataset('x', (self.nx,), dtype=np.float)
            dz  = f.create_dataset('z', (self.nz,), dtype=np.float)
            dwn = f.create_dataset('write_num', (), dtype=np.int)
            dt1 = f.create_dataset('t1',        (), dtype=np.float)
            dt2 = f.create_dataset('t2',        (), dtype=np.float)
            ddt = f.create_dataset('dt',        (), dtype=np.float)
            di1 = f.create_dataset('i1',        (), dtype=np.int)
            di2 = f.create_dataset('i2',        (), dtype=np.int)
            dn  = f.create_dataset('n',         (), dtype=np.int)
            #if ju.mpirank==0:
            dwn[()] = self.write_num
            dt1[()] = self.t1
            dt2[()] = self.t2
            ddt[()] = self.dt
            di1[()] = self.i1
            di2[()] = self.i2
            dn[()]  = self.n
            dx[self.xslice]=self.x[:]
            dz[self.zslice]=self.z[:]

           
            for a in self.nm:
                d = f.create_dataset(a, self.size, dtype=np.float)
                d[self.slices]=self.var[a]
        self.write_num+=1
        


    
  
