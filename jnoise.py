import numpy as np 
import jutil as ju
import h5py
import math
import copy
from dedalus import public as de

class jnoise:
    def __init__(self,d,k,lmbd,A,domain):

        self.lmbd=lmbd     # Decay time for evolving noise field
        self.d=d           # type of noise field e.g. xyz
        self.k=k           # decay constants for spatial shaping
        self.A=A           # amplitude of noise field 

        self.cslice  = domain.dist.coeff_layout.slices(scales=1)
        self.lcshape = domain.dist.coeff_layout.local_shape(scales=1)
        self.gcshape = domain.dist.coeff_layout.global_shape(scales=1)
        self.gslice  = domain.dist.grid_layout.slices(scales=1)
        self.lgshape = domain.dist.grid_layout.local_shape(scales=1)
        self.ggshape = domain.dist.grid_layout.global_shape(scales=1)
        self.f       = domain.new_field()
        self.dim     = domain.dim
    
        self.L  = np.ndarray((self.dim,),dtype=np.float64);
        self.N  = np.ndarray((self.dim,),dtype=np.int64);
        self.CN = np.ndarray((self.dim,),dtype=np.int64);
        self.w=list()
        self.df=list()
        self.gtype=list()

        for j in range(self.dim):
            basis=domain.bases[j]
            self.gtype.append(ju.grid_type(basis.grid))
            self.L[j]=basis.interval[1]-basis.interval[0]
            if   self.gtype[j]=='Fourier':   self.w.append(basis.wavenumbers[self.cslice[j]])
            elif self.gtype[j]=='SinCos':    self.w.append(basis.wavenumbers[self.cslice[j]])
            elif self.gtype[j]=='Chebyshev': self.w.append(basis.elements[self.cslice[j]]*math.pi/self.L[j])
            else: ju.logger.info('Unknown basis type {} '.format(self.gtype[j]))
            
            self.N[j]  = basis.base_grid_size
            self.CN[j] = self.w[j].size

        if self.dim==2:
            self.w[0] = np.reshape(self.w[0],(self.CN[0],1))
            self.w[1] = np.reshape(self.w[1],(1,self.CN[1]))
        else:
            self.w[0] = np.reshape(self.w[0],(self.CN[0],1,1))
            self.w[1] = np.reshape(self.w[1],(1,self.CN[1],1))
            self.w[2] = np.reshape(self.w[2],(1,1,self.CN[2]))

        for j in range(self.dim):
            if k[j]>0: x=k[j]*np.exp(-k[j]*np.abs(self.w[j]))
                #self.df[j][self.w[j]==0]=0
                #
            else: x=np.ones(self.w[j].shape,dtype=np.float64)
            x[self.w[j]==0]=0    # knock out constant term
            sx=np.sqrt(ju.jsum(x**2))
            self.df.append(x/sx) # maybe should be /2 for fourier


        if self.dim==2:
            if self.gtype[0]=='Fourier' and self.gtype[1]=='Fourier': self.df[0]=self.df[0]/np.sqrt(2)
            if self.gtype[0]=='SinCos'  and self.gtype[1]!='Fourier': self.df[0]=self.df[0]*2
        else:
            if self.gtype[0]=='Fourier' and self.gtype[1]=='Fourier' and self.gtype[2]=='Fourier': self.df[0]=self.df[0]/np.sqrt(2)
            if self.gtype[0]=='SinCos'  and self.gtype[1]=='Fourier' and self.gtype[2]!='Fourier': self.df[0]=self.df[0]*np.sqrt(2)
            if self.gtype[0]=='Fourier' and self.gtype[1]=='SinCos'  and self.gtype[2]!='Fourier': self.df[0]=self.df[0]*np.sqrt(2)
            if self.gtype[0]=='SinCos'  and self.gtype[1]=='SinCos'  and self.gtype[2]=='Fourier': self.df[0]=self.df[0]*np.sqrt(2)
            if self.gtype[0]=='SinCos'  and self.gtype[1]=='SinCos'  and self.gtype[2]!='Fourier': self.df[0]=self.df[0]*np.sqrt(2)*2
                     
        if not 'x' in d:
            self.lcshape[0]=1
            self.cslice[0]=slice(0,1)
            self.f.meta['x']['constant']=True
            self.df[0]=1
            self.N[0] =1
            self.gshape[0] = 1
            self.lgshape[0] = 1
        else:
            if self.gtype[0]=='SinCos': self.f.meta['x']['parity']=1
          
        if self.dim==3:
            if not 'y' in d:
                self.lcshape[1]=1
                self.cslice[1]=slice(0,1)
                self.f.meta['y']['constant']=True
                self.df[1]=1
                self.gshape[1] = 1
                self.lgshape[1] = 1
            else:
                if self.gtype[1]=='SinCos': self.f.meta['y']['parity']=1
             
        if not 'z' in d:
            self.lcshape[-1]=1
            self.cslice[-1]=slice(0,1)
            self.f.meta['z']['constant']=True
            self.df[-1]=1
            self.gshape[-1] = 1
            self.lgshape[-1] = 1
        else:
            if self.gtype[-1]=='SinCos': self.f.meta['z']['parity']=1
 
        self.ggsize = self.ggshape.prod()
        self.lgsize = self.lgshape.prod()
        self.gcsize = self.gcshape.prod()
        self.lcsize = self.lcshape.prod()

        #ju.print_vars(self)
        
        if self.dim==3: self.f['c'] = self.df[0]*self.df[1]*self.df[2]
        else:           self.f['c'] = self.df[0]*self.df[1]          

        self.irms = self.rms()

        self.f.set_scales(1)
        self.f['c'] = self.randc()
        self.normalize()
        for j in range(self.dim): ju.logger.info('df[{}].size: {}'.format(j,self.df[j].shape))
        #for j in range(self.dim): ju.logger.info('w[{}].shape: {}'.format(j,self.w[j].shape))
        ju.logger.info('self.lcshape: {}'.format(self.lcshape))
        

        
        
        #ju.logger.info('uE {:9.2e}, dnoise {:8.4f}, NYZ {:8.4f}'.format(uE,dnoise,np.sqrt(NY**2+NZ**2)))
        #ju.logger.info('uE {:9.2e}, EQ {:8.4f}, s+b {:8.4f}, Fx {:8.4f}, Ifd {:9.2e}, NN {:8.4f}'.format(uE,EQ/(U*IA*H*W),Isb,state_ev['Fx'],ju.jint(fd),NN))


        #JA.savexyz('fu/fu-{:05d}.hdf5'.format(solver.iteration),'u',fu['g'],AA)
        #JA.savexyz('fd/fd-{:05d}.hdf5'.format(solver.iteration),'d',fd['g'],AA)

    def rms(self):
        self.f.set_scales(1)
        r=np.sqrt(ju.jsum(self.f['g']**2)/self.ggsize)
        return(r)

    def normalize(self):
        self.f['g']=copy.copy(self.f['g']*self.A*np.sqrt(ju.jsum(self.f['g']**2)/self.ggsize))
        
    def print(self):
        vv=list(set(vars(self)) - {'w','df','f'}) 
        for v in vv: ju.logger.info('{}: {}'.format(v,eval('self.'+v)))
        for j in range(self.dim):
            ju.logger.info('w[{}].shape:   {} [{:7.2f} {:7.2f}]'.format(j, self.w[j].shape,   ju.jmin(self.w[j]  ), ju.jmax(self.w[j]  )))
            ju.logger.info('df[{}].shape:  {} [{:7.2f} {:7.2f}]'.format(j, self.df[j].shape,  ju.jmin(self.df[j] ), ju.jmax(self.df[j] )))
        ju.logger.info('f[g].shape: {} [{:7.2f} {:7.2f}]'.format(        self.f['g'].shape, ju.jmin(self.f['g']), ju.jmax(self.f['g'])))
        ju.logger.info('rms: {}'.format(self.rms()))


    def save(self,fnm,JA):
        f.set_scales(1)
        JA.savexyz(fnm,'f',f['g'],1)

    def randc(self):
        if self.dim==3: c = self.A*self.df[0]*self.df[1]*self.df[2]*np.random.normal(size=self.lcshape)
        else:           c = self.A*self.df[0]*self.df[1]*           np.random.normal(size=self.lcshape)
        return(c)
    
    def int(self,dt):
        (a,b)=ju.noisedt(dt,self.lmbd,1)
        self.f.set_scales(1)
        self.f['c'] = a*self.f['c']+b*self.randc()
        #ju.logger.info('jnoise.int: rms {:7.2f}'.format(self.rms()))
    
        
