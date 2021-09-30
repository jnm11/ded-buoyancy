import numpy as np
import jutil as ju
import copy as copy

global nmfloat
global nmbool 
global nmint
global nmstr


nmfloat=('Xmin','Xmax','FB','hab','meanu','Length','Width','Height','Velocity','Radius','Time','Angle','Gravity','Buoyancy','Sediment','dimple','dimplewy','dimplewz','ddT','gmin','gmax','Umin','Umax','mindt','maxdt','Force','forceb','forces','clipu','Re','Scb','Scs','ddSc','Peb','Pes','Ped','AA','AAS','AAJ','VMEM','wu','wb','hu','hb','U1','Qu','SV','bV','SVh','bVh','Noise','noiseT','noiseL','xa','xb','r0','x0','x1','x2','x3','x4','x5','x6','x7','xn1','xn2','xn3','xn4','Wx','Wy','Wz','Wr','dtcheck','dtjcheck','dtb','dtu','dtv','dtw','dtp','dtnoise','dtstats','dtx','dty','dtz','dtxy','dtxz','dtyz','dtxyz','dtmomb','dtavrg','dtleft','dtright','dtpm','dtgc','dtfpid','dtforce','dtslice','fpid','PIDL','PIDD','PIDI','PIDP','PIDT','PIDX','PIDIT','PIDDD','PIDST','PIDS1','PIDS2','inoise','bnoise','snoise','slicex','slicey','slicez','dtj','dtjx','dtjy','dtjz','dtjr','dtjyz','dtjxz','dtjxyz','dtja','dtjar','dtjax','dtjay','dtjaz','dtjaxy','dtjayz','dtjaxz','dtjaxyz','dtjb','dtjd','dtjs','dtjp','dtju','dtjv','dtjw','dtjWx','dtjWy','dtjWz','dtjS','dtjE','dtjaS','dtjaE','dtjsev','dtjavar','dtjdivyz','Fu','Fv','Fw','dnoise','minB','maxB','pmss','pmsl','Skp','Ski','Skg','sinitial','fuspecial1','fuspecial2','topT','wdivxl','wdivxr','topddM','dbT','dbws','ddws','ddlogT','dblogT')

nmbool=('cm','PIDG','signals','parallel','avrg','fpidu','fpidv','fpidw','fpidb','fpids','clipB','clipS','fbmult','fsmult','fbmin','fsmin','fbmax','fsmax','Conservative','db','oddg','dtjadv','dtjsw','dtjd1','dtjd2','pmzr','pmzt','pmzrt','topdrt','topd','topr','topdd','topu','ful','fur','wul','wvl','wur','wvr','wwl','wwr','rescale','fixI','pmIu','dd','ddserf','ddq','ddscl','intdx','divxI','topq','topdivr','dddd','dbu','dbz','dbserf','dbdirect','ddserf','ddlogh','dblogh','ddD','dbD','dbin')   

nmint=('Nx','Ny','Nz','Forcing','m1','m2','ddj','ddk','dbj','dbk')

nmstr=('noised','scheme','Inlet','Tx','Ty','Tz','lbc','rbc','slbc','srbc','blbc','brbc','bntype','fuspecial','extx','exty','extz','WallTime','sType')  


#rm -rf ~/pm/qpmf0 ; mpiexec -n 10 ded_gc-1.31.py --preset qpmf0 pm/qpmf0
#rm -rf ~/gc/ccle1 ; mpiexec -n 10 ded_gc-1.31.py --preset cccle(pcle1 gc/ccle1
#rm -rf ~/gc/qgcf3 ; mpiexec -n 10 ded_gc-1.31.py --preset qgcf3 gc/qgcf3
#rm -rf ~/gc/qgcf5 ; mpiexec -n 4  ded_gc-1.31.py --preset qgcf5 gc/qgcf5
#rm -rf ~/gc/qgcf6 ; mpiexec -n 10 ded_gc-1.31.py --preset qgcf6 gc/qgcf6
#rm -rf ~/gc/qlr1  ; mpiexec -n 4 ded_gc-1.31.py --preset qlr1 gc/qlr1
def merge_param(p,q):
    for a in q: p[a]=q[a]
    return(p)

def preset(p,s):
    ju.logger.info('Loading preset {}'.format(s))
    if   s == 'pmf0':   p=pmf0(p)  # plume forcing 0
    elif s == 'pmf1':   p=pmf1(p)  # plume forcing 1
    elif s == 'pmf2':   p=pmf2(p)  # plume forcing 2
    elif s == 'pmf3':   p=pmf3(p)  # plume forcing 3
    elif s == 'pmf4':   p=pmf4(p)  # plume forcing 4
    elif s == 'pmf6':   p=pmf6(p)  # plume forcing 6
    elif s == 'pmf7':   p=pmf7(p)  # plume forcing 7
    elif s == 'pmf8':   p=pmf8(p)  # plume forcing 8
    elif s == 'gcf0':   p=gcf0(p)  # gravity gravity no forcing
    elif s == 'gcf1':   p=gcf1(p)  # gravity gravity forcing 1
    elif s == 'gcf2':   p=gcf2(p)  # gravity gravity forcing 2
    elif s == 'gcf3':   p=gcf3(p)  # gravity gravity forcing 3
    elif s == 'gcf4':   p=gcf4(p)  # gravity gravity forcing 4             
    elif s == 'gcf5':   p=gcf5(p)  # gravity gravity forcing 5           
    elif s == 'gcf6':   p=gcf6(p)  # gravity gravity forcing 6           
    elif s == 'gcf7':   p=gcf7(p)  # gravity gravity forcing 7           
    elif s == 'le':     p=le(p)    # Lock exchange gravity current
    elif s == 'ccle':   p=ccle(p)  # CC Lock Exchange gravity current
    else:
        ju.logger.error("preset: unknown {}".format(s))
        ju.abort_run()

    return(p)

# Parameters do not need to be specified
def param_default():    
    # Coordinate types
    p={}
    p['scheme']    = np.str('RK443')
    p['WallTime']      = np.str('99:00:00')
    p['Tx']        = np.str('Fourier')
    p['Ty']        = np.str('Fourier')
    p['Tz']        = np.str('Fourier')
    p['maxdt']     = 5e-2
    p['mindt']     = 1e-4
    p['Force']     = 100
    p['forceb']    = 1
    p['forces']    = 1
    p['Forcing']   = 0
    p['Height']=None
    p['Length']=None
    p['Width']=None

    p['Nx']=None
    p['Ny']=None
    p['Nz']=None
    
    p['Wx']=None
    p['Wy']=None
    p['Wz']=None

    # Record parameters
    p['dtstats']   = 0.1  
     
    # Noise forcing paramters
    #p['noised']  = np.str('xz') # Noise direction
    #p['noiseL']  = 1    # Laplacian soothing, 1/L^2
    #p['noiseT']  = 1    # Inverse correlation time 1/T

    p['Angle']         = 0
    
    p['U']         = 1

    # General parameters    
    p['Scb']       = 1
    p['Scs']       = 1
    p['Gravity']         = 1
    p['Time']         = 1000
    p['AA']        = 1.5
    p['AAJ']       = 1
    p['AAS']       = 1

    p['Buoyancy']               = 1
    p['brbc']            = 'dz'
    p['blbc']            = 'dz'
    p['srbc']            = 'dz'
    p['slbc']            = 'dz'
    
    # Program control
    p['dtjcheck']        = 30
    p['signals']         = ju.bool(True)
    p['parallel']        = ju.bool(True)

    p['hu']      = 1
    p['hb']      = 1

    p['pmsl']    = 0.1
    return(p)

def set_PID(p):
    #PID controller for gravity current front
    p['PIDL']      = 1
    p['PIDP']      = 0.5          
    p['PIDI']      = 0.1
    p['PIDG']      = ju.bool(True)
    p['PIDD']      = 0.05
    p['PIDT']      = 0.01
    p['PIDX']      = 20    
    p['PIDS1']     = 5    
    p['PIDS2']     = 7    
    p['gmin']      = 0.05
    p['gmax']      = 30
    p['Umin']      = 0.05
    p['Umax']      = 2
    return(p)

def gc(p): 
    p['lbc']     = np.str('slip')
    p['rbc']     = np.str('slip')
    p['Tz']      = np.str('Cheb')
    p['Angle']       = 0
    p['Re']      = 2000 
    p['Time']       = 100 
    p['Length']       = 24
    p['Width']       = 2
    p['Height']       = 2
    p['Gravity']       = 1.5
    #p['noised']  = np.str('xz')
    #p['noiseL']  = 1
    #p['noiseT']  = 1
    #p['xn1']     = 0
    #p['xn2']     = 0
    #p['xn3']     = 0
    #p['xn4']     = 0 
    return(p)
    
def pm(p):
   #plume parameters
    p['Inlet']   = np.str('gaussian')
    p['Re']      = 200
    p['Radius']  = 0.5   # radius of inlet
    p['Buoyancy']       = 1     # inlet velocity
    p['Time']       = 1000
    p['Angle']       = 90
    p['Tx']      = np.str('Fourier')
    p['Ty']      = np.str('Fourier')
    p['Tz']      = np.str('Fourier')
    p['Length']       = 16
    p['Width']       = 4
    p['Height']       = 4
    p['noised']  = 'yz'
    p['dtjb']    = 1     
    p['dtjyz']   = 0.1     
    p['dtstats'] = 0.1   
    p['topu']    = 1
    p['topddM']  = 1
    p['wdivxl']  = 4
    p['wdivxr']  = 1
    p['topT']    = 1
    p['ddSc']     = 1
  
    return(p)

def pmf0(p):
    p=pm(p)
    p['Forcing'] = 0
    p['x0']      = 1.0 # [ 0 x0] negative div(u)
    p['x1']      = 1.5 # [x0 x1] zero velocity
    p['x2']      = 2.5 # [x1 x2] positive div(u)
    p['x3']      = 3.0 # [x2 x3] constant u
    return(p)

def pmf1(p): # Set all parameters for a quick plume simulation with forcing 1
    p=pm(p)
    p['Forcing'] = 1
    p['Angle']       = -180
    p['Tz']      = np.str('Cheb')
    p['Nx']      = 64
    p['Ny']      = 64
    p['Nz']      = 256
    p['Length']       = 4
    p['Width']       = 4
    p['Height']       = 16
    p['x0']      = 0.5 # [ 0 x0] negative div(u)
    p['x1']      = 1.0 # [x0 x1] zero velocity
    p['x2']      = 1.5 # [x1 x2] positive div(u_
    p['x3']      = 2.0 # [x2 x3] constant u
    return(p)

def pmf2(p):
    p=pm(p)
    p['Forcing'] = 2
    p['x0']      = 0.5 # 0.5 [0.0 0.5] gaussian 1    
    p['x1']      = 1.0 # 1.0 [0.5 1.0] switch region 
    p['x2']      = 1.5 # 1.5 [1.0 1.5] switch region 
    p['x3']      = 2.0 # 2.0 [1.5 2.0] gaussian 2    
    return(p)

def pmf3(p):
    p=pm(p)
    p['Forcing'] = 3
    p['x0']      = 0.50
    p['x1']      = 0.50
    p['x2']      = 1.75
    p['x3']      = 2.00
    p['x4']      = 0.00
    p['x5']      = 1.00
    p['x6']      = 1.50
    p['x7']      = 2.50
    p['xn1']     = 1.50 # Noise forcing region
    p['xn2']     = 2.0
    p['xn3']     = 2.0
    p['xn4']     = 2.5
    return(p)

def pmf4(p):
    p=pm(p)
    p['Forcing'] = 4
    p['x0']      = 1.0
    p['x1']      = 1.0
    p['x2']      = 2.5
    p['x3']      = 3.0
    p['x4']      = 0.0
    p['x5']      = 1.0
    p['x6']      = 2.0
    p['x7']      = 3.0
    p['xn1']     = 2.5 # Noise forcing region
    p['xn2']     = 3.0
    p['xn3']     = 3.0
    p['xn4']     = 3.5
    return(p)

def pmf6(p):
    p=pm(p)
    p['Forcing'] = 6
    p['x0']      = 1.0
    p['x1']      = 1.5
    p['x2']      = 2.0
    p['x3']      = 3.0
    p['x4']      = 0.0
    p['x5']      = 1.0
    p['x6']      = 2.0
    p['x7']      = 3.0
    p['Length']       = 30
    p['Width']       = 20
    p['Height']       = 20
    p['Inlet']   = np.str('circle')
    p['Radius']  = 1     # radius of inlet
    p['pmss']    = 0.1
    p['pmsl']    = 5
    p['Wx']      = 0.1
    p['Wy']      = 0.1
    p['Wz']      = 0.1
    p['Re']      = 500
    p['dtjb']    = 0.1
    p['dtjayz']  = 1
    p['dtja']    = 20
    p['Buoyancy']       = 1
    return(p)

def pmf7(p):
    p=pm(p)
    p['Forcing']  = 7
    p['Tx']       = 'SinCos'
    p['Ty']       = 'Fourier'
    p['Tz']       = 'Fourier'
    p['Length']        = 30
    p['Width']        = 20
    p['Height']        = 20
    p['Inlet']    = np.str('circle')
    p['Radius']   = 0.5 # radius of inlet
    p['pmss']     = 0.1
    p['pmsl']     = 5
    p['Re']       = 500
    p['dtjb']     = 0.1
    p['dtjayz']   = 1
    p['dtja']     = 20
    p['dtjsev']   = 1
    p['Buoyancy']        = 1
    p['r0']       = 0.5
    p['ddT']      = 1
    p['dbT']      = 1
    p['dd']     = True
    p['ddu']    = False
    p['Wr']       = None
    p['hu']       = None
    p['hb']       = None
    p['x0']       = None 
    p['x1']       = None 
    p['x2']       = None 
    p['x3']       = None 
    p['x4']       = None
    p['x5']       = None 
    p['x6']       = None 
    p['x7']       = None    
    p['Gravity']        = 1
    p['pmzr']     = None
    p['pmzt']     = True
    p['pmzrt']    = None
    p['wul']      = None
    p['wur']      = None
    p['wvl']      = True
    p['wur']      = None
    p['wwl']      = None
    p['wwr']      = None
    p['topd']     = None
    p['topr']     = None
    p['topu']     = True
    p['oddg']     = False
    p['ddscl']  = None
    p['wdivxl']   = 4
    p['wdivxr']   = 1
    return(p)

def pmf8(p):
    p=pm(p)
    p['Forcing']  = 8
    p['Tx']       = 'SinCos'
    p['Ty']       = 'Fourier'
    p['Tz']       = 'Fourier'
    p['Length']   = 30
    p['Width']    = 20
    p['Height']   = 20
    p['Inlet']    = np.str('circle')
    p['Radius']   = 0.25 # radius of inlet
    p['pmss']     = 0.1
    p['pmsl']     = 5
    p['Re']       = 200
    p['dtjb']     = 0.1
    p['dtjayz']   = 1
    p['dtja']     = 20
    p['dtjsev']   = 1
    p['Buoyancy'] = 1
    p['r0']       = None
    p['ddT']      = None
    p['dbT']      = None
    p['dd']       = None
    p['Wr']       = None
    p['hu']       = None
    p['hb']       = None
    p['x0']       = None 
    p['x1']       = None 
    p['x2']       = None 
    p['x3']       = None 
    p['x4']       = None
    p['x5']       = None 
    p['x6']       = None 
    p['x7']       = None    
    p['Gravity']        = 1
    p['pmzr']     = None
    p['pmzt']     = True
    p['pmzrt']    = None
    p['wul']      = None
    p['wur']      = None
    p['wvl']      = True
    p['wur']      = None
    p['wwl']      = None
    p['wwr']      = None
    p['topd']     = None
    p['topr']     = None
    p['topu']     = None
    p['oddg']     = None
    p['topT']     = None
    p['ddSc']      = None
    p['ddscl']  = None
    p['wdivxl']   = None
    p['wdivxr']   = None
    return(p)
       
#round(cumprod([1 1.2500    1.2500    1.2800    1.2500    1.2500    1.2800    1.2500    1.2500    1.2800    1.2500    1.2500    1.2650    1.2648    1.2500    1.2500    1.2800    1.2500    1.2500    1.2800    1.2500])*256*3/2/2)*2
#round(cumprod([1 1.2500    1.2500    1.2800    1.2500    1.2500    1.2800    1.2500    1.2500    1.2800    1.2500    1.2500    1.2650    1.2648    1.2500    1.2500    1.2800    1.2500    1.2500    1.2800    1.2500])*256*2)*2
#def setN(p,n,m):
#    if m==1:  Nx=np.array([384,480,600,768,960,1200,1536,1920,2400,3072,3840,4800,6072,7680,9600,12000,15360,19200,24000,30720,38400],np.int64)
#    else:     Nx=np.array([256,320,400,512,640, 800,1024,1280,1600,2048,2560,3200,4048,5120,6400, 8000,10240,12800,16000,20480,25600],np.int64))
#    p['Nx'] = Nx[n-1]
#    p['Ny'] = 2*np.int64(  p['Nx']*p['Width']/p['Length']/2)  
#    p['Nx'] = 2*np.int64(3*p['Nx']*p['Height']/p['Length']/4)  
#    p['Wx'] = p['Length']/p['Nx']*6
#    p['Wy'] = p['Width']/p['Ny']*6
#    p['Wz'] = p['Height']/p['Nz']*6
 
def le(p): # Lock exchange gravity current
    p['Forcing'] = 0
    p['Tx']      = np.str('SinCos')
    p['Ty']      = np.str('Fourier')
    p['Tz']      = np.str('Cheb')
    p['Angle']       = 0
    p['U']       = None
    p['U1']      = None
    p['inoise']  = 0.1
    p['Gravity']       = 1
    p['Length']       = 40 
    p['Width']       = 2
    p['Height']       = 2
    p['hu']      = 2
    p['hb']      = 2
    p['xb']      = 10 
    p['lbc']     = np.str('noslip')
    p['rbc']     = np.str('noslip')
    p['Time']       = 50
    p['Re']      = 200
    return(p)

def ccle(p):     # CC experiment
    p=le(p)
    p['Re']      = 5000
    p['Length']       = 30 
    p['Width']       = 2.5
    p['Height']       = 1
    p['hu']      = 1
    p['hb']      = 1
    p['Time']       = 50
    p['xb']      = 10
    p['lbc']     = np.str('noslip')
    p['rbc']     = np.str('slip')
    return(p)

def gcf3(p):     # forcing 3
    p=gc(p)
    p['Forcing'] = 3
    p['x0']      = 0.0
    p['x1']      = 0.5
    p['x2']      = 1.0
    p['x3']      = 2.0
    p['x4']      = 0.5 # start of buoyancy increase
    p['x5']      = 1.0 # end   of buoyancy increase
    p['x6']      = 2.0 # start of buoyancy forcing decrease
    p['x7']      = 3.0 # end   of buoyancy forcing decrease   
    p['U1']      = 0 
    return(p)

def gcf5(p): #  5
    p=gc(p)
    p['Forcing'] = 5
    p['x0']      = 0.0 # Region were velocity is forced to -U
    p['x1']      = 0.5
    p['x2']      = 1.0
    p['x3']      = 1.5
    p['x4']      = 1.0 # start of buoyancy increase
    p['x5']      = 1.5 # end   of buoyancy increase
    p['x6']      = 2.0 # start of buoyancy forcing decrease
    p['x7']      = 2.5 # end   of buoyancy forcing decrease   
    return(p)

def gcf6(p): # Set gravity current forcing to 6 
    p=gc(p)
    p['Forcing'] = 6
    # x0-x1 is constant velocity x1 to x3 is streamline adjustmenat  x2 to x3 ia weight decrease
    # Buoyancy is x4 to x7
    p['x0'] = 0.00
    p['x1'] = 0.25
    p['x2'] = 1.00
    p['x3'] = 1.25
    p['x4'] = 0.50
    p['x5'] = 1.00
    p['x6'] = 2.00
    p['x7'] = 2.50    
    return(p)

def gcf7(p): # Set gravity current forcing to 7 
    p=gc(p)
    p['Forcing'] = 7
    p['Tx']      = 'SinCos'
    p['Length']  = 20
    p['Width']   = 1
    p['Height']  = 2
    p['hu']      = 0.8
    p['Gravity'] =  1
    p['dtja']    = 20
    p['dtjayz']  = 1
    p['dtjay']   = 10
    p['dtjb']    = 1
    p['dtjsev']  = 1
    p['dtjavar'] = 1
    p['Force']   = 50
    p['maxdt']   = 0.1
    p['ddT']     = 1
    p['dbT']     = 1
    p['dd']      = True
    p['db']      = True
    p['ddj']     = 10
    p['ddk']     = 2
    p['dbj']     = 10
    p['dbk']     = 2
    p['ful']     = True
    p['fur']     = True
    return(p)

def gcf63(p): # Set gravity current forcing to 6 but make it behave like forcing 3 
    p=gc(p)
    p['Forcing'] = 6
    # x0 to x3 control velocity
    # Buoyancy is x4 to x7
    # x0-x1 is constant velocity x1 to x3 is streamline adjustment  x2 to x3 is weight decrease
    p['x0'] = 0.00
    p['x1'] = 0.50
    p['x2'] = 1.50
    p['x3'] = 2.50
    p['x4'] = 0.50
    p['x5'] = 1.00
    p['x6'] = 2.00
    p['x7'] = 2.50    
    return(p)

def qgcf5st(p): # forcing 5 with step in PIDX
    p=qgcf5(p)
    p['PIDST']   = 25
    p['PIDS1']   = 5
    p['PIDS2']   = 7
    p['PIDP']    = 1.0          
    p['PIDI']    = 0.5
    p['PIDD']    = 0.1
    p['U']       = 0.3
    return(p)

def name_type(args,restart):
    typs=('gc','bg','pm')
    name=''
    p={}
    if ju.mpirank==0:
        name =  np.str(args['NAME'])
        if not restart: p=param_default()
        if len(name)>2:
            for a in typs:
                if name[0:3]==a +'/':  p['sType'] = np.str(a)
            name=name[3:]
        try:
            ju.logger.info('Type: {}'.format(p['sType']))
            ju.logger.info('Name: {}'.format(name))
        except:
            ju.logger.info('name_type: Unknown simulation type "{}"'.format(name))
            ju.set_error('type')
    ju.check_abort()
    p=ju.comm.bcast(p, root=0)
    name=ju.comm.bcast(name,root=0)
    return(p,name)


def jint(s):
    if s=='None': r=None
    else:         r=np.int64(s)
    return(r)

def jfloat(s):
    if s=='None': r=None
    else:         r=np.float64(s)
    return(r)

def jstring(s):
    if s=='None': r=None
    else:         r=np.str(s)
    return(r)

def jbool(a):
    if   a=="None":  x=None
    elif a=="False": x=np.bool(False)
    elif a=="0":     x=np.bool(False)
    elif a=="F":     x=np.bool(False)
    elif a=="f":     x=np.bool(False)
    elif a=="True":  x=np.bool(True)
    elif a=="1":     x=np.bool(True)
    elif a=="T":     x=np.bool(True)
    elif a=="f":     x=np.bool(True)
    else:
        logger.error("jbool: Boolean argument must be True, False, 0 or 1")
        set_abort('jbool')
    return(x)
           
def read_args(p,args):
    global nmfloat
    global nmbool
    global nmint
    global nmstr

    for a in nmfloat:
        try:
            if args['--' + a]: p[a] =  jfloat(args['--' + a])
        except:
            ju.logger.info('{}  should be a float not "{}"'.format(a,args['--' + a]))

    for a in nmbool:
        try:
            if args['--' + a]: p[a] =  jbool(args['--' + a])
        except:
            ju.logger.info('{}  should be a boolean not "{}"'.format(a,args['--' + a]))

    for a in nmint:
        if args['--' + a]:
            try: p[a] =  jint(args['--' + a])
            except: ju.logger.info('{} should be an int not "{}"'.format(a,args['--' + a]))

    for a in nmstr:
        try:
            if args['--' + a]: p[a] =  jstring(args['--' + a])
        except:
            ju.logger.info('{}  should be an string not "{}"'.format(a,args['--' + a]))
            
    #ju.logger.info('read_args: p={}'.format(p))
    return(p)

def tobool(a):
    if   a is None:  return(None)
    elif a==True:    return(True)
    elif a==False:   return(False)
    elif a==1:       return(True)
    elif a==0:       return(False)
    elif a=="True":  return(True)
    elif a=="False": return(False)
    elif a=="T":     return(True)
    elif a=="F":     return(False)
    else:            return(None)

def remove_None(p):
    b=copy.copy(p)
    for a in p:
        if b[a] is None: del(b[a])
    return(b)

def remove_zeros(p,v):
    b=copy.copy(p)
    for a in v:
        if a in p:
            if b[a]==0: del(b[a])
    return(b)

def checkp(p,a,b,s):
    if a in p:
        if p[a] is not None:
            if p[a] not in b:
                ju.logger.error('validate: ' + s + ' {}'.format(p[a]))
                ju.abort_run()


def jpop(p,a):
    for b in a:
        if b in p: p.pop(b)

def write_options(p,fn):
    global nmfloat
    global nmbool
    global nmint
    global nmstr
    nmstr   = list(set(nmstr)   & set(p.keys()))
    nmint   = list(set(nmint)   & set(p.keys()))
    nmbool  = list(set(nmbool)  & set(p.keys()))
    nmfloat = list(set(nmfloat) & set(p.keys()))

    if ju.mpirank==0:
        f=fn.open('w')
        for a in nmstr:   f.write('--{} {} '.format(a,p[a]))
        for a in nmint:   f.write('--{} {} '.format(a,p[a]))
        for a in nmbool:  f.write('--{} {} '.format(a,p[a]))
        for a in nmfloat: f.write('--{} {} '.format(a,p[a]))
        f.write('\n')
        f.close()

def validate(p):
    global nmfloat
    global nmbool
    global nmint
    global nmstr

    # Old variables to delete
    jpop(p,('f7s','wdiffdom','ddeven','AB','alpha','adddudx','dddf','gclf','gcrf','topdudx','ddU','dU1','dUX'))
    
    # Old variables to rename
    b = (
        ('dddd',     'diffdom'),
        ('AAJ',      'AAj'),
        ('Noise',    'noise'),
        ('Force',    'force'),
        ('Forcing',  'forcing'),
        ('Inlet',    'inlet'),
        ('Radius',   'radius'),
        ('Length',   'length'),
        ('Width',    'width'),
        ('Height',   'height'),
        ('Radius',   'radius'),
        ('Angle',    'angle'),
        ('WallTime', 'time'),
        ('Buoyancy', 'B'),
        ('Sediment', 'S'),
        ('Time',     'T'),
        ('Velocity', 'U'),
        ('Radius',   'R'),
        ('Length',   'L'),
        ('Width',    'W'),
        ('Height',   'H'),
        ('Radius',   'Radius'),
        ('Angle',    'q'),
        ('Gravity',  'g'),
        ('dtjadv',   'dtjudel'),
        ('dtjsw,'    'dtjdiff'),
        ('ddws',     'ddivws'),
        ('dd',       'ddiv'),
        ('ddserf',   'ddivserf'),
        ('ddq',      'ddivq'),
        ('ddscl',    'ddivscl'),
        ('dddd',     'ddddiv'),
        ('ddserf',   'ddivserf'),
        ('ddj',      'ddivj'),
        ('ddk',      'ddivk'),     
        ('ddT',      'dT'),     
         )

    for a in b:
        if a[1] in p:
            if a[0] not in p: p[a[0]]      = p.pop(a[1])
            else:             p.pop(a[1])
            
    
    if 'bminmax' in p:
        p['fbmax']=p.pop('bminmax')
        p['fbmin']=p['fbmax']
    if 'sminmax' in p:
        p['fsmax']=p.pop('sminmax')
        p['fsmin']=p['fsmax']
        
    if 'ww' in p: p['wwl']=p.pop('ww')
    p=remove_zeros(p,('Forcing','fuspecial'))
                   
    #if 'dd' in p: p['dd']=np.bool(p['dd'])
    
    # Deal with old parameters
    if 'ck_time' in p:
        p['dtcheck']=p['ck_time']
        p.pop('ck_time')

    if 'Sc' in p:
        p['Scb']=p['Sc']
        p['Scs']=p['Sc']
        p.pop('Sc')

    if 'dT' in p:
        p['dbT']=p['dT']
        p['ddT']=p['dT']
        p.pop('dT')
        
    if 'ubc' in p:
        p['rbc']=p['ubc']
        p.pop('ubc')

              
    # Convert old style bytes into str
    nmstr   = list(set(nmstr)   & set(p.keys()))
    for a in nmstr:
        #ju.logger.info('  {:10s}: {} {}'.format(a,type(p[a]),p[a]))
        if not isinstance(p[a],str): p[a]=p[a].decode('utf-8')
        #ju.logger.info('  {:10s}: {} {}'.format(a,type(p[a]),p[a]))

    if 'Peb' in p: p['Scb']=p.pop('Peb')/p['Re']
    if 'Pes' in p: p['Scs']=p.pop('Pes')/p['Re']
    if 'Ped' in p: p['ddSc']=p.pop('Ped')/p['Re']

    if not 'avrg' in p and 'dtavrg' in 'p': p['avrg'] = p['dtavrg']>0

    if p['mindt'] <=0:
        ju.logger.info('ded_type:validate: mindt must be positive')
        p['mindt']=1e-7

    if p['sType']=='pm':
        jpop(p,('Wx','Wy','Wz','hb','hu'))
        if 'topT' not in p: p['topT'] = 1 
        if 'r0'   not in p: p['r0']   = 0.5
    else: jpop(p,('pmsl','pmss','Radius','Inlet'))

    if p['sType']=='gc':
        if 'gmin' not in p: p['gmin'] = 0.05
        if 'gmax' not in p: p['gmax'] = 30
        if 'Umin' not in p: p['Umin'] = 0.05
        if 'Umax' not in p: p['Umax'] = 2
        if 'Buoyancy' in p and p['Tz']=='Cheb':
            if 'blbc' not in p: p['blbc'] = 'dz'
            if 'brbc' not in p: p['brbc'] = 'dz'
        if 'Sediment' in p and p['Tz']=='Cheb':
            if 'slbc' not in p: p['slbc'] = 'dz'
            if 'srbc' not in p: p['srbc'] = 'dz'
            
    else: jpop(p,('PIDP','PIDL','PIDI','PIDG','PIDD','PIDT','PIDX','PIDIT','PIDDD','PIDST','PIDS1','PIDS2'))

    fuspecial = get_param(p,'fuspecial')
    forcing   = get_param(p,'Forcing')
    
    if forcing is None and fuspecial is None: jpop(p,('Force','forceb','forces','x0','x1','x2','x3','x4','x5','x6','x7','fuspecial'))
    
     
  
    # Validate input parameters
    if p['sType']=='plume': p['sType']='pm'

    checkp(p,'lbc',('slip','noslip'),'Unknown left  boundary condition')
    checkp(p,'rbc',('slip','noslip'),'Unknown right boundary condition')
    checkp(p,'Tx',('Fourier','SinCos'),'Unknown x grid type')
    checkp(p,'Ty',('Fourier','SinCos'),'Unknown y grid type')
    checkp(p,'Tz',('Fourier','SinCos','Cheb'),'Unknown z grid type')
    checkp(p,'Inlet',('gaussian','circle'),'Unknown inlet type')
    checkp(p,'noised',('xy','xz','yz'),'Unknown noised')

    
    if p['Tz'] != 'Cheb': jpop(p,('slbc','srbc','blbc','brbc'))
        
    if p['Nx'] is not None and p['Length'] is not None: dx=p['Length']/p['Nx']
    else: dx=None
    if p['Ny'] is not None and p['Width'] is not None: dy=p['Width']/p['Ny']
    else: dy=None
    if p['Nz'] is not None and p['Height'] is not None: dz=p['Height']/p['Nz']
    else: dz=None

    if 'Tx' in p:
        if p['Tx']=='Cheb' and dx is not None: dx=dx*3/2
    if 'Ty' in p:
        if p['Ty']=='Cheb' and dy is not None: dy=dy*3/2
    if 'Tz' in p:
        if p['Tz']=='Cheb' and dz is not None: dz=dz*3/2
    
    if dx==None and dy==None and dz==None:
        ju.logger.info('At least one dimension must have a size and grid count')
        ju.logger.info('Nx {}, Ny {}, Nz, {}'.format(p['Nx'],p['Ny'],p['Nz']))
        ju.logger.info('L  {}, W  {}, H,  {}'.format(p['Length'],p['Width'],p['Height']))
        ju.logger.info('dx {}, dy {}, dz, {}'.format(dx,dy,dz))
        ju.abort_run()
 
    if p['Width'] is None:
        if p['Ny'] is None:
            p['Width']  = 1
            p['Ny'] = 1
        else:
            if   dx is not None: p['Width']  = p['Ny']*p['Length']/p['Nx']
            elif dz is not None: p['Width']  = p['Ny']*p['Height']/p['Nz']
            ju.logger.info('ded_type.validate: Setting W {:7.4f}'.format(p['Width']))

    if p['Height'] is None:
        if p['Nz'] is None:
            p['Height']  = 1
            p['Nz'] = 1
        else:
            if   dx is not None: p['Height']  = p['Ny']*p['Length']/p['Nx']
            elif dy is not None: p['Height']  = p['Nz']*p['Height']/p['Nz']
            ju.logger.info('ded_type.validate: Setting H {:7.4f}'.format(p['Height']))
           

 
    L=p['Length']
    W=p['Width']
    H=p['Height']
    
    if p['Nx'] is None:
        if   dy is not None: p['Nx']=np.int64(L/dy)
        elif dz is not None: p['Nx']=np.int64(L/dz)

    dx=p['Length']/p['Nx']

    if p['Ny'] is None:
        if   p['Ty']=="Cheb": p['Ny']=np.int64(3/2*W/dx)
        else:                 p['Ny']=np.int64(    W/dx)
 
    if p['Nz'] is None:
        if   p['Tz']=="Cheb": p['Nz']=np.int64(3/2*H/dx)
        else:                 p['Nz']=np.int64(    H/dx)

    dim=np.int(p['Nx']>1) + np.int(p['Ny']>1) +np.int(p['Nz']>1)  
    if dim==1 and ju.mpisize>1:
        ju.logger.info('Dedalus cannot parallelise 1d problems')
        ju.set_error('1d problem')

    if p['Nx']>1: Nx=p['Nx']
    else:
        jpop(p,('Wx','Tx'))
        p['Length']=1.0
        p['Nx']=1

    if p['Ny']>1: Ny=p['Ny']
    else:
        jpop(p,('Wy','Ty'))
        p['Width']=1.0
        p['Ny']=1

    if p['Nz']>1: Nz=p['Nz']
    else:
        jpop(p,('Wz','Tz'))
        p['Height']=1.0
        p['Nz']=1

             
    if p['sType']=='bg':
        p['Ny']=1
        p['Nz']=1
        p['Width']=1
        p['Height']=1

    ju.check_abort()
    
    p=remove_zeros(p,('bV','SV','bVh','SVh','S','clipu','dtcheck','dtjcheck','dtb','dtu','dtv','dtw','dtp','dtnoise','dtstats','dtx','dty','dtz','dtxy','dtxz','dtyz','dtxyz','dtmomb','dtavrg','dtleft','dtright','dtpm','dtgc','dtfpid','dtforce','dtslice','fpid','PIDX','inoise','bnoise','snoise','slicex','slicey','slicez','dtj','dtjx','dtjy','dtjz','dtjr','dtjyz','dtjxz','dtjxyz','dtja','dtjar','dtjax','dtjay','dtjaz','dtjaxy','dtjayz','dtjaxz','dtjaxyz','dtjb','dtjd','dtjs','dtjp','dtju','dtjv','dtjw','dtjWx','dtjWy','dtjWz','dtjS','dtjE','dtjaS','dtjaE','dtjsev','dtjavar','dnoise','Noise','wdivxl','wdivxr'))
    p=remove_None(p)

    if 'Gravity'  not in p: jpop(p,('Angle',))
    #if 'Forcing'  not in p: jpop(p,('hu',))
    if 'Buoyancy' not in p: jpop(p,('blbc','brbc','Scb','Peb','bV','clipB','fbmin','fbmax','forceb'))
    if 'Sediment' not in p: jpop(p,('slbc','srbc','Scs','Pes','SV','clipS','fsmin','fsmax','forces'))
    if 'PIDX'     not in p: jpop(p,('PIDG','PIDD','PIDI','PIDP','PIDT','PIDIT','PIDDD','PIDST','PIDS1','PIDS2','Umax','Umin','gmin','gmax'))

    jpop(p,('fpid','fpidb','fpids','fpidu','fpidv','fpidw','xn1','xn2','xn3','xn4','VMEM','alpha','avrg'))
    if 'fbmult' in p:
        if p['fbmult']: jpop(p,('fbmin','fbmax'))

    nmstr   = list(set(nmstr)   & set(p.keys()))
    nmint   = list(set(nmint)   & set(p.keys()))
    nmbool  = list(set(nmbool)  & set(p.keys()))
    nmfloat = list(set(nmfloat) & set(p.keys()))
    missing = list(set(p.keys()) - set(nmstr) - set(nmint) - set(nmbool) - set(nmfloat))

    for a in missing: ju.logger.info('Missing: {} {}'.format(a,p[a]))
    if len(missing)>0: quit()
    for a in nmstr  : p[a] = np.str(p[a])
    for a in nmint  : p[a] = np.int64(p[a])
    for a in nmbool : p[a] = tobool(p[a])
    for a in nmfloat: p[a] = np.float64(p[a])




    for a in ('AAS', 'AAJ','AA'):
        if a not in p: p[a] = None

    if p['AAS'] is None: p['AAS'] = 1
    if p['AAJ'] is None: p['AAJ'] = 1
    if p['AA']  is None: p['AA']  = 3/2

    topu  = get_bool_param(p,'topu')
    topr  = get_bool_param(p,'topr')
    topd  = get_bool_param(p,'topd')
    topq  = get_bool_param(p,'topq')
    topdd = get_bool_param(p,'topdd')
    dd    = get_bool_param(p,'dd')
    db    = get_bool_param(p,'db')
    ds    = get_bool_param(p,'ds')

    if dd and 'ddT' not in p: p['ddT']=1
    if db and 'dbT' not in p: p['dbT']=1
    if ds and 'dsT' not in p: p['dsT']=1
     
    if p['sType']=='pm':
        if np.int(topu)+np.int(topr)+np.int(topdd)+np.int(topd) +np.int(topq) !=1:
            ju.logger.info('Exactly one of topu:{} topr:{} topdd:{} topdq:{} must be specified'.format(topu,topr,topd,topdd,topq))
            ju.abort_run()
        if not topdd: jpop(p,'topddM')
        if not dd:  jpop(p,('wdivxl','wdivxr','divxI'))
        if not topd and not topu: jpop(p,'topT')
    else: jpop(p,('wdivxl','wdivxr','ddSc','topT','topu','topr','topd','topdd'))


    ju.logger.info('')
    ju.logger.info('String parameters')
    for a in nmstr: ju.logger.info('  {:10s}: {}'.format(a,p[a]))
    ju.logger.info('')
    ju.logger.info('Integer parameters')
    for a in nmint: ju.logger.info('  {:10s}: {}'.format(a,p[a]))
    ju.logger.info('')
    ju.logger.info('Boolean parameters')
    for a in nmbool: ju.logger.info('  {:10s}: {}'.format(a,p[a]))
    ju.logger.info('')
    ju.logger.info('Floating point parameters')
    for a in nmfloat: ju.logger.info('  {:10s}: {:11.5f}'.format(a,p[a])) 
    ju.logger.info('')

    return(p)

 
def get_param(p,n,d=None):
    if n in p: s=p[n]
    else: s=d
    return(s)

def get_bool_param(p,n):
    if n in p: s=p[n]
    else: s=False
    return(s)

