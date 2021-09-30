# Set the output data to writing buoyancy every 0.1 time units
# scale increases the resolution of the output datafiles
import jutil as ju
from polytrope.tools.checkpointing import Checkpoint

def fdtmomb(s,problem):
    s.add_task("integ(b,   )", scales=1, name='b')
    s.add_task("integ(b*x  )", scales=1, name='bx')
    s.add_task("integ(b*z  )", scales=1, name='bz')
    s.add_task("integ(b*u*x)", scales=1, name='bux')
    s.add_task("integ(b*u*z)", scales=1, name='buz')
    s.add_task("integ(b*w*x)", scales=1, name='bwx')
    s.add_task("integ(b*w*z)", scales=1, name='bwz')
    s.add_task("integ(b*x*x)", scales=1, name='bxx')
    s.add_task("integ(b*z*z)", scales=1, name='bzz')
    s.add_task("integ(b*x*z)", scales=1, name='bxz')
    if 'v' in problem.variables:
        s.add_task("integ(b*y  )", scales=1, name='by')
        s.add_task("integ(b*y*y)", scales=1, name='byy')
        s.add_task("integ(b*x*y)", scales=1, name='bxy')
        s.add_task("integ(b*y*z)", scales=1, name='byz')
        s.add_task("integ(b*u*y)", scales=1, name='buy')
        s.add_task("integ(b*w*y)", scales=1, name='bwy')
        s.add_task("integ(b*v*x)", scales=1, name='bvx')
        s.add_task("integ(b*v*y)", scales=1, name='bvy')
        s.add_task("integ(b*v*z)", scales=1, name='bvz')
    return(s)
        

def fdtint(s,problem,d1,d2,d3):
    v=problem.variables
    if 'u' not in v:
        if d1=='x': d1=None
        if d2=='x': d2=None
        if d3=='x': d3=None

    if 'v' not in v:
        if d1=='y': d1=None
        if d2=='y': d2=None
        if d3=='y': d3=None

    if 'w' not in v:
        if d1=='z': d1=None
        if d2=='z': d2=None
        if d3=='z': d3=None
        
    if d1 is None and d2 is None and d3 is None:
        s1=""
        s2=""
    else:
        s1 = "integ("
        s2 = ""
        if d1 is not None: s2 += ",'"+d1+"'"
        if d2 is not None: s2 += ",'"+d2+"'"
        if d3 is not None: s2 += ",'"+d3+"'"
        s2 += ")"

    #ju.logger.info('fdtint d1:{} d2:{} d3:{}, "{}"'.format(d1,d2,d3,s1 + "f" + s2))

    if 'SR' in problem.substitutions: s.add_task(s1 + "SR" + s2, scales=1, name='S')
    if 'E' in problem.substitutions: s.add_task(s1 + "E"  + s2, scales=1, name='E')
    if 'b' in v:              s.add_task(s1 + "b  "     + s2, scales=1, name='b')
    if 's' in v:              s.add_task(s1 + "s  "     + s2, scales=1, name='s')
    if 'u' in v:              s.add_task(s1 + "u  "     + s2, scales=1, name='u')
    if 'v' in v:              s.add_task(s1 + "v  "     + s2, scales=1, name='v')
    if 'w' in v:              s.add_task(s1 + "w  "     + s2, scales=1, name='w')
    if 'p' in v:              s.add_task(s1 + "p  "     + s2, scales=1, name='p')
    if 'b' in v:              s.add_task(s1 + "b*b"     + s2, scales=1, name='bb')
    if 's' in v:              s.add_task(s1 + "s*s"     + s2, scales=1, name='ss')
    if 'u' in v:              s.add_task(s1 + "u*u"     + s2, scales=1, name='uu')
    if 'v' in v:              s.add_task(s1 + "v*v"     + s2, scales=1, name='vv')
    if 'w' in v:              s.add_task(s1 + "w*w"     + s2, scales=1, name='ww')
    if 'u' in v and 'b' in v: s.add_task(s1 + "u*b"     + s2, scales=1, name='ub')
    if 'u' in v and 's' in v: s.add_task(s1 + "u*s"     + s2, scales=1, name='us')
    if 'v' in v and 'b' in v: s.add_task(s1 + "v*b"     + s2, scales=1, name='vb')
    if 'v' in v and 's' in v: s.add_task(s1 + "v*s"     + s2, scales=1, name='vs')
    if 'w' in v and 's' in v: s.add_task(s1 + "w*b"     + s2, scales=1, name='wb')
    if 'w' in v and 's' in v: s.add_task(s1 + "w*s"     + s2, scales=1, name='ws')
    if 'u' in v and 'v' in v: s.add_task(s1 + "u*v"     + s2, scales=1, name='uv')
    if 'u' in v and 'w' in v: s.add_task(s1 + "u*w"     + s2, scales=1, name='uw')
    if 'v' in v and 'w' in v: s.add_task(s1 + "v*w"     + s2, scales=1, name='vw')
    if 'v' in v and 'w' in v: s.add_task(s1 + "vorx**2" + s2, scales=1, name='Ex')
    if 'u' in v and 'w' in v: s.add_task(s1 + "vory**2" + s2, scales=1, name='Ey')
    if 'u' in v and 'v' in v: s.add_task(s1 + "vorz**2" + s2, scales=1, name='Ez')
    return(s)

def fdtavrg2(s,problem):
    v=problem.variables
    
    if 'v' in v:
        s1="integ("
        s2=",'y')"
    else:
        s1=""
        s2=""
    a=list(set(v) & set(('b','s','u','v','w','p','ab','as','au','aw','ap','abu','abw','asu','asw','auu','avv','aww','auw','avw','auv','aur','arr','aru')))
    #a=list(set(v) & set(('b','s','as','au','aw','ap','abw','asu','asw','auw','avw','auv','aur','arr','aru')))
    for b in a:
        ts=str(s1+b+s2)
        ju.logger.info('s.add_task("{:12s}", scales=1, name="{:3s}")'.format(ts,b))
        s.add_task(ts, scales=1, name=b)
        #ju.logger.info(type(s1+b+s2))
        #ju.logger.info(type(b))
        #ju.logger.info(type('ab'))
    return(s)
            
def fdtavrg(s,problem):
    v=problem.variables
    
    if 'v' in v:                                                               
        if 's'   in v: s.add_task("integ(s,  'y')", scales=1, name='s')      
        if 'b'   in v: s.add_task("integ(b,  'y')", scales=1, name='b')                   
        if 'u'   in v: s.add_task("integ(u,  'y')", scales=1, name='u')                   
        if 'v'   in v: s.add_task("integ(v,  'y')", scales=1, name='v')                   
        if 'w'   in v: s.add_task("integ(w,  'y')", scales=1, name='w')                       
        if 'p'   in v: s.add_task("integ(p,  'y')", scales=1, name='p')                       
        if 'ab'  in v: s.add_task("integ(ab, 'y')", scales=1, name='ab')       
        if 'as'  in v: s.add_task("integ(as, 'y')", scales=1, name='as')       
        if 'au'  in v: s.add_task("integ(au, 'y')", scales=1, name='au')       
        if 'aw'  in v: s.add_task("integ(aw, 'y')", scales=1, name='aw')       
        if 'ap'  in v: s.add_task("integ(ap, 'y')", scales=1, name='ap')       
        if 'abu' in v: s.add_task("integ(abu,'y')", scales=1, name='abu')      
        if 'abw' in v: s.add_task("integ(abw,'y')", scales=1, name='abw')      
        if 'asu' in v: s.add_task("integ(asu,'y')", scales=1, name='asu')      
        if 'asw' in v: s.add_task("integ(asw,'y')", scales=1, name='asw')      
        if 'auu' in v: s.add_task("integ(auu,'y')", scales=1, name='auu')      
        if 'avv' in v: s.add_task("integ(avv,'y')", scales=1, name='avv')      
        if 'aww' in v: s.add_task("integ(aww,'y')", scales=1, name='aww')      
        if 'auw' in v: s.add_task("integ(auw,'y')", scales=1, name='auw')      
        if 'avw' in v: s.add_task("integ(avw,'y')", scales=1, name='avw')      
        if 'auv' in v: s.add_task("integ(auv,'y')", scales=1, name='auv')      
        if 'aur' in v: s.add_task("integ(aur,'y')", scales=1, name='aur')      
        if 'arr' in v: s.add_task("integ(arr,'y')", scales=1, name='arr')      
        if 'aru' in v: s.add_task("integ(aru,'y')", scales=1, name='aru')      
    else:                                                                      
        if 's'   in v: s.add_task("s",   scales=1, name='s')                       
        if 'b'   in v: s.add_task("b",   scales=1, name='b')                                    
        if 'u'   in v: s.add_task("u",   scales=1, name='u')                                    
        if 'w'   in v: s.add_task("w",   scales=1, name='w')                                    
        if 'p'   in v: s.add_task("p",   scales=1, name='p')                                    
        if 'ab'  in v: s.add_task("ab",  scales=1, name='ab')                  
        if 'au'  in v: s.add_task("au",  scales=1, name='au')                  
        if 'aw'  in v: s.add_task("aw",  scales=1, name='aw')                  
        if 'ap'  in v: s.add_task("ap",  scales=1, name='ap')                  
        if 'abu' in v: s.add_task("abu", scales=1, name='abu')                 
        if 'abw' in v: s.add_task("abw", scales=1, name='abw')                 
        if 'abs' in v: s.add_task("asu", scales=1, name='asu')                 
        if 'abs' in v: s.add_task("asw", scales=1, name='asw')                 
        if 'auu' in v: s.add_task("auu", scales=1, name='auu')                 
        if 'aww' in v: s.add_task("aww", scales=1, name='aww')                 
        if 'auw' in v: s.add_task("auw", scales=1, name='auw')                 
        if 'aur' in v: s.add_task("aur", scales=1, name='aur')                 
        if 'arr' in v: s.add_task("arr", scales=1, name='arr')                 
        if 'aru' in v: s.add_task("aru", scales=1, name='aru')                 
       
    return(s)
    
        
def fdtlr(s,problem,lr):
    v=problem.variables
    if 'v' in problem.variables:
        s1="integ("+lr+"("
        s2="),'y')"
    else:
        s1=lr+"("
        s2=")"
    #ju.logger.info('fdtlr {}: {}'.format(lr,s1+"b"+s2))
    if 'b'  in v: s.add_task(s1+"b"+s2,         scales=1, name='b')
    if 's'  in v: s.add_task(s1+"s"+s2,         scales=1, name='s')
    if 'u'  in v: s.add_task(s1+"u"+s2,         scales=1, name='u')
    if 'v'  in v: s.add_task(s1+"v"+s2,         scales=1, name='v')
    if 'w'  in v: s.add_task(s1+"w"+s2,         scales=1, name='w')
    if 'p'  in v: s.add_task(s1+"p"+s2,         scales=1, name='p')

    if 'u'  in v: s.add_task(s1+"u*u"+s2,       scales=1, name='uu')
    if 'v'  in v: s.add_task(s1+"v*v"+s2,       scales=1, name='vv')
    if 'w'  in v: s.add_task(s1+"w*w"+s2,       scales=1, name='ww')
    if 'u'  in v and 'v' in v: s.add_task(s1+"u*v"+s2, scales=1, name='uv')
    if 'u'  in v and 'w' in v: s.add_task(s1+"u*w"+s2, scales=1, name='uw')
    if 'v'  in v and 'w' in v: s.add_task(s1+"v*w"+s2, scales=1, name='vw')
    if 'u'  in v and 'b' in v: s.add_task(s1+"u*b"+s2, scales=1, name='ub')
    if 'u'  in v and 's' in v: s.add_task(s1+"u*s"+s2, scales=1, name='us')

    if 'bz' in v: s.add_task(s1+"bz"+s2,        scales=1, name='bdz')
    if 'sz' in v: s.add_task(s1+"sz"+s2,        scales=1, name='sdz')
    if 'uz' in v: s.add_task(s1+"uz"+s2,        scales=1, name='udz')
    if 'vz' in v: s.add_task(s1+"vz"+s2,        scales=1, name='vdz')
    if 'wz' in v: s.add_task(s1+"wz"+s2,        scales=1, name='wdz')
    if 'p'  in v: s.add_task(s1+"dz(p)"+s2,     scales=1, name='pdz')
    if 'u'  in v: s.add_task(s1+"dz(u*u)"+s2,   scales=1, name='uudz')
    if 'v'  in v: s.add_task(s1+"dz(v*v)"+s2,   scales=1, name='vvdz')
    if 'w'  in v: s.add_task(s1+"dz(w*w)"+s2,   scales=1, name='wwdz')

    if 'u'  in v: s.add_task(s1+"dz(dz(u))"+s2, scales=1, name='udzz')
    if 'v'  in v: s.add_task(s1+"dz(dz(v))"+s2, scales=1, name='vdzz')
    if 'w'  in v: s.add_task(s1+"dz(dz(w))"+s2, scales=1, name='wdzz')
    if 'u'  in v and 'b' in v: s.add_task(s1+"dz(u*b)"+s2, scales=1, name='ubdz')
    if 'u'  in v and 's' in v: s.add_task(s1+"dz(u*s)"+s2, scales=1, name='usdz')
    if 'u'  in v and 'v' in v: s.add_task(s1+"dz(u*v)"+s2, scales=1, name='uvdz')
    if 'u'  in v and 'w' in v: s.add_task(s1+"dz(u*w)"+s2, scales=1, name='uwdz')
    if 'v'  in v and 'w' in v: s.add_task(s1+"dz(v*w)"+s2, scales=1, name='vwdz')
    return(s)

def fdtfpid(s,problem):
    v=problem.variables
    if 'fpidu' in problem.variables: s.add_task("fpidu", scales=1, name='fpidu')
    if 'fpidv' in problem.variables: s.add_task("fpidv", scales=1, name='fpidv')
    if 'fpidw' in problem.variables: s.add_task("fpidw", scales=1, name='fpidw')
    if 'fpidb' in problem.variables: s.add_task("fpidb", scales=1, name='fpidb')
    if 'fpids' in problem.variables: s.add_task("fpids", scales=1, name='fpids')
    return(s)

def fdtavrggc(s,problem):
    if 'v' in problem.variables:
        s.add_task("integ(left(u),           'y')", scales=1, name='Lu')
        s.add_task("integ(left(ab         ), 'y')", scales=1, name='Lab')
        s.add_task("integ(left(dx(ap)     ), 'y')", scales=1, name='Lapdx')
        s.add_task("integ(left(dx(auu)    ), 'y')", scales=1, name='Lauudx')
        s.add_task("integ(left(dz(auw)    ), 'y')", scales=1, name='Lauwdz')
        s.add_task("integ(left(dx(dx(au)) ), 'y')", scales=1, name='Laudxx')
        s.add_task("integ(left(dz(dz(au)) ), 'y')", scales=1, name='Laudzz')
        
        s.add_task("integ(right(u),          'y')", scales=1, name='Ru')
        s.add_task("integ(right(ab        ), 'y')", scales=1, name='Rab')
        s.add_task("integ(right(dx(auu)   ), 'y')", scales=1, name='Rauudx')
        s.add_task("integ(right(dz(auw)   ), 'y')", scales=1, name='Rauwdz')
        s.add_task("integ(right(dx(ap)    ), 'y')", scales=1, name='Rapdx')
        s.add_task("integ(right(dx(dx(au))), 'y')", scales=1, name='Raudxx')
        s.add_task("integ(right(dz(dz(au))), 'y')", scales=1, name='Raudzz')
        
        s.add_task("integ(ab,            'y','z')", scales=1, name='ab')
        s.add_task("integ(dz(ap    ),    'y','z')", scales=1, name='apdz')
        s.add_task("integ(dx(auw   ),    'y','z')", scales=1, name='auwdx')
        s.add_task("integ(dz(aww   ),    'y','z')", scales=1, name='awwdz')
        s.add_task("integ(dx(dx(aw)),    'y','z')", scales=1, name='awdxx')
    else:
        s.add_task("left(u),        ", scales=1, name='Lu')
        s.add_task("left(ab         ", scales=1, name='Lab')
        s.add_task("left(dx(ap)     ", scales=1, name='Lapdx')
        s.add_task("left(dx(auu)    ", scales=1, name='Lauudx')
        s.add_task("left(dz(auw)    ", scales=1, name='Lauwdz')
        s.add_task("left(dx(dx(au)) ", scales=1, name='Laudxx')
        s.add_task("left(dz(dz(au)) ", scales=1, name='Laudzz')
        
        s.add_task("right(u),       ", scales=1, name='Ru')
        s.add_task("right(ab        ", scales=1, name='Rab')
        s.add_task("right(dx(auu)   ", scales=1, name='Rauudx')
        s.add_task("right(dz(auw)   ", scales=1, name='Rauwdz')
        s.add_task("right(dx(ap)    ", scales=1, name='Rapdx')
        s.add_task("right(dx(dx(au))", scales=1, name='Raudxx')
        s.add_task("right(dz(dz(au))", scales=1, name='Raudzz')
    
        s.add_task("integ(ab,        'z')", scales=1, name='ab')
        s.add_task("integ(dz(ap    ),'z')", scales=1, name='apdz')
        s.add_task("integ(dx(auw   ),'z')", scales=1, name='auwdx')
        s.add_task("integ(dz(aww   ),'z')", scales=1, name='awwdz')
        s.add_task("integ(dx(dx(aw)),'z')", scales=1, name='awdxx')
       
    return(s)



def fdtgc(s,problem):
    if 'v' in problem.variables:
        s.add_task("integ(integ(b,'y')**2,     'z')", scales=1, name='bs')
        s.add_task("integ(integ(p,'y')**2,     'z')", scales=1, name='ps')
        s.add_task("integ(integ(u,'y')**2,     'z')", scales=1, name='us')
        s.add_task("integ(b*z,             'y','z')", scales=1, name='bz')
        s.add_task("integ(p*z,             'y','z')", scales=1, name='pz')
        s.add_task("integ(u*z,             'y','z')", scales=1, name='uz')
        s.add_task("integ(vory**2,         'y','z')", scales=1, name='Ey')
        s.add_task("integ(vory*u,          'y','z')", scales=1, name='qwy')
    else:
        s.add_task("integ(b**2,     'z')", scales=1, name='bs')
        s.add_task("integ(p**2,     'z')", scales=1, name='ps')
        s.add_task("integ(u**2,     'z')", scales=1, name='us')
        s.add_task("integ(b*z,      'z')", scales=1, name='bz')
        s.add_task("integ(p*z,      'z')", scales=1, name='pz')
        s.add_task("integ(u*z,      'z')", scales=1, name='uz')
        s.add_task("integ(vory**2,  'z')", scales=1, name='Ey')
        s.add_task("integ(vory*u,   'z')", scales=1, name='qwy')
    return(s)

def fdtgcp(s,problem):
    s.add_task("integ(             PIN(integ(u,'y')),'z')", scales=1, name='us')
    s.add_task("integ(integ(b,'y')*PIN(integ(u,'y')),'z')", scales=1, name='bus')
    s.add_task("integ(             POS(integ(u,'y')),'z')", scales=1, name='up')
    s.add_task("integ(integ(b,'y')*POS(integ(u,'y')),'z')", scales=1, name='bup')
    s.add_task("integ(           z*POS(integ(u,'y')),'z')", scales=1, name='zup')
    s.add_task("integ(             POS(integ(w,'y')),'z')", scales=1, name='wp')
    s.add_task("integ(integ(b,'y')*POS(integ(w,'y')),'z')", scales=1, name='bwp')
    s.add_task("integ(           z*POS(integ(w,'y')),'z')", scales=1, name='zwp')
    #    s.add_task("integ(POS(u),          'y','z')", scales=1, name='up')
    #    s.add_task("integ(POS(u)*b,        'y','z')", scales=1, name='bup')
    #    s.add_task("integ(POS(u)*z,        'y','z')", scales=1, name='zup')
    #    s.add_task("integ(POS(w),          'y','z')", scales=1, name='wp')
    #    s.add_task("integ(POS(w)*b,        'y','z')", scales=1, name='bwp')
    #    s.add_task("integ(POS(w)*z,        'y','z')", scales=1, name='zwp')
    return(s)

def fdtslice(s,problem,d,e):
    v=problem.variables
    if 'u' in v: s.add_task('interp(u, ' + d + '={})'.format(e), scales=1, name='u')
    if 'v' in v: s.add_task('interp(v, ' + d + '={})'.format(e), scales=1, name='v')
    if 'w' in v: s.add_task('interp(w, ' + d + '={})'.format(e), scales=1, name='w')
    if 'b' in v: s.add_task('interp(b, ' + d + '={})'.format(e), scales=1, name='b')
    if 's' in v: s.add_task('interp(s, ' + d + '={})'.format(e), scales=1, name='s')
    return(s)

def fdtpm(s,problem):
    s.add_task("integ(b,            'y','z')", scales=1, name='b')
    s.add_task("integ(b*(y*y+z*z),  'y','z')", scales=1, name='brr')
    s.add_task("integ(u,            'y','z')", scales=1, name='u')
    s.add_task("integ(u*(y*y+z*z),  'y','z')", scales=1, name='urr')
    s.add_task("integ(p,            'y','z')", scales=1, name='p')
    s.add_task("integ(p*(y*y+z*z),  'y','z')", scales=1, name='prr')
    s.add_task("integ(b*u,          'y','z')", scales=1, name='bu')
    s.add_task("integ(b*u*(y*y+z*z),'y','z')", scales=1, name='burr')
    s.add_task("integ(y*v,          'y','z')", scales=1, name='vy')
    s.add_task("integ(y*v*(y*y+z*z),'y','z')", scales=1, name='vyrr')
    s.add_task("integ(z*w,          'y','z')", scales=1, name='wz')
    s.add_task("integ(z*w*(y*y+z*z),'y','z')", scales=1, name='wzrr')
    s.add_task("integ(u*u,          'y','z')", scales=1, name='uu')
    s.add_task("integ(u*u*(y*y+z*z),'y','z')", scales=1, name='uurr')
    s.add_task("integ(u*v*y,        'y','z')", scales=1, name='uvy')
    s.add_task("integ(u*w*z,        'y','z')", scales=1, name='uwz')
    return(s)

def fdtforce(s,problem):
    v=problem.parameters
    if 'psi'    in v: s.add_task("psi",     scales=1, name='psi')
    if 'psix'   in v: s.add_task("psix",    scales=1, name='psix')
    if 'psiy'   in v: s.add_task("psiy",    scales=1, name='psiy')
    if 'psiz'   in v: s.add_task("psiz",    scales=1, name='psiz')
    if 'fdx'    in v: s.add_task("fdx",     scales=1, name='fdx')
    if 'fdy'    in v: s.add_task("fdy",     scales=1, name='fdy')
    if 'fdz'    in v: s.add_task("fdz",     scales=1, name='fdz')
    if 'fb'     in v: s.add_task("fb",      scales=1, name='fb')
    if 'fs'     in v: s.add_task("fs",      scales=1, name='fs')
    if 'fu'     in v: s.add_task("fu",      scales=1, name='fu')
    if 'fv'     in v: s.add_task("fv",      scales=1, name='fv')
    if 'fw'     in v: s.add_task("fw",      scales=1, name='fw')
    if 'wd'     in v: s.add_task("wd",      scales=1, name='wd')
    if 'wb'     in v: s.add_task("wb",      scales=1, name='wb')
    if 'ws'     in v: s.add_task("ws",      scales=1, name='ws')
    if 'wu'     in v: s.add_task("wu",      scales=1, name='wu')
    if 'wv'     in v: s.add_task("wv",      scales=1, name='wv')
    if 'ww'     in v: s.add_task("ww",      scales=1, name='ww')
    if 'fd'     in v: s.add_task("fd",      scales=1, name='fd')
    if 'wnoise' in v: s.add_task("wnoise",  scales=1, name='wnoise')
    if 'fpidu'  in v: s.add_task("fpidu",   scales=1, name='fpidu')
    if 'fpidv'  in v: s.add_task("fpidv",   scales=1, name='fpidv')
    if 'fpidw'  in v: s.add_task("fpidw",   scales=1, name='fpidw')
    if 'fd' in v: 
        ts='fd'
        if 'fu' in v: ts += '-dx(fu)'
        if 'fv' in v: ts += '-dx(fv)'
        if 'fw' in v: ts += '-dx(fw)'
        if len(ts)>2: s.add_task(ts,scales=1, name='div')
    if len(s.tasks)==0: s=None
    return(s)
 
def set_analysis_tasks(ddir,solver,problem,max_writes,param,mode):
    var=problem.variables
    dtb =       0
    dts =       0
    dtu =       0
    dtv =       0
    dtw =       0
    dtp =       0
    dtnoise =   0
    dtx =       0
    dty =       0
    dtz =       0
    dtxy =      0
    dtxz =      0
    dtyz =      0
    dtxyz =     0
    dtmomb =    0
    dtavrg =    0
    dtleft =    0
    dtright =   0
    dtpm =      0
    dtgc =      0
    dtforce =   0
    dtstats =   0
    dtslice =   0
    dtcheck =   0
    dtfpid  =   0
    slicex  = None
    slicey  = None
    slicez  = None
    
    if 'dtcheck' in param: dtcheck =  param['dtcheck']
    if 'dts'     in param: dts =      param['dts']
    if 'dtb'     in param: dtb =      param['dtb']
    if 'dtu'     in param: dtu =      param['dtu']
    if 'dtv'     in param: dtv =      param['dtv']
    if 'dtw'     in param: dtw =      param['dtw']
    if 'dtp'     in param: dtp =      param['dtp']
    if 'dtnoise' in param: dtnoise =  param['dtnoise']
    if 'dtx'     in param: dtx =      param['dtx']
    if 'dty'     in param: dty =      param['dty']
    if 'dtz'     in param: dtz =      param['dtz']
    if 'dtxy'    in param: dtxy =     param['dtxy']
    if 'dtxz'    in param: dtxz =     param['dtxz']
    if 'dtyz'    in param: dtyz =     param['dtyz']
    if 'dtxyz'   in param: dtxyz =    param['dtxyz']
    if 'dtmomb'  in param: dtmomb =   param['dtmomb']
    if 'dtavrg'  in param: dtavrg =   param['dtavrg']
    if 'dtleft'  in param: dtleft =   param['dtleft']
    if 'dtright' in param: dtright =  param['dtright']
    if 'dtpm'    in param: dtpm =     param['dtpm']
    if 'dtgc'    in param: dtgc =     param['dtgc']       
    if 'dtfpid'  in param: dtfpid =   param['dtfpid']       
    if 'dtslice' in param: dtslice =  param['dtslice']       
    if 'dtforce' in param: dtforce =  param['dtforce']       
    if 'dtstats' in param: dtstats =  param['dtstats']       
    if 'slicex'  in param: slicex  =  param['slicex']       
    if 'slicey'  in param: slicey  =  param['slicey']       
    if 'slicez'  in param: slicez  =  param['slicez']       
     
    analysis_tasks = []
 
    m0=max_writes
    m1=max_writes*10
    m2=max_writes*100
    m3=max_writes*1000

    ju.logger.info('Analysis tasks cadences')

    parallel = param['parallel']
    
    if dtcheck>0:
        ju.logger.info('  dtcheck:  {:3.0f} min'.format(dtcheck))
        checkpoint = Checkpoint(ddir)
        checkpoint.set_checkpoint(solver, wall_dt=dtcheck*60, mode=mode,parallel=parallel)
     
    if dtgc>0:
        ju.logger.info('  dtgc:      {:6.3f}'.format(dtgc))
        s = solver.evaluator.add_file_handler(ddir / 'gc',    sim_dt=dtgc, parallel=parallel,max_writes=m2)
        analysis_tasks.append(fdtgc(s,problem))

    if dtslice>0:
        ju.logger.info("  dtslice:  {:7.3f}".format(dtslice))
        if slicex is not None and 'u' in var:
            ju.logger.info("  slicex:   {:7.3f}".format(slicex))
            s = solver.evaluator.add_file_handler(ddir / 'slicex',    sim_dt=dtslice, parallel=parallel,max_writes=m2)
            analysis_tasks.append(fdtslice(s,problem,'x',slicex))
        if slicey is not None and 'v' in var:
            ju.logger.info("  slicey:   {:7.3f}".format(slicey))
            s = solver.evaluator.add_file_handler(ddir / 'slicey',    sim_dt=dtslice, parallel=parallel,max_writes=m2)
            analysis_tasks.append(fdtslice(s,problem,'y',slicey))
        if slicez is not None and 'w' in var:
            ju.logger.info("  slicez:   {:7.3f}".format(slicez))
            s = solver.evaluator.add_file_handler(ddir / 'slicez',    sim_dt=dtslice, parallel=parallel,max_writes=m2)
            analysis_tasks.append(fdtslice(s,problem,'z',slicez))

    if dtpm>0:
        ju.logger.info('  dtpm:      {:6.3f}'.format(dtpm))
        s = solver.evaluator.add_file_handler(ddir / 'pm',    sim_dt=dtpm, parallel=parallel,max_writes=m2)
        analysis_tasks.append(fdtpm(s,problem))
 
    if dtforce>0:
        ju.logger.info('  dtforce:   {:6.3f}'.format(dtforce))
        s = solver.evaluator.add_file_handler(ddir / 'force', sim_dt=1000000, parallel=parallel,max_writes=1,mode='overwrite')
        s=fdtforce(s,problem)     
        if s is not None: analysis_tasks.append(s)

    if dtmomb>0 and 'b' in var:
        ju.logger.info('  dtmomb:    {:6.3f}'.format(dtmomb))
        s = solver.evaluator.add_file_handler(ddir / 'momb',  sim_dt=dtmomb, parallel=parallel,max_writes=m3)
        analysis_tasks.append(fdtmomb(s,problem))

    if dtxyz>0:
        ju.logger.info('  dtxyz:     {:6.3f}'.format(dtxyz))
        s = solver.evaluator.add_file_handler(ddir / 'xyz',   sim_dt=dtxyz, parallel=parallel,max_writes=m3)
        analysis_tasks.append(fdtint(s,problem,'x','y','z'))
        
    if dtxy>0:
        ju.logger.info('  dtxy:      {:6.3f}'.format(dtxy))
        s = solver.evaluator.add_file_handler(ddir / 'xy',     sim_dt=dtxy, parallel=parallel,max_writes=m2)
        analysis_tasks.append(fdtint(s,problem,'x','y',None))

    if dtxz>0:
        ju.logger.info('  dtxz:      {:6.3f}'.format(dtxz))
        s = solver.evaluator.add_file_handler(ddir / 'xz',       sim_dt=dtxz, parallel=parallel,max_writes=m2)
        analysis_tasks.append(fdtint(s,problem,'x','z',None))

    if dtyz>0:
        ju.logger.info('  dtyz:      {:6.3f}'.format(dtyz))
        s = solver.evaluator.add_file_handler(ddir / 'yz',       sim_dt=dtyz, parallel=parallel,max_writes=m2)
        analysis_tasks.append(fdtint(s,problem,'y','z',None))
         
    if dtx>0:
        ju.logger.info('  dtx:       {:6.3f}'.format(dtx))
        s = solver.evaluator.add_file_handler(ddir / 'x',         sim_dt=dtx, parallel=parallel,max_writes=m1)
        analysis_tasks.append(fdtint(s,problem,'x',None,None))
            
    if dty>0:
        ju.logger.info('  dty:       {:6.3f}'.format(dty))
        s = solver.evaluator.add_file_handler(ddir / 'y',         sim_dt=dty, parallel=parallel,max_writes=m1)
        analysis_tasks.append(fdtint(s,problem,'y',None,None))
    
    if dtz>0:
        ju.logger.info('  dtz:       {:6.3f}'.format(dtz))
        s = solver.evaluator.add_file_handler(ddir / 'z',         sim_dt=dtz, parallel=parallel,max_writes=m1)
        analysis_tasks.append(fdtint(s,problem,'z',None,None))

    if dtright>0:
        ju.logger.info('  dtright:   {:6.3f}'.format(dtright))
        s = solver.evaluator.add_file_handler(ddir / 'right', sim_dt=dtright, parallel=parallel,max_writes=m2)
        analysis_tasks.append(fdtlr(s,problem,'right'))
       
    if dtleft>0:
        ju.logger.info('  dtleft:    {:6.3f}'.format(dtleft))
        s = solver.evaluator.add_file_handler(ddir / 'left',   sim_dt=dtleft, parallel=parallel,max_writes=m2)
        analysis_tasks.append(fdtlr(s,problem,'left'))
        
    if dtavrg>0:
        ju.logger.info('  dtavrg:    {:6.3f}'.format(dtavrg))
        s = solver.evaluator.add_file_handler(ddir / 'avrg',   sim_dt=dtavrg, parallel=parallel,max_writes=m0)
        analysis_tasks.append(fdtavrg(s,problem))  
        #s = solver.evaluator.add_file_handler(ddir / 'avrg2',   sim_dt=dtavrg, parallel=parallel,max_writes=m0)
        #analysis_tasks.append(fdtavrg2(s,problem))  
       
    if dtfpid>0: 
        ju.logger.info('  dtfpid:    {:6.3f}'.format(dtpid))
        s = solver.evaluator.add_file_handler(ddir / 'fpid',   sim_dt=dtfpid, parallel=parallel,max_writes=1)
        analysis_tasks.append(fdtfpid(s,problem))
    
    if dtnoise and 'noisef' in problem.variables:
        ju.logger.info('  dtnoise:   {:6.3f}'.format(dtnoise))
        s = solver.evaluator.add_file_handler(ddir / 'noise',  sim_dt=dtnoise, parallel=parallel,max_writes=1)
        s.add_task("noisef", scales=1, name='noisef')
        s.add_task("noisef2",scales=1, name='noisef2')
        s.add_task("noisex", scales=1, name='noisex')
        s.add_task("noisey", scales=1, name='noisey')
        s.add_task("noisez", scales=1, name='noisez')
        s.add_task("wnoise", scales=1, name='wnoise')
        analysis_tasks.append(s)
        

    if dtu>0 and 'u' in problem.variables:
        ju.logger.info('  dtu:       {:6.3f}'.format(dtu))
        s = solver.evaluator.add_file_handler(ddir / 'u',          sim_dt=dtu, parallel=parallel,max_writes=1)
        s.add_task("u", scales=1, name='u')
        analysis_tasks.append(s)
        
    if dtv>0 and 'v' in problem.variables:
        ju.logger.info('  dtv:       {:6.3f}'.format(dtv))
        s = solver.evaluator.add_file_handler(ddir / 'v',          sim_dt=dtv, parallel=parallel,max_writes=1)
        s.add_task("v", scales=1, name='v')
        analysis_tasks.append(s)
    
    if dtw>0 and 'w' in problem.variables:
        ju.logger.info('  dtw:       {:6.3f}'.format(dtw))
        s = solver.evaluator.add_file_handler(ddir / 'w',          sim_dt=dtw, parallel=parallel,max_writes=1)
        s.add_task("w", scales=1, name='w')
        analysis_tasks.append(s)
    
    if dtp>0 and 'p' in problem.variables:
        ju.logger.info('  dtp:       {:6.3f}'.format(dtp))
        s = solver.evaluator.add_file_handler(ddir / 'p',          sim_dt=dtp, parallel=parallel,max_writes=1)
        s.add_task("p", scales=1, name='p')
        analysis_tasks.append(s)
        
    if dtb>0 and 'b' in problem.variables:
        ju.logger.info('  dtb:       {:6.3f}'.format(dtb))
        s = solver.evaluator.add_file_handler(ddir / 'b',          sim_dt=dtb, parallel=parallel,max_writes=1)
        s.add_task("b", scales=1, name='b')
        analysis_tasks.append(s)

    if dts>0 and 's' in problem.variables:
        ju.logger.info('  dts:       {:6.3f}'.format(dts))
        s = solver.evaluator.add_file_handler(ddir / 's',          sim_dt=dtb, parallel=parallel,max_writes=1)
        s.add_task("s", scales=1, name='s')
        analysis_tasks.append(s)
           
    if     dtstats>0: ju.logger.info('  dtstats:   {:6.3f}'.format(dtstats))


    #ju.logger.info('  dir(analysis_tasks):   {}'.format(dir(analysis_tasks)))

    
    if False:
        for j in range(len(analysis_tasks)):
            task=analysis_tasks[j]
            if task is None:
                ju.logger.info('  task {:2d}: None'.format(j))
                ju.abort=True
                ju.check_abort()
                ju_test_abort()
            else:
                ju.logger.info('  task {:2d}: dt:{} path:{}'.format(j,task.sim_dt,task.base_path))      
                #for v in vars(task): print('{}: {}'.format(v,eval('task.'+v)))
                #ju.logger.info('  task.current_path     {}'.format(task.current_path))    
                #ju.logger.info('  task.iter             {}'.format(task.iter))             
                #ju.logger.info('  task.max_size         {}'.format(task.max_size))         
                #ju.logger.info('  task.max_writes       {}'.format(task.max_writes))       
                #ju.logger.info('  task.parallel         {}'.format(task.parallel))         
                #ju.logger.info('  task.sim_dt           {}'.format(task.sim_dt))
                #ju.logger.info('  task.total_write_num  {}'.format(task.total_write_num))  
                #ju.logger.info('  task.wall_dt          {}'.format(task.wall_dt))
                #ju.logger.info('  len(task.tasks)            {}'.format(len(task.tasks)))              
                for j in range(len(task.tasks)):
                    t=task.tasks[j]
                    ju.logger.info('    {:2d}: {:4s}: {}'.format(j,t['name'],t['operator']))             
                    #   ju.logger.info('      t[scales]:   {}'.format(t['scales']))              
                    #    ju.logger.info('  t            {}'.format(t))              
                    #    ju.logger.info('  type(t)      {}'.format(type(t)))              
                    #     ju.logger.info('  t[layout]:   {}'.format(t['layout']))              
                    # ju.logger.info('  task.tasks.scales     {}'.format(task.tasks.scales))              
                    # ju.logger.info('  task.tasks.max_writes {}'.format(task.tasks.max_writes))              

   
    return(analysis_tasks)


def save_final(ddir,solver,parallel,dt):
    ju.logger.info("Final step to save state")
    final_checkpoint = Checkpoint(ddir, checkpoint_name='final')
    final_checkpoint.set_checkpoint(solver, wall_dt=0.1, mode="overwrite", parallel=parallel)
    solver.step(dt)

