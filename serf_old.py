                  
# Symmetric error functions with zero mean and range of 1
# Along with h and w derivatives for curve fitting
def sserfz(xx,hh,w,HH):
    #[g gh gw]=serfz(x,h,w,H) zero mean symmetric error function  C+ A*(erf((h - y)/w) + erf((h + y)/w))
    # where C is such that int(g,y=0..H)=0
    # and A is such that g(H)-g(0)=1 # g=f/d
    
    h=hh/w
    H=HH/w
    x=xx/w
    

    px1  = np.exp(-(h - x)**2)/spi/H
    px2  = np.exp(-(h + x)**2)/spi/H
    pH1  = np.exp(-(h - H)**2)/spi/H
    pH2  = np.exp(-(h + H)**2)/spi/H
    fx1  = sp.erf(h - x)
    fx2  = sp.erf(h + x)
    fH1  = sp.erf(h - H)
    fH2  = sp.erf(h + H)
    rfh  = sp.erf(h)
    xph  = np.exp(-h**2)/spi/H
    
    d   =  fH2+fH1-2*rfh
    f   =  fx1+fx2-(1-h/H)*fH1-(1+h/H)*fH2 + pH1-pH2
    fw  =  -2*H*((h - x)*px1+(h + x)*px2)+pH1-pH2
    fh  =  (2*H*(px1+px2)+(fH1-fH2)/H)
    dw  =  2*H*((H-h)*pH1-(H+h)*pH2 + 2*h*xph)
    dh  =  2*H*(pH1+pH2 -2*xph)
    g   =  f/d
    gw  =  (fw/d-f*dw/d**2)/w
    gh  =  (fh/d-f*dh/d**2)/w
    return(g,gh,gw)


def fit_serfzL(x,y,I,h,w,H):
    # [L a Lh Lw Lhh Lww Lhw] = fit_serfzL(x,y,h,w,H) return residual and derivatives

    (f,fh,fw)=serfz(x,h,w,H)
    
    yy  = (I*y*y).sum()
    yf  = (I*y*f).sum()
    ff  = (I*f*f).sum()
    yfh = (I*y*fh).sum()
    ffh = (I*f*fh).sum()
    yfw = (I*y*fw).sum()
    ffw = (I*f*fw).sum()
    
    fhfh = (I*fh*fh).sum()
    fhfw = (I*fh*fw).sum()
    fwfw = (I*fw*fw).sum()
    
    a=yf/ff
    L= (yy-yf**2/ff)/2
    
    Lh = yf/ff**2*(yf*ffh-ff*yfh)
    Lw = yf/ff**2*(yf*ffw-ff*yfw)
    
    Lhh = (-ff**2*yfh**2 + 4*ff*ffh*yf*yfh + ff*fhfh*yf**2 - 4*ffh**2*yf**2)/ff**3
    Lww = (-ff**2*yfw**2 + 4*ff*ffw*yf*yfw + ff*fwfw*yf**2 - 4*ffw**2*yf**2)/ff**3
    Lhw = (-yfw*yfh*ff**2 + yf*(fhfw*yf + 2*ffh*yfw + 2*ffw*yfh)*ff - 4*yf**2*ffh*ffw)/ff**3
    return(L,a,Lh,Lw,Lhh,Lww,Lhw)


# Fit zero average symmetric error function using LM
def fit_serfz(x,y,I,h=None,w=None,H=None):
    # x coordinate
    # y function to fit
    # I integration weights sum(I)=H
    # h initial guess for h
    # w initial guess for w

    if h is None: h=H/2
    if w is None: w=H/10

    maxj=100
    maxk=100

    tole=1e-10
    tolp=1e-10*H
    (L,a,Lh,Lw,Lhh,Lww,Lhw) = fit_serfzL(x,y,I,h,w,H)
    #print('L {:}:a {:}:'.format(L,a))
 

    bh=0
    bw=0
    it=1
    for j in range(maxj):
        Lo=L
        e = np.sqrt(Lh**2+Lw**2)
        for k in range(maxk):
            Jhh = Lhh+bh
            Jww = Lww+bw
            D=Jhh*Jww - Lhw*Lhw
            dh=(Lhw*Lw-Jww*Lh)/D
            dw=(Lhw*Lh-Jhh*Lw)/D
            qh=h+dh
            qw=min(w*10,max(w/10,w+dw))
            it=it+1
            (L,a,Lh,Lw,Lhh,Lww,Lhw) = fit_serfzL(x,y,I,qh,qw,H)
            ndp=np.sqrt(dh**2+dw**2)

            if L<Lo or ndp<=tolp: 
                bh=bh/2
                bw=bw/2
                break
            
            if bh==0 or bw==0:
                bh=abs(Lhh)
                bw=abs(Lww)
            else:
                bh=bh*2
                bw=bw*2

        h=qh
        w=qw
        if e<tole or ndp<tolp: 
            break
    Y=a*serfz(x,h,w,H)[0]
    return(h,w,a,Y)



# Fit general symmetric error function using LM
# Works by removing the average value and then fitting the zero mean function
def fit_serf(x,y,I,h=None,w=None,H=None):
    # x coordinate
    # y function to fit
    # I integration weights sum(I)=H
    # h initial guess for h
    # w initial guess for w
    # H system height

    b = (I*y).sum()/H

    (h,w,a,Y)=fit_serfz(x,y-b,I,h,w,H)
    return(h,w,a,b,b+Y)

def serf11(x,h,w): #even erf function such that serf1(0)=1
    X=x/w
    H=h/w
    hh=H**2
    xx=X**2
    ww=w*w
    e = 2*np.exp(-hh-xx)/(spi*w)
    if h<w*1e-5:
        hhhh=hh*hh
        xxxx=xx*xx
        y = ((45 + (6*xxxx - 8*xx)*hhhh + 30*xx*hh)*exp(-xx))/45
        I = spi*w*(hhhh + 30*hh + 90)/180
        yw = 4*(hhhh*xxxx + (-16/3*hhhh + 5*hh)*xx + 4*hhhh - 10*hh + 15/2)*xx*exp(-xx)/(15*w)
        yww = 8*(hhhh*xx**3 + (-65/6*hhhh + 5*hh)*xxxx + (28*hhhh - 55/2*hh + 15/2)*xx - 14*hhhh + 25*hh - 45/4)*xx*exp(-xx)/(15*ww)
    else:
        eh=sp.erf(H);
        e0= 2*H*np.exp(     -H**2)/(w*spi*eh);
        e1=(H-X)*np.exp(-(H-X)**2)/(w*spi*eh);
        e2=(H+X)*np.exp(-(H+X)**2)/(w*spi*eh);
        y  = (sp.erf(H+X)+sp.erf(H-X))/(2*eh);
        yw  = e0*y-e1-e2 ;
        yww = - 2*e0*(e1+e2) - 2*( ((H-X)**2-1)*e1 + ((H+X)**2-1)*e2)/w+ 2*e0*((hh - 1)/w + e0)*y;
        I = h/eh      
    return(y,I,yw)

def fit_serf1_w(x,y,W,h,w=None,B=None):
    #[w b]=fit_serf_w(x,y,W,h,w) % Fit a symmetric error function with fixed h erf1(0)=1 erf1(x)->0 
    #        erf((h - x)/w) + erf((h + x)/w)
    # y =  B ------------------------------
    #             2*erf(h/w)
    # x coordinate
    # W is integration weight
    # w is initial starting guess or None
    # B is function scaling or zero
    if B is None: B=1; 
    if w is None: w = h/2
    maxj=10;
    n=np.inf
    dw=0
    for j in range(maxj):
        on=n
        ow=w
        for k in range(maxj):
            w=ow+dw;
            (f,I,fw)=erf1(x,h,w);
            e     = f-y/B;
            n=(W*e**2).sum();
            if n<on:  break;
            dw=dw/2;
        Iew   = (W*fw  *  e).sum()
        Iww   = (W*fw  * fw).sum()
        dw    = - Iew/Iww;
        if abs(dw)<w*1e-10:  break
        dw = min(w/2,max(-w/2,dw));
    (Y,I,fw)=erf1(x,h,w);
    return (B*Y,w)
 

# basic serf function  f=erf((h - x)/w)/2 + erf((h + x)/w)/2



#def serfI(xx,hh,w,HH,I=None): # even erf function such that int(serf(x)dx) = 1
#                           # Care is taken so that it is well defined in the limit that h=0
#    h=hh/w
#    H=HH/w
#    x=xx/w
#    
#    spi=np.sqrt(math.pi)
#    if h < 1e-5:
#        eH=np.exp(-H)/sp.erf(H)/spi
#        ex=np.exp(-x)/sp.erf(H)/spi
#        y=2/w*ex*(1 + ((2*eH*H)/3 + (2*x**2)/3 - 1/3)*h**2 + ((4*H**2*eH**2)/9 + ((12*H**3 + 40*H*x**2 - 38*H)*eH)/90 + (2*x**4)/15 - (2*x**2)/5 + 1/10)*h**4)
#    else:
#        if I is None: I = w*((H+h)*erf(H + h)-(H-h)*erf(H-h) - exp(-(H-h)**2)/spi + exp(-(H + h)**2)/spi)
#        y=(erf(x+h)+erf(h-x))/I
#        xn=h-x
#        xp=h+x
#        en=exp(-xn**2)/spi
#        ep=exp(-xn**2)/spi
#        fw = 2*en*xn - 2*ep*xp
#        fh = 2*ep + 2*en
#        fww = -4*xp**3*ep/w + 4*ep*xp/w + 4*xn**3*en/w - 4*en*xn/w
#        fhh = -4*ep*xp/w + 4*en*xn/w
#        fhw = 4*xp**2*ep/w - 2*ep/w + 4*xn**2*en/w - 2*en/w
#        
#        
#        
#        expH=np.exp(-H)
#        
#        
#    I=w/spi*((H-h)*sp.erf(H-h) + (H + h)*erf(H + h) - (exp(-(H-h)**2) - exp(-(H + h)**2))
#             
#
#    
#    a=dict()
#    if h<w*1e-5:
#        H = h/w
#        X = x/w
#        e = 2*np.exp(-H**2-X**2)/(spi*w)
#        y = e*(1+(2*X**2*(1/3)+2/3)*H**2+(2/15*(X**4)+4/15*(X**2)+4/15)*H**4)
#        a['fw']  = e*(2*X**2-1+(4/3*(X**4)-2*X**2)*H**2+(4/15*(X**6)-2/3*(X**4))*H**4)/w
#        a['fh']  = e*((4*X**2*(1/3)-2/3)*H+(8/15*(X**4)-4/15*(X**2)-4/15)*H**3)/w
#        a['fww'] = e*( 4*X**4-10*X**2+2+(8/3*(X**6)-44/3*(X**4)+16*X**2-2)*H**2+(8/15*(X**8)-68/15*(X**6)+28/3*(X**4)-4*X**2)*H**4)/w**2
#        a['fhh'] = e*( 4*X**2*(1/3)-2/3+(8/5*(X**4)-52/15*(X**2)+8/15)*H**2)/w**2
#        a['fhw'] = e*( (-8*X**2+2+8/3*(X**4))*H+16*X**2*(X**4-5*X**2+15/4)*H**3*(1/15))/w**2
#    else: 
#        y = (sp.erf((x+h)/w)-sp.erf((x-h)/w))/(2*h)
#        e = np.exp(-(h**2+x**2)/w**2)
#        c = e*np.cosh(2*x*h/w**2)
#        s = e*np.sinh(2*x*h/w**2)
#        a['fw']  = -(2*(c*h-s*x))/(spi*h*w**2)
#        a['fh']  = 2*c/(spi*h*w)-y/h
#        a['fww'] = -(4*(c*h**3-3*h**2*s*x-c*(w**2-3*x**2)*h+s*w**2*x-s*x**3))/(spi*w**5*h)
#        a['fhh'] = 2*y/h**2-(4*(c*h**2+c*w**2-h*s*x))/(spi*h**2*w**3)
#        a['fhw'] = (4*(c*h**3+c*h*x**2-2*h**2*s*x-(1/2)*s*w**2*x))/(spi*w**4*h**2)
#    return(y,a)

def serf(x,h,w,H): #even erf function such that int(serf(x)dx) = 1
                 # Care is taken so that it is well defined in the limit that h=0
    spi=np.sqrt(math.pi)
    a=dict()
    if h<w*1e-5:
        H = h/w
        X = x/w
        e = 2*np.exp(-H**2-X**2)/(spi*w)
        y = e*(1+(2*X**2*(1/3)+2/3)*H**2+(2/15*(X**4)+4/15*(X**2)+4/15)*H**4)
        a['fw']  = e*(2*X**2-1+(4/3*(X**4)-2*X**2)*H**2+(4/15*(X**6)-2/3*(X**4))*H**4)/w
        a['fh']  = e*((4*X**2*(1/3)-2/3)*H+(8/15*(X**4)-4/15*(X**2)-4/15)*H**3)/w
        a['fww'] = e*( 4*X**4-10*X**2+2+(8/3*(X**6)-44/3*(X**4)+16*X**2-2)*H**2+(8/15*(X**8)-68/15*(X**6)+28/3*(X**4)-4*X**2)*H**4)/w**2
        a['fhh'] = e*( 4*X**2*(1/3)-2/3+(8/5*(X**4)-52/15*(X**2)+8/15)*H**2)/w**2
        a['fhw'] = e*( (-8*X**2+2+8/3*(X**4))*H+16*X**2*(X**4-5*X**2+15/4)*H**3*(1/15))/w**2
    else: 
        y = (sp.erf((x+h)/w)-sp.erf((x-h)/w))/(2*h)
        e = np.exp(-(h**2+x**2)/w**2)
        c = e*np.cosh(2*x*h/w**2)
        s = e*np.sinh(2*x*h/w**2)
        a['fw']  = -(2*(c*h-s*x))/(spi*h*w**2)
        a['fh']  = 2*c/(spi*h*w)-y/h
        a['fww'] = -(4*(c*h**3-3*h**2*s*x-c*(w**2-3*x**2)*h+s*w**2*x-s*x**3))/(spi*w**5*h)
        a['fhh'] = 2*y/h**2-(4*(c*h**2+c*w**2-h*s*x))/(spi*h**2*w**3)
        a['fhw'] = (4*(c*h**3+c*h*x**2-2*h**2*s*x-(1/2)*s*w**2*x))/(spi*w**4*h**2)
    return(y,a)

# Find index of nearest value in an array using bisection
def serf_hwu(x,y,WW,h,w,U1):
   W=WW/WW.sum()
   f=serf(x,h,w)[0];
   F=(W*f).sum()
   f     = f-F
   my=(W*y).sum();
   ny=my+(U1-my)/f[0]*f
   return(ny)



def fit_serf_w(x,y,WW,h,w=None):
    #[w b my]=fit_serf_w(x,y,W,h,w) % Fit a symmetric error function with fixed h
    # Directly solves for my and b and then uses Newton-Raphson iteration to find w
    # The fitted function is
    # f = my + b*(serf(x,h,w)-1/H);
    #
    #
    
    W=WW/WW.sum()
    my=(W*y).sum();
    y=y-my;
    if w is None: w = (x.max()-x.min())/20
    maxj=10;
    for j in range(maxj):
        (f,a)=serf(x,h,w);
        F=(W*f).sum()
        f     = f-F
        fw    = a['fw']  - (W*a['fw']).sum()
        fww   = a['fww'] - (W*a['fww']).sum()
        #Iyy  = (W*y   *  y).sum()
        Ify   = (W*f   *  y).sum()
        Iff   = (W*f   *  f).sum()
        Ifyw  = (W*fw  *  y).sum()
        Iffw  = (W*fw  *  f).sum()
        Ifyww = (W*fww *  y).sum()
        Iffww = (W*fww *  f).sum()
        Ifwfw = (W*fw  * fw).sum()
        b   = Ify/Iff
        #b = np.sqrt(Iyy/Iff)# Match u**2   
        if j==maxj: break
        #dw = Ifyw/Ifwfw/b
        dw=-Ify*Iff*(Iff*Ifyw-Iffw*Ify)/((Ify*Ifyww+Ifyw**2)*Iff**2-(Iffww*Ify+Ifwfw*Ify+4*Iffw*Ifyw)*Ify*Iff+4*Iffw**2*Ify**2)
        #dw=-Ify*Iff*(Iff*Ifyw-Iffw*Ify)/(           Ifyw**2 *Iff**2-(          Ifwfw*Ify+4*Iffw*Ifyw)*Ify*Iff+4*Iffw**2*Ify**2)
        #dw=-        (Iff*Ifyw-Iffw*Ify)/(                                     -Ifwfw*Ify+  Iffw*Ifyw)
        if np.abs(dw)<w*1e-10: break
        w=w+dw
        #logger.info('j {:} w {:} dw {:}'.format(j,w,dw))
    ny=my+b*f
    f=serf(0,h,w)[0]-F
    U1=my+b*f
    
    return (w,U1,ny)
