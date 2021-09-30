def gcf7(param,problem)


         ddiv     = dp.get_bool_param(param,'ddiv')
         ddivj    = dp.get_param(param,'ddivj')
         ddivk    = dp.get_param(param,'ddivk')
         dbj      = dp.get_param(param,'dbj')
         dbk      = dp.get_param(param,'dbk')


               # Forces only by adding and removing fluid Tx should be sincos

        nfd = domain.new_field()  # div(u)
        ju.set_parity(nfd,Tx,Ty,Tz,1,1,1)
        nfd.set_scales(AA)

        if Tx=='SinCos':
            UL=0 # velocity is zero on left of domain
            if ddivj is None:
                mm=2
                (wdr,wdrs) = FFTAx.dfsl(xx-L,mm)           # Forcing for right divergence
                (wdl,wdls) = FFTAx.dfsl(xx,mm)             # Forcing for left divergence
                (wdl,wdls) = FFTAx.dfsl(xx,mm)             # Sampling region for left divergence
                wur        = FFTAx.f(   xx-L,mm)[wdrs,...] # Forcing for right velocity
                wul        = FFTAx.f(   xx,mm)[wdls,...]   # Forcing for left velocity
                wdl=wdl.reshape((wdl.shape[0],1))
                sdl=sdl/(wdl.sum()*dxx)
                sdls=wdls
            else:
                ju.logger.info('ddivj {} ddivk {}'.format(ddivj,ddivk))
                (wdls,wdl,wul,Nl,dwdl)=FFTAx.deltajks(xx.flatten(),  ddivj,0,ddivk)
                (wdrs,wdr,wur,Nr,dwdr)=FFTAx.deltajks(xx.flatten()-L,ddivj,0,ddivk)
                Nwdl = (wdl.sum()*dxx)
                Nwdr = (wdr.sum()*dxx)
                wdl  = wdl/Nwdl
                wdr  = wdr/Nwdr
                dwdl = dwdl/Nwdl
                dwdr = dwdr/Nwdr
                wdl  = wdl.reshape((wdl.shape[0],)+xx.shape[1:])
                wdr  = wdr.reshape((wdr.shape[0],)+xx.shape[1:])
                dwdl = dwdl.reshape((dwdl.shape[0],)+xx.shape[1:])
                dwdr = dwdr.reshape((dwdr.shape[0],)+xx.shape[1:])
                wul  = wul.reshape((wul.shape[0],)+xx.shape[1:])
                wur  = wur.reshape((wur.shape[0],)+xx.shape[1:])
                sdl  = wdl/(wdl.sum()*dxx)
                sdls = wdls
 
            if dbj is None:
                (wbl,wbls) = FFTAx.dfsl(xx,4)       # Forcing for left buoyancy
                (wbr,wbrs) = FFTAx.dfsl(L-xx,4)     # Forcing for right buoyancy
                (sbl,sbls) = FFTAx.dfsl(xx,4)
                sbl=sbl.reshape((sbl.shape[0],1))/(dxx*sbl.sum())
                fbx=FFTAx.heaviside(L/2-xx,2)[0]
            else:
                (wbls,wbl)=FFTAx.deltajks(xx.flatten()  ,dbj,0,dbk)[0:2]
                (wbrs,wbr)=FFTAx.deltajks(xx.flatten()-L,dbj,0,dbk)[0:2]
                wbl  = wbl/wbl.max()
                wbr  = wbr/wbr.max()
                wbl  = wbl.reshape((wbl.shape[0],)+xx.shape[1:])
                wbr  = wbr.reshape((wbr.shape[0],)+xx.shape[1:])
                sbl  = wbl/(wbl.sum()*dxx)
                sbls = wbls
                fbx  = FFTAx.heaviside(L/2-xx,2)[0]
        else:
            UL=-1 # Velocity is -1 on left of domain
            tu   = gx.new_field(Px=0)
            wul  = FFTAx.Hjk(ddivj,AA*ddivk)
            nwul = wul.shape[0]
            wdls = slice(0,nwul)
            #tu['g'] = wul[np.int64(np.minimum(nwul-1,np.abs(np.mod(np.arange(Nxx)+Nxx/2,Nxx-1)-Nxx/2)))]
            Nxx2=np.int64(Nxx/2)
            tu['g'][0:nwul] = wul
            tu['g'][nwul:Nxx2] = 1
            tu['g'][Nxx2:Nxx2+nwul] = np.flip(wul)
            tu['g'][Nxx2+nwul:] = 0
            
            wdl  = tu.differentiate('x')['g'][wdls]
            dwdl = tu.differentiate('x','x')['g'][wdls]
            wdr  = None
            wur  = None
            dwdr = None
            sdl  = wdl/(wdl.sum()*dxx)
            sdls = wdls
            sdl  = sdl.reshape((sdl.shape[0],)+xx.shape[1:])
            # Buoyancy sampling
            sbl  = wdl/(wdl.sum()*dxx)
            sbl  = sbl/(sbl.sum()*dxx)
            sbls = slice(2*nwul,3*nwul)
            
            # Buoyancy forcing weight
            wbl  = np.concatenate((wul,np.ones_like(wul),np.flip(wul)))
            wbls = slice(0,3*nwul)
            
            wbr = None
               
            # Buoyancy forcing function
            fbx  = np.zeros((Nxx,))
            fbx[     1*nwul:     2*nwul] = wul # Blend between densities
            fbx[     2*nwul:Nxx2+1*nwul] = 1
            fbx[Nxx2+1*nwul:Nxx2+2*nwul] = np.flip(wul)

            fbx  =  fbx.reshape(( fbx.shape[0],)+xx.shape[1:])
            wbl  =  wbl.reshape(( wbl.shape[0],)+xx.shape[1:])
            sbl  =  sbl.reshape(( sbl.shape[0],)+xx.shape[1:])
            wdl  =  wdl.reshape(( wdl.shape[0],)+xx.shape[1:])
            wul  =  wul.reshape(( wul.shape[0],)+xx.shape[1:])
            dwdl = dwdl.reshape((dwdl.shape[0],)+xx.shape[1:])

            #dtuu  = tu.differentiate('x')['g']
            #ddtuu = tu.differentiate('x','x')['g']
            #tuu=tu['g']
            #ju.write_hdf5_vars('force/tuu.hdf5',('tuu','dtuu','ddtuu','nwul'),locals())
            #quit()
            
      
        
        #fuzz = ju.smoothwf(zz,nfz,hu-Wz,hu+Wz)
        #fuzz = U1/U-fuzz/(gIzz*fuzz).sum()*(1+U1/U)*H broken for parallel
        UR = 1-4*ju.smoothwf(zz,nfz,0,H)
        if ddiv:
            sev['divz']   = UR.flatten()
            gdw                = np.zeros((Nzz,), dtype=np.float64)
            gww                = np.zeros((Nzz,), dtype=np.float64)
            avar['dwdz']       = np.zeros((Nzz,), dtype=np.float64)
            avar['ww']          = np.zeros((Nzz,), dtype=np.float64)

        #ju.logger.info('wbr.shape {} wbrs {}'.format(wbr.shape,wbrs))
        #ju.logger.info('wbl.shape {} wbls {}'.format(wbl.shape,wbls))
        nwb.set_scales(AA)
        if wbl is not None: nwb['g'][wbls,...]  = wbl/wbl.max()
        if wbr is not None: nwb['g'][wbrs,...] += wbr/wbr.max()
                    

        if db:
            gbp            = np.zeros((Nzz,), dtype=np.float64)
            sev['db'] = np.zeros((Nzz,), dtype=np.float64)
            avar['db']     = np.zeros((Nzz,), dtype=np.float64)
            if yycolor: yycomm.Gatherv(sendbuf=hbzz, recvbuf=(sev['db'], zzcounts), root=0)
            sbp=sev['db']
            sev['db'] = ju.comm.bcast(sev['db'],0)
            hbzz=((sev['db'][zzslices]).reshape(szzzl))
        else:
            nfb.set_scales(AA)
            nfb['g'] = hbzz*FFTAx.heaviside(L/2-xx,2)[0]


        if pwvl or pwvr or pwwl or pwwr:
            nww = domain.new_field()  # div(u)
            ju.set_parity(nww,Tx,Ty,Tz,1,1,1)
            nww.set_scales(AA)
            nww['g']=0
            if pwul or pwvl or pwwl: nww['g'][wdls,...] += wdl/wdl.max()
            if pwur or pwvr or pwwr: nww['g'][wdrs,...] += wdr/wdr.max()
            nww['g'] = ju.clip01(nww['g'])
            if pwwl or pwwr:
                if pwul or pwur: problem.substitutions['wu'] = 'ww'
                if pwvl or pwvr: problem.substitutions['wv'] = 'ww'
            elif pwul or pwur:
                nwu=nww
                del(nww)
                if pwvl or pwvr: problem.substitutions['wv'] = 'wu'
            else:
                nwv=nww
                del(nww)
                
        if wdl is not None: nfd['g'][wdls,...] = wdl*UR
        if wdr is not None: nfd['g'][wdrs,...] = wdr
            
        if pful or pfur:
            rfu = domain.new_field()  
            ju.set_parity(rfu,Tx,Ty,Tz,-1,1,1)
            rfu.set_scales(AA)

        if pful:
            ful=wul*wdl
            fulv = -dwdl/Re
        if pfur:
            if wur  is not None: fur  = wur*wdr
            if dwdr is not None: furv = -dwdr/Re
        if Tx=='SinCos':
            if pwul: rfu['g'][wdls,...] = ful*sev['divz'][zzslices].reshape(szzzl)**2+fulv*sev['divz'][zzslices].reshape(szzzl)
            if pwur: rfu['g'][wdrs,...] = fur*U**2+furv*U
        else:
            if pwul: rfu['g'][wdls,...] = ((wul*sev['divz'][zzslices].reshape(szzzl)-U)*wdl-dwdl/Re)*sev['divz'][zzslices].reshape(szzzl)

                        UR=sev['divz'][zzslices].reshape(szzzl)
        DU=U*(UR-UL)
        UU=U*UL+wul*DU
        if pful: problem.parameters['rfu']['g'][wdls,...] = dU*(wdl*UU-dwdl/Re)
        if pfur: problem.parameters['rfu']['g'][wdrs,...] = fur*U**2-dwdr/Re*U
 
        # This forcing approach forces the horizontal velocity with a constant value, -U,
        # The Buoyancy must be large enough to overcome this
        # Force the velocity and buoyancy to zero in a region behind the moving current
 
