function ded_pm_test_f7(nm,trg)
if nargin==0
  nm='pm/f7/01';
  nm='pm/test/02';
end
dd=[ded_dedalus_data_dir '/'];

if nargin<2
  trg=[0 inf];
end

c1=[1 1 1];
c2=[0 0 1];
c3=0.7*[1 1 1];

ded_force_fx(nm);
q=ded_read_param(nm);
c=ded_coord(nm);
fx=ded_hdf(nm,'force/fx');
wb=ded_hdf(nm,'force/wb');
fr=ded_hdf(nm,'force/force_s1');
[ayzc ayz]  = ded_read_javrg(nm,'ayz',trg,true');
avar = ded_read_g(nm,'avar',[],trg);
sev  = ded_read_g(nm,'sev', [],trg);
%a    = ded_read_javrg(nm,'a',trg,'combines');
a=[];

x   = c.x;
y   = c.y;
z   = c.z; 
xx  = c.Ax;
yy  = c.Ay;
zz  = c.Az;
dxx = c.dAx;
exx=[xx;xx+q.L];

if ~isempty(avar)
  fnavar=fieldnames(avar);
  mavar={};
  for j=1:length(fnavar)
    sz=size(avar.(fnavar{j}));
    if length(sz)==3
      if all(sz==[c.NAy c.NAz length(avar.t)])
        mavar{end+1}=fnavar{j}; 
      end
    end
  end
end
if ~isempty(sev)
  fnsev=fieldnames(sev);
  msev={};
  for j=1:length(fnsev)
    sz=size(sev.(fnsev{j}));
    if length(sz)==3
      if all(sz==[c.NAy c.NAz length(sev.t)])
        msev{end+1}=fnsev{j}; 
      end
    end
  end
end


ayzd = [];
if ~isempty(ayz)
  if isfield(ayz,'dudx');
    ayzd = cat(2,ayz.dudx);
  end
end
if ~isempty(wb)
  wb=squeeze(min(min(wb.wb)));
  wb=fftresample(wb,c.NAx,1);
end
if ~isempty(fx)
  Bwxx  = fx.Bwxx;
  figure;
  subplot(3,2,1);plot(xx(1+fx.wdls),fx.wul,'s-');ylabel('wul');axis('tight');
  subplot(3,2,2);plot(xx(1+fx.wdls),fx.wdl,'s-');ylabel('wdl');axis('tight');
  subplot(3,2,3);plot(xx(1+fx.wdrs),fx.wur,'s-');ylabel('wur');axis('tight');
  subplot(3,2,4);plot(xx(1+fx.wdrs),fx.wdr,'s-');ylabel('wdr');axis('tight');
  subplot(3,2,5);plot(xx(1+fx.sdrs),fx.sdr,'s-');ylabel('sur');axis('tight');
  subplot(3,2,6);plot(xx(1+fx.sdrs),fx.sdr,'s-');ylabel('sur');axis('tight');
else
  Bwxx=[];
end
fyz=ded_read_hdf([dd nm '/force/fyz.hdf5']);
nwdr=[];
rr=[];
wrr=[];
if ~isempty(fyz)
  if isfield(fyz,'nwdr'); nwdr = squeeze(fyz.nwdr); end
  if isfield(fyz,'rr');   rr   = squeeze(fyz.rr);   end
  if isfield(fyz,'wrr');  wrr  = squeeze(fyz.wrr);  end
end


if ~isempty(wb) & ~isempty(fx)
  figure;clf
  ewb=[wb;flip(wb)];
  if isfield(fx,'wdivx')
    ewdivx=[fx.wdivx;flip(fx.wdivx)];
    rgxmin=min(xx(fx.wdivx>.999));
    rgxmax=max(xx(fx.wdivx>.999));
  else
    rgxmin=min(xx(wb<.999));
    rgxmax=max(xx(wb<.999));
    ewdivx=NaN*exx;
  end
  for j=1:3
    subplot(3,1,j);
    h=plot(exx,ewdivx,'s-',exx,ewb);
    legend(h,'wdivx','wb','location','best')
    if j==1
      ax=[0 q.L 0 1];
      title('wb and wdivx');
    elseif j==2
      if isfield(q,'wdivxl');title(sprintf('wdivxl %4.1f',q.wdivxl));end
      axis([0 rgxmin 0 1]);
    else
      if isfield(q,'wdivxr');title(sprintf('wdivxr %4.1f',q.wdivxr));end
      axis([rgxmax q.L 0 1]);
    end
  end
end

% $$$ if ~isempty(fr)
% $$$   if isfield(fr,'wu')
% $$$     wv=squeeze(mean(mean(fr.wu)));
% $$$   elseif isfield(fr,'wv')
% $$$     wv=squeeze(mean(mean(fr.wv)));
% $$$   elseif isfield(fr,'ww')
% $$$     wv=squeeze(mean(mean(fr.ww)));
% $$$   end
% $$$   wvl=[1*wv(1:q.Nx/2);0*wv(q.Nx/2+1:end)];
% $$$   wvr=[1*wv(1:q.Nx/2);1*wv(q.Nx/2+1:end)];
% $$$ else 
wvl=exp(-((c.x    )/(5*c.dx)).^2);
wvr=exp(-((q.L-c.x)/(5*c.dx)).^2);
wvl=wvl/sum(wvl);
wvr=wvr/sum(wvr);
%end
wvl=reshape(wvl,[1 1 c.Nx])/sum(wvl);
wvr=reshape(wvr,[1 1 c.Nx])/sum(wvr);
z=reshape(c.z,[c.Nz 1]);
y=reshape(c.y,[1 c.Ny]);
r=sqrt(z.^2+y.^2);
ur=[];
br=[];
if ~isempty(a)
  vrl=(sum(wvl.*a.v,3).*y+sum(wvl.*a.w,3).*z)./r;
  vrr=(sum(wvr.*a.v,3).*y+sum(wvr.*a.w,3).*z)./r;
  pl=sum(wvl.*a.p,3);
  pr=sum(wvr.*a.p,3);
  ul=sum(wvl.*a.u,3);
  ur=sum(wvr.*a.u,3);
  bl=sum(wvl.*a.b,3);
  br=sum(wvr.*a.b,3);
  clear('a');
  figure;
  subplot(2,2,1);mesh(c.y,c.z,vrl);xlabel('t');title('vr bottom');
  subplot(2,2,2);mesh(c.y,c.z,vrr);xlabel('t');title('vr top');
  subplot(2,2,3);mesh(c.y,c.z,pl);xlabel('t');title('p bottom');
  subplot(2,2,4);mesh(c.y,c.z,pr);xlabel('t');title('p top');
end

if length(mavar)>0
  figure;
  spsz=ceil(sqrt(length(mavar)));
  spsz(2)=ceil(length(mavar)/spsz(1));
  for j=1:length(mavar)
    f=avar.(mavar{j})(:,:,end);
    subplot(spsz(1),spsz(2),j);
    mesh(c.Ay,c.Az,f);xlabel('t');title(sprintf('avar: %s [%8.4f %8.4f]',mavar{j},min(f(:)),max(f(:))));
    axis('tight');
  end
end

if length(msev)>0
  figure;
  spsz=ceil(sqrt(length(msev)));
  spsz(2)=ceil(length(msev)/spsz(1));
  for j=1:length(msev)
    f=sev.(msev{j})(:,:,end);
    subplot(spsz(1),spsz(2),j);
    mesh(c.Ay,c.Az,f);xlabel('t');title(sprintf('sev: %s [%8.4f %8.4f]',msev{j},min(f(:)),max(f(:))));
    axis('tight');
  end
end

oddg=false;
if ~isempty(Bwxx)
  if isfield(q,'oddg'); oddg=q.oddg ;end
  if oddg
    eBwxx=[Bwxx;flip(Bwxx)];
  else
    eBwxx=[Bwxx;-flip(Bwxx)];
  end
  
  figure;
  subplot(2,1,1);plot(exx,eBwxx,'-s');ylabel('fb');
  subplot(2,1,2);plot(exx,[wb;flip(wb)],'-s');ylabel('wb');
end



figure;
nsp=~isempty(ayzd) +isfield(sev,'divx') + isfield(avar,'dudx')+isfield(avar,'q')+isfield(avar,'u');
nn=0;
if ~isempty(ayzd) 
  nn=nn+1;subplot(nsp,1,nn);
  line(xx,0*xx,'color',c3)
  gu_plotg(c.Jx,ayzd,c1,c2);        
  axis([0 q.L median([0;ayzd(ayzd<=0)/5]) 2*median([0;ayzd(ayzd>=0)]) ]);ylabel('ayzd');  
end
if isfield(sev,'divx')
  nn=nn+1;subplot(nsp,1,nn);
  line(xx,0*xx,'color',c3)
  gu_plotg(xx,sev.divx,c1,c2);        
  axis('tight');ylabel('divx');  
end
if isfield(avar,'dudx')
  nn=nn+1;subplot(nsp,1,nn);
  line(xx,0*xx,'color',c3)
  gu_plotg(xx,avar.dudx,c1,c2);
  axis('tight');ylabel('dudx');;
end
if isfield(avar,'q')
  nn=nn+1;subplot(nsp,1,nn);
  line(xx,0*xx,'color',c3)
  gu_plotg(xx,avar.q,c1,c2);
  axis([0 q.L -inf -0.1*min(avar.q(:))]);ylabel('q');
end
if isfield(avar,'u')
  nn=nn+1;subplot(nsp,1,nn);
  line(xx,0*xx,'color',c3)
  gu_plotg(xx,avar.u,c1,c2); 
  axis('tight');ylabel('u');
  xlabel('x');
end



figure;
dxx=xx(2)-xx(1);
nsp=isfield(sev,'divx')+isfield(avar,'u')+isfield(sev,'qm')+isfield(avar,'topq')+isfield(sev,'topr');   
k=0;
if isfield(sev,'divx') 
  Idivx = sum(sev.divx,1)*dxx; 
  k=k+1;subplot(nsp,1,k);plot(sev.t,Idivx); 
  axis('tight');
  ylabel('Idivx');  
end;
if isfield(avar,'u')    
  Iu    = sum(avar.u,1)*dxx; 
  k=k+1;subplot(nsp,1,k);plot(avar.t,Iu);    
  axis('tight');
  ylabel('Iu'); 
end
if isfield(sev,'qm')    
   k=k+1;subplot(nsp,1,k);plot(sev.t,sev.qm);    
  axis('tight');
  ylabel('qm'); 
end
if isfield(sev,'topr')    
  k=k+1;subplot(nsp,1,k);plot(sev.t,sev.topr);    
  axis('tight');
  ylabel('topr'); 
end
if isfield(avar,'topq')    
  k=k+1;subplot(nsp,1,k);plot(avar.t,avar.topq);    
  axis('tight');
  ylabel('topq'); 
end;

xlabel('t');


fns=ded_get_fn(nm,'divyz');
for j=1:length(fns)
  if j==1;figure;end
  a=ded_read_hdf(fns{j});
  mesh(a.y,a.z,a.divyz);
  xlabel('y');
  xlabel('z');
  drawnow
end




fnf=['~/' nm '/force/rfu-00001.hdf5'];
if isfile(fnf);
  f=ded_read_hdf(['~/' nm '/force/rfu-00001.hdf5']);
  f=squeeze(sum(sum(f.rfu)))*c.dAy*c.dAz;
  If = pr_int(f,c.dAx,1,1,-1);
  f  = fftresample(f,c.Nx,-1);
  If = fftresample(If,c.Nx, 1);
else
  f=0;
  If=0;
end



if ~isempty(a)
  u  = cat(2,ayzc.u);
  d  = cat(2,ayzc.d);
  ud = cat(2,ayzc.ud);
  uu = cat(2,ayzc.uu);
  p  = cat(2,ayzc.p);
  p  = p-p(1,:);
end

Prho=-1;

if oddg
  odd=ded_hdf(nm,'force/oddcx');
  q.g=q.g*fftresample(odd.oddcx,c.NJx,-1);
  Pb=1;
  Ps=1;
else
  Pb=-1;
  Ps=-1;
end
Pu=-1;
Pp=1;
Pd=1;

rho=0;
if ~isempty(a)
  if isfield(a,'b'); rho =rho+q.B*q.g.*cat(2,ayzc.b); end
  if isfield(a,'s'); rho =rho+q.S*q.g.*cat(2,ayzc.s); end
  
  
  Iud    = pr_int( ud,c.dx,1,1,  Pu*Pd);
  duudx  = pr_diff(uu,c.dx,1,1,0,Pu*Pu);
  dpdx   = pr_diff( p,c.dx,1,1,0,Pp);
  ududx  = duudx - ud;
  Iududx = uu    - Iud;
  dudx   = pr_diff( u,c.dx,1,1,0,Pu);
  duddx  = pr_diff( u,c.dx,2,1,0,Pu);
  Irho   = pr_int(rho,c.dx,1,1,  Prho);

  figure;
  subplot(1,5,1);plot(c.x,p);                             xlabel('p')
  subplot(1,5,2);plot(c.x,p - Irho);                      xlabel('-Irho')
  subplot(1,5,3);plot(c.x,p - Irho+ Iududx);              xlabel('+Idudx')
  subplot(1,5,4);plot(c.x,p - Irho+ Iududx-If);           xlabel('-If')
  subplot(1,5,5);plot(c.x,p - Irho+ Iududx-If-dudx/q.Re); xlabel('-nu')
  suptitle('Integrated x momentum balance');
  figure;
  subplot(1,5,1);plot(c.x,dpdx);                           xlabel('p')
  subplot(1,5,2);plot(c.x,dpdx - rho);                     xlabel('-Irho')
  subplot(1,5,3);plot(c.x,dpdx - rho+ ududx);              xlabel('+Idudx')
  subplot(1,5,4);plot(c.x,dpdx - rho+ ududx-f);            xlabel('-If')
  subplot(1,5,5);plot(c.x,dpdx - rho+ ududx-f-duddx/q.Re); xlabel('-nu')
  suptitle('x momentum balance');
end

return;

A="--preset pmf7 --Nx 128 --Re 200 --radius 0.25 -L 4 -W 2 -H 2 --dtjb 0.1 --dtjsev 0.1 --dtjavar 0.1 --topd True --dtforce 100 --pmss None"
B="--reset --pfn test/01 --rfn test/01"

mpiexec -n 8 ded_gc.py $A --Ty Fourier --Tz Fourier --oddg False  -wul False --wur False pm/test/01
mpiexec -n 8 ded_gc.py $B --wur True pm/test/02
mpiexec -n 8 ded_gc.py $B --wur True --wul True pm/test/03



a=ded_read_hdf('~/pm/test/asahi/13/a/a-00010.hdf5');
c=ded_coord('pm/test/asahi/13');
z=reshape(c.z,[c.Nz 1 1]);
y=reshape(c.y,[1 c.Ny 1]);
r=sqrt(y.^2+z.^2);
ur = (ayzc.v.*y+ayzc.w.*z)./max(realmin,r);



return;

nm='pm/test/aino/03';
c=ded_coord(nm);
p=ded_read_param(nm);
dd=[ded_dedalus_data_dir '/'];
fa=ded_read_hdf([dd nm '/force/fyz.hdf5']);

[zz yy]=ndgrid(c.Az,c.Ay);
rr=max(realmin,sqrt(zz.^2+yy.^2));
q=randn(1);

vr=-q./rr;
ny=yy./rr;
nz=zz./rr;
vy=ny.*vr;
vz=nz.*vr;
S1=sum(sum(fayzc.nry.*vy+fayzc.nrz.*vz));
I=sum(fayzc.nwdr(:))*c.dAy*c.dAz;
disp([S1/q I]);





