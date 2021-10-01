function  [fs c fx]=ded_gc_test_f7_2(nm)
%nm='gc/test/f7/14';
c=ded_coord(nm);
p=ded_read_param(nm);
fx=ded_read_hdf(['~/' nm '/force/fx.hdf5']);
fb=ded_read_hdf(['~/' nm '/force/fb.hdf5']);
fd=ded_read_hdf(['~/' nm '/force/fd.hdf5']);
wb=ded_read_hdf(['~/' nm '/force/wb.hdf5']);
av=ded_read_g(nm,'avar');
sv=ded_read_g(nm,'sev');
fns=ded_get_fn(nm,'final',[],'state');fs=ded_read_hdf(fns{end});
if ndims(fs.b)==3
  fnm=fieldnames(fs)
  for j=1:length(fnm)
    if ndims(fs.(fnm{j}))==3
      fs.(fnm{j})=squeeze(mean(fs.(fnm{j}),2));
    end
  end
end



figure(1);clf;

subplot(2,1,1);hold('on');h=[];l={};
if isfield(fx,'wdl'); h(end+1) = plot(c.Ax(1+fx.wdls),fx.wdl/max(fx.wdl),'r-s');l{end+1}='wdl';end;
if isfield(fx,'wul'); h(end+1) = plot(c.Ax(1+fx.wuls),fx.wul/max(fx.wul),'b-^');l{end+1}='wul';end;
if isfield(fx,'ful'); h(end+1) = plot(c.Ax(1+fx.fuls),fx.ful/max(fx.ful),'g->');l{end+1}='ful';end;
if isfield(fx,'sdl'); h(end+1) = plot(c.Ax(1+fx.sdls),fx.sdl/max(fx.sdl),'k-v');l{end+1}='sdl';end;
axis('tight');a1=axis;
legend(h,l,'location','best');

subplot(2,1,2);hold('on');h=[];l={};
if isfield(fx,'wbl'); h(end+1) = plot(c.Ax(1+fx.wbls),fx.wbl/max(fx.wbl),'r-s');l{end+1}='wbl';end;
if isfield(fx,'fbl'); h(end+1) = plot(c.Ax(1+fx.fbls),fx.fbl/max(fx.fbl),'b-v');l{end+1}='fbl';end;
if isfield(fx,'sbl'); h(end+1) = plot(c.Ax(1+fx.sbls),fx.sbl/max(fx.sbl),'k-v');l{end+1}='sbl';end;
axis('tight');a2=axis;
legend(h,l,'location','best');
subplot(2,1,1);set(gca,'xlim',[0 max(a1(2),a2(2))]);
subplot(2,1,2);set(gca,'xlim',[0 max(a1(2),a2(2))]);


if ~isfield(fs,'divz')
  fs.divz=fx.GUR;
end

if isfield(av,'ww') | isfield(av,'dwdz') | isfield(sv,'divz')
  figure(2);clf;
  if isfield(av,'ww');subplot(3,1,1);plot(c.Az,sqrt(av.ww));ylabel('ww');end
  if isfield(av,'dwdz');subplot(3,1,2);plot(c.Az,av.dwdz);ylabel('dwdz');end
  if isfield(sv,'divz');subplot(3,1,3);plot(c.Az,sv.divz);ylabel('divz');end
end
  

ww=ichebintw(c.NAz);
w=ichebintw(c.Nz);
figure(7);clf
fs.divzc=ichebf2c(fs.divz,1);
d=mean(ichebc2f(fs.divzc(1:c.Nz),1),2);

if isfield(fs,'db');
  fs.dbc=ichebf2c(fs.db,1);
  fs.dbb=mean(ichebc2f(fs.dbc(1:c.Nz),1),2);
  b=fs.dbb;
else
  b=max(fs.b(:,fx.sbls(end)),[],2);
end

[pp dd]=fit_erf(c.z,d);disp(sprintf('U=[%6.2f %6.3f], hu=%6.3f wu=%6.4f',pp(1)-pp(2),pp(1)+pp(2),pp(4),1/pp(3)));
[pp bb]=fit_erf(c.z,b);disp(sprintf('B=[%6.2f %6.3f], hb=%6.3f wb=%6.4f',pp(1)-pp(2),pp(1)+pp(2),pp(4),1/pp(3)));
disp(sprintf('B=%6.2f, hu=%6.3f Qu=%6.4f Qb=%6.4f %6.4f',w*b,w*(d>0),w*max(0,d),w*(b.*d),w*max(0,b.*d)));
subplot(2,1,1);plot(c.z,dd,c.z,d);ylabel('ddiv erf');
subplot(2,1,2);plot(c.z,bb,c.z,fs.b(:,fx.sbls+1)); ylabel('b erf');



%dc;ded_gc_2db(nm,{'b','psi'},'y')

f=findmax(max(abs(diff(fs.w,1,2)),[],2));

figure(3);clf;
subplot(4,1,1);plot(c.Jx(1:20),fs.w(f,1:20),'-s');ylabel('w');xlabel('x');
subplot(4,1,2);plot(c.Jz,fs.w(:,1:10));ylabel('w');xlabel('z');
subplot(4,1,3);plot(c.Jx(1:20),fs.p(f,1:20),'-s');ylabel('p');xlabel('x');
subplot(4,1,4);plot(c.Jz,fs.p(:,1:10));ylabel('p');xlabel('z');


P=fs.p+fs.u.^2/2;

rg=fx.wdls+1;
d=fd.fd;

fs.divzc=ichebf2c(fs.divz,1);
divz=ichebc2f(fs.divzc(1:c.Nz),1);
d(:,fx.wdls+1)=ichebc2f(fs.divzc(1:c.Nz),1).*fx.wdl';

figure(5);clf;

phs=fs.p+p.g*ichebintf(fs.b,c.dimz,[],p.H,1,true);phs=phs-phs(end,:);
%subplot(2,1,1);plot(c.Jz,phs(:,rg));xlabel('z');ylabel('phs');
%subplot(2,1,2);plot(c.Jx(rg),fs.p(:,rg(1))-fs.p(:,rg));xlabel('x');ylabel('p');
subplot(2,1,1);plot(c.Jz,phs);xlabel('z');ylabel('phs');
subplot(2,1,2);plot(c.Jx,phs);xlabel('x');ylabel('phs');
%subplot(2,1,2);plot(c.Jx,fs.p(1,:)-fs.p);xlabel('x');ylabel('p');


if false
  k=pi/(c.dx*p.ddivj*p.ddivk)
  sx=min(pi/2,c.x*k);
  wdl=cos(sx).^p.ddivj;
  S=sum(wdl*c.dAx);
  wdl=wdl/S;
  idl=(8/3+cos(sx).^4+4*cos(sx).^2*(1/3)).*sin(sx)/(5*k)/S;
  
  plot(c.x(rg),fs.u(:,rg)-divz.*idl(rg)');
  
  plot(c.x(ix),fx.wul.*fx.wdl,'s-',c.x(ix),fx.ful,'s-');
end
fx.UL=double(fx.UL)

Pux=fs.u_parity(1);
dimx=c.dimx;

if true
  wudl  = zeros(1,c.NAx);wudr  = zeros(1,c.NAx);
  ful   = zeros(1,c.NAx);fur   = zeros(1,c.NAx);
  wdl   = zeros(1,c.NAx);wdr   = zeros(1,c.NAx);
  dwdl  = zeros(1,c.NAx);dwdr  = zeros(1,c.NAx);
  
  if isfield(fx,'wudl'); wudl( 1+fx.wdls)  = fx.wudl; end
  if isfield(fx,'dwdl'); dwdl(1+fx.wdls)   = fx.dwdl; end
  if isfield(fx, 'wdl');  wdl( 1+fx.wdls)  = fx.wdl;  end
  if isfield(fx, 'ful');  ful( 1+fx.fuls)  = fx.ful;  end
  
  if isfield(fx,'wudr'); wudr( 1+fx.wdrs)  = fx.wudr; end
  if isfield(fx,'dwdr'); dwdr(1+fx.wdrs)   = fx.dwdr; end
  if isfield(fx, 'wdr');  wdr( 1+fx.wdrs)  = fx.wdr;  end
  if isfield(fx, 'fur');  fur( 1+fx.furs)  = fx.fur;  end
  
  %keyboard
  %dful=pr_diff(ful,c.dAx,2,1,[],Pux);
  %plot(c.Ax,max(0,dful-wdl));
   
   
  DU=p.V*(divz-fx.UL);
  UL=p.V*fx.UL+DU*ful;
  UR=p.V*fur;
  fL=UL.^2/2-DU*wdl/p.Re;
  FL=(wdl*p.V*fx.UL + wudl.*DU  - dwdl/p.Re).*DU;
  FR=(                wudr.*p.V - dwdr/p.Re).*p.V;

  if false
    % check the calculation of the body force
    rg=1:fx.wdls(end)+5;
    x=c.Ax(rg);
    dfL=pr_diff(fL,c.dAx,c.dimx,1,[],Pux^2);
    rfu=ded_read_hdf('force/rfu.hdf5');
    plot(x,dfL(:,rg)-rfu.rfu(:,rg),x,FL(:,rg)-rfu.rfu(:,rg));
    
    
    
    %plot(IFL(1,fx.wdls)-fL(1,fx.wdls));
    %F=zeros(c.Nz,c.NAx);
    %IF=zeros(c.Nz,c.NAx);
    
    IFL=pr_int(FL,c.dAx,c.dimx,[],Pux);
    IFR=pr_int(FR,c.dAx,c.dimx,[],Pux);
    
  IFL=fftresample(IFL,[c.Nz c.Nx], [0 Pux]);
  IFR=fftresample(IFR,[c.Nz c.Nx], [0 Pux]);
  ful=fftresample(ful,[c.Nz c.Nx], [0 Pux]);
  wdl=fftresample(wdl,[c.Nz c.Nx], [0 Pux]);
  end
  figure(6);clf;
  rg=1:fx.wdls(end)+1;
  UL=fftresample(UL,[c.Nz c.Nx], [0 Pux]);
  subplot(2,1,1);plot(c.x(rg),fs.u(:,rg)-UL(:,rg));xlabel('x');ylabel('du');
  subplot(2,1,2);plot(c.z,    fs.u(:,rg)-UL(:,rg));xlabel('z');ylabel('du');
  figure(9);clf;
  subplot(2,1,1);plot(c.x(rg),fs.w(:,rg));xlabel('x');ylabel('dw');
  subplot(2,1,2);plot(c.z,    fs.w(:,rg));xlabel('z');ylabel('dw');

  %  rg=find(IFL<0.5-3e-5);
  
% $$$   figure(6);clf;
% $$$   subplot(2,1,2);plot(c.z,fs.u(:,rg)./ful(rg)',c.z,divz);ylabel('u');
% $$$   xlabel('z');
% $$$   
% $$$   figure(8);clf;
% $$$   
  if false
  
  if isfield(fx,'ful');
    Iful=fx.ful'.^2/2;
  else;
    Iful=0;
  end
  if isfield(fx,'fur');
    Ifur=fx.fur'.^2/2;
    Ifur=Ifur(1)-Ifur;
  else
    Ifur=0;
  end

  IF=zeros(c.Nz,c.NAx);
  IF(:,1+fx.wdls)=divz.^2*Iful-divz*fx.wdl'/p.Re;
  IF(:,1+fx.wdls(end):end)=repmat(IF(:,1+fx.wdls(end)),1,c.NAx-fx.wdls(end));
  if isfield(fx,'wdr');
    wdr=fx.wdr(1)-fx.wdr';
    IF(:,1+fx.wdrs)=IF(:,1+fx.wdls(end))-repmat(+p.V^2*Ifur-p.V*wdr/p.Re,c.Nz,1);
  end
  DF=pr_diff(IF,c.dAx,2,1,[],1);
  plot(DF(1,:)-F(1,:));
  
  IF=fftresample(IF,[c.Nz c.Nx],[0 1]);
  %IF=pr_int(F,c.dx,2,1,-1); % parity seems wrong
  P=fs.u.^2/2+fs.p-IF;
  plot(c.x,P);
  end
  end
return;

L=5;
N=100;
x=(0.5:N-0.5)*L/N;
k=1;
dx=x(2)-x(1);
f=sin(k*x/L*pi);
F1=pr_int(f,dx,2,1,-1); % parity seems wrong
F2=-L/(k*pi)*(cos(k*x/L*pi)-1);
plot(x,F1,x,F2);

[fs c fx]=ded_gc_test_f7_2('gc/01');
d=ded_read_hdf('force/fd.hdf5');d=d.fd;
%Pux=0;
%
%d=ded_resample(fd.fd,1/p.AA,Pux);
fu=ded_read_hdf('force/fu.hdf5');

W=ded_read_hdf('force/W.hdf5')
subplot(2,1,1);plot(c.Ax,fu.fu(end,:),c.x,fs.u(end,:),c.x,W.W(end,:));
subplot(2,1,2);plot(c.Ax,fu.fu(  1,:),c.x,fs.u(  1,:),c.x,W.W(1,:));


[fs c fx]=ded_gc_test_f7_2('gc/01');
d=ded_read_hdf('force/fd.hdf5');d=d.fd;
f=ded_read_hdf('force/rfu.hdf5');f=f.rfu;
w=ichebintw(size(d,1));
ww=ichebintw(size(f,1));

Idu=w*(d.*fs.u);
If = ww*f;
plot(c.x,Idu,'-s',c.Ax,If,'^-');
title([sum(c.dx*Idu) sum(c.dAx*If)]);

f=zeros(c.NAx,1);f(fx.sdls+1)=fx.sdl;
sdl=fftresample(f,c.Nx,Pux^2);
e=[sum(fx.sdl*c.dAx) sum(sdl*c.dx)];
   
 


plot(c.x,fs.p(end,:),c.Ax,a.rfu(end,:))

wul = zeros(c.NAx,1);
ful = zeros(c.NAx,1);
wul(1+fx.wuls)= fx.wul;
ful(1+fx.fuls)= fx.ful;
wul=fftresample(wul,c.Nx,Pux^2)';
ful=fftresample(ful,c.Nx,Pux)';
u=fx.UL+ful.*(divz-fx.UL);

imagesc(wul.*(fs.u-u));
F=w*(wul.*(u-fs.u));
plot(c.x,F);
title(sum(c.dx*F));

plot(c.x,w*fs.u,c.x,w*u);

