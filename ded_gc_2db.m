function ded_gc_2db(nm,typs,ft,pp)
%ded_gc_2db('gc/f6/mg/4000/1000/4608/22608',{'b','psi'},'ay');
if nargin<2
  typs=[];
end
if nargin<4
  pp=struct;
end

fns=ded_get_fn(nm,ft);

if isempty(fns)
  return;
end

if ft(1)=='a'
  ddt=true;
else
  ddt=false;
  dt=1;
end

if ~isfield(pp,'clip'); pp.clip=struct; end
if ~isfield(pp,'W');    pp.W=1024;      end
if ~isfield(pp,'ar');   pp.ar=1;        end

clipnm=fieldnames(pp.clip);

if isempty(typs)
  typs={'b'};
end
nmin=nm;
p=ded_read_param(nm); 
c=ded_coord(nm);
p.name=nm;
%b=ded_read_mom(nm);
d=ded_read_g(nm,'yz',[]);

if any(ft=='y')
  scly=1/p.W;
else
  scly=1;
end

a=ded_read_hdf(fns{1});
%  a=ded_read_javrg(nm,'a',trg,true);%ded_read_yavg(nm);
if isfield(a,'t1');  a.t=(a.t1+a.t2)/2; end

ss=ded_read_stats(nm);
ltyps=typs;
if any(strcmp('omegapsi',typs)) | any(strcmp('ppsi',typs)) 
  ltyps=cat(2,ltyps,{'u','v','w'});
end
if any(strcmp('ppsi',typs)) 
  ltyps=cat(2,ltyps,{'p'});
end

dim=find(size(a.b)==length(c.Jz));

U=p.V;
if ~isfield(p,'B')
  p.B=1;
end

tol=1e-3;
switch(p.version)
  case '1.30'
    p.B=1;
    if isfield(d,'u')
      d.u   =   d.u/(p.V);
    end
    d.b   =   d.b/(p.B);
    d.uu  =  d.uu/(p.V);
    d.e   =   NaN;
    d.s   =   NaN;
    a.b   = a.b/p.B;
    a.u   = a.u/p.V;
  case '1.31'
    if ~isempty(g)
      n=min(length(d.t),length(s.t));
      g.t=g.t(1:n);
      dnm=setdiff(fieldnames(g),{'t','x','y','z','nx','ny','nz','nt','sz'});
      for j=1:length(dnm)
        g.(dnm{j})=g.(dnm{j})(:,1:n);
      end   
      Zhbz  = max(0,2*g.bz)./max(tol,d.b);
      Zhbb  =        d.b.^2./max(tol,g.bs);
      d.t=d.t(1:n);
      dnm=setdiff(fieldnames(d),{'t','x','y','z','nx','ny','nz','nt','sz'});
      for j=1:length(dnm)
        d.(dnm{j})=d.(dnm{j})(:,1:n);
      end   
    end
    if ~isempty(d)
      Zhb   =     max(0,d.b)/(p.W*p.B);
      d.ev  =  d.bw./max(tol,d.b);
      d.ev  = d.ev/max(abs(d.ev(:)));
      d.bu  =  d.bu/(p.W*p.hu*p.V1/2);
      d.bu  =  d.bu/max(abs(d.bu(:)));
      %d.u   =   d.u/(p.W*p.hu*p.V1);
      d.uu  =  d.uu/(p.W*p.hu*p.V );
      %    a.b   =   a.b/(p.W*p.B);
      %a.u   =   a.u/(p.W*p.hu*p.V);
      d.E   =   d.E/max(d.E(:));
      d.S   =   d.S/max(d.S(:));
    end
end
if ~isempty(d)
  UU = -mean(d.u/(p.W*p.H));
  [Ut Uf] = unique(d.t);
  UU=UU(Uf);
else
  Ut=a.t;
  UU=repmat(p.V,size(Ut));
end

if strcmp(p.Tx,'SinCos')
  uparity=-1;
  vparity=1;
  wparity=1;
else
  uparity=0;
  vparity=0;
  wparity=0;
end


nz=round(2*p.H/p.L*p.Nx);
at=a.t;
a.H=p.H;
a.x=c.Jx;
a.z=c.Jz;
if ddt;  dt=a.dt; end
scl=scly/dt;

fvel=cell_match_str(typs,{'bpsi','div','omegapsi','ppsi','psi'});
if sum(fvel)>0
  if isfield(a,'u') & isfield(a,'w') 
    [a.psi a.omega a.div]=ded_psi(a.u*scl,a.w*scl,uparity,wparity,p.H,c.dJx);
  else
    typs={typs{find(~fvel)}};
    disp(sprintf('gc_2d.m: No velocity data available',nm));
  end
end

a=ded_zgrid(a,[],[],[],[],[],p.H,1);

a.t=at;
z=a.z;
x=a.x;

S=pp.W/p.L;
if isfield('hu',p)
  HT=min(5*p.hu,p.H);
else
  HT=p.H;
end
T=p.T;


ntyp=length(typs);

k=0;
for j=1:ntyp
  t=typs{j};
  switch(t)
    case {'b','s'}
      k=k+1;
      fnm{k}=t;
      fld{k}=min(1,max(0,a.(t)*scl));
      cnts{k}=(1:2:19)/20;
      ncnts(k)=length(cnts{k});
      fmax(k)=1;
      fmin(k)=0;
      clim{k}=[0 p.B];
    case {'u','w','p'}
      k=k+1;
      fnm{k}=t;
      fld{k}=a.(t)*scl;
      prc=percentile(fld{k}(:),[1 99]);
      fmin(k)=floor(prc(1)*10)/10;
      fmax(k)=ceil(prc(2)*10)/10;
      clim{k}=[fmin(k) fmax(k)];
      if p.V==0
        absu=max(abs([fmin(k) fmax(k)]));
        cnts{k}=linspace(-absu/2,absu/2,6);
      else
        cnts{k}=[0 -0.9*p.V -0.95*p.V -0.99*p.V -0.999*p.V -0.9999*p.V];
      end
      ncnts(k)=length(cnts{k});
    case {'psi','omega','bpsi','phi','omegapsi','ppsi','div'}
      k=k+1;
      fnm{k}=t;
      switch(t)
        case 'bpsi'
          fld{k}=min(1,max(0,a.b*scl/p.B));
          fld2{k}=a.psi;
          ncnts(k)=25;
        case 'psi'
          fld{k}=a.psi;
          ncnts(k)=25;
        case 'omegapsi'
          fld{k}=omega;
          fld2{k}=a.psi;
          ncnts(k)=25;
        case 'ppsi'
          fld{k}=a.p;
          fld2{k}=a.psi;
          ncnts(k)=25;
        case 'omega'
          fld{k}=a.omega;
          ncnts(k)=0;
        case 'div'
          fld{k}=a.div;
          ncnts(k)=0;
      end
      prc=percentile(fld{k}(:),[0 1 99 100]);
      fmin(k)=floor(prc(2)*10)/10;
      fmax(k)=ceil(prc(3)*10)/10;
      clim{k}=prc([1 4]);
      if ncnts(k)>0
        if p.V>0
          cnts{k}='percentile';
        else
          cnts{k}='linear';
        end
      else
        cnts{k}=[];
      end
      cnts{k}=(min(a.psi(:)):p.V/10:max(a.psi(:)))';
    otherwise
      error(sprintf('ded_gc_2d: Unknown type "%s" requested',t));
  end
end
nfld=length(fld);

fh=gcf;
clf;
ah=jsubplot([1 ntyp],[0.05 0.15/ntyp],[0.01 0.1/ntyp],[0.01 0.05/ntyp]);
set(fh,'units','pixels');
pf=get(fh,'position');
set(fh,'position',[pf(1:2),round(S*p.L),ntyp*round(pp.ar*S*HT)]);

xp=([1:ceil(p.L)-1])';


vs=0.4;vt=vs;
bs=0.4;
NULLxz=repmat(NaN,[nz,p.Nx]);
NULLx =repmat(NaN,[p.Nx,1]);

if isempty(a)
  return;
end
for k=1:nfld
  axes(ah(k));
  ih(k)=imagesc(c.Jx,z,NULLxz,[ 0   1  ]); 
  set(gca,'clim',clim{k},'ydir','normal','xlim',[0 p.L],'ylim',[0 HT],'ytick',[],'xtick',[]);
  set(gca,'dataaspectratio',[pp.ar 1 1]);
  ap=get(gca,'position');
  icb(k)=colorbar('position',[0.95 ap(2) 0.02 ap(4)],'ytick',clim{k});
  ch{k}(:,1)=line(zeros(2,ncnts(k)),zeros(2,ncnts(k)),'color',[1 1 1]);
  ch{k}(:,2)=line(zeros(2,ncnts(k)),zeros(2,ncnts(k)),'color',[1 1 1]);
  switch(typs{k})
    case 'u'
      cmap([fmin(k) -p.V*[1.01 0.99] 0.01*p.V*[-1 1] fmax(k)]);
    case 'b'
      cmap([0 0.01 0.99 1]*p.B);
    case 's'
      cmap([0 1]*p.S);
    case 'psi'
      cmap([fmin(k) -0.01 0.01 fmax(k)]);
  end
end
axes(ah(1));
ht=text(p.L/200,HT,'aaaa','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1],'BackgroundColor',[0 0 0 0.5]);
set(ah(end),'xtick',0:p.L);

%axes(a1)
%line(repmat([xp'-bs/2 xp'+bs/2],2,1),repmat([0;HT],1,2*length(xp)),'color',[0.7 0 0]);
%hr=line(repmat(xp,1,p.Nz)+NaN,c.Jz,'color',[0 0 1],'linewidth',2);
%axes(a2);
%line(repmat([xp'-vt xp' xp'+vt],2,1),repmat([0;HT],1,3*length(xp)),'color',[0.7 0 0]);
%hp=line(repmat(xp,1,p.Nz)+NaN,c.Jz,'color',[0 0 1],'linewidth',2);
if 0
  axes(a3);cla;
  set(a3,'xlim',[0 p.L],'ylim',[-1 1],'box','on');
  hh1=line(c.Jx,NULLx,'color',[1 0 0],'LineStyle','-');
  hh2=line(c.Jx,NULLx,'color',[0 0 1],'LineStyle','-');
  hh3=line(c.Jx,NULLx,'color',[0 1 0]/2,'LineStyle','-');
  hh4=line(c.Jx,NULLx,'color',[0 1 1]/3,'LineStyle','-');
  hh5=line(c.Jx,NULLx,'color',[1 1 0]/3,'LineStyle','-');
  line([0 p.L],[0 0],'color',0.7*[1 1 1]);
  line([0;p.L],[[1;1] [-1;-1] [1;1]/2],'color',0.7*[1 1 1],'LineStyle','--');
  legend([hh1 hh2 hh3 hh4 hh5],{'b','z','bb','bu','ev'},'location','eastoutside');
  if isempty(g)
    mina=min([d.b(:);d.p(:)]);
    maxa=max([d.b(:);d.p(:)]);
  else
    mina=min([d.b(:);d.p(:);d.bu(:);d.ev(:);g.bz(:)]);
    maxa=max([d.b(:);d.p(:);d.bu(:);d.ev(:);g.bz(:)]);
  end
  rga=(mina+maxa)/2+1.05*(maxa-mina)/2*[-1 1];
end


%tt=unique([a.t;d.t]);
tt=unique([a.t]);

x=c.Jx;
fao=0;
fbo=0;
ts=ded_interp_stats(ss,p,tt);
nfns=length(fns);
for j=1:nfns
  a=ded_read_hdf(fns{j});
  if isfield(a,'t1');  a.t=(a.t1+a.t2)/2; end
  if ddt;  dt=a.dt; end
  scl=scly/dt;
  set(ht,'string',a.t);
  disp(sprintf('%3i/%3i %6.2f',j,nfns,a.t));
  for k=1:nfld
    switch(typs{k})
      case('psi')
        if ndims(a.u)==3 ; a.u=squeeze(mean(a,u,2));a.w=squeeze(mean(a,w,2));end
        F=ded_psi(scl*a.u,scl*a.w,uparity,wparity,p.H,c.dJx);
      otherwise
        F=scl*a.(typs{k});
    end
    if ndims(F)==3 ; F=squeeze(mean(F,2));end
    F=ichebc2f(ichebf2c(F,dim),dim,z,p.H);
    if isfield(pp.clip,typs{k})
      minF=pp.clip.(typs{k})(1);
      maxF=pp.clip.(typs{k})(2);
      F=min(maxF,max(minF,F));
    else
      minF=min(F(:));
      maxF=max(F(:));
      if minF==maxF;maxf=minF+eps*max(eps,abs(minF));end;
    end
    axes(ah(k));
    set(ih(k),'CData',F);
    set(ah(k),'clim',[minF maxF]);
    set(icb(k),'ytick',unique([minF 0 maxF]));
    if ischar(cnts{k})
      switch(cnts{k})
        case 'range'
          cc=linspace(min(F(:)),max(F(:)),ncnts(k))
        case 'righteven'
          cc=interp(z,F(:,end,j),linspace(z(1),z(end),ncnts(k)),'linear');
        case 'rightlin'
          cc=linspace(F(1,end),F(end,end),ncnts(k));
        case 'linear'
          cc=linspace(cc(1),cc(2),ncnts(k));
        case 'percentile'
          cc=percentile(F(:),linspace(1,99,ncnts(k)));
        case 'otherwise'
          cc=percentile(F(:),linspace(1,99,ncnts(k)));
      end
    else
      cc=cnts{k};
    end
    for i=1:min(length(cc),ncnts(k))
      nnc=size(ch{k},2);
      [cx,cy]=longest_contours(x,z,F,cc(i),nnc);
      if nnc==1
        set(ch{k}(i),'Xdata',cx,'Ydata',cy);
      else
        for ii=length(cx)+1:nnc;set(ch{k}(i,ii),'Xdata',[],'Ydata',[]);end
        for ii=1:length(cx);set(ch{k}(i,ii),'Xdata',cx{ii},'Ydata',cy{ii});end
      end
    end
  end
  drawnow;
  
  %  pause(0.1);
end

return

function b=cmap(cc)


sc=cc-cc(1);
sc=sc/max(sc);
dc=diff(sc);
ctol=100*eps;
f=linspace(cc(1),cc(end),max(512,round(10/min(dc))));
c=[];
x=[];
n=64;
k=0;
for j=1:length(cc)-1
  if dc<=0.01
    m=2;
    c=zeros(m,3);
  else
    k=k+1;
    m=n;
    switch(k)
      case(1); c=[c;parula(n)]; 
      case(2); c=[c;hot(n)]; 
      case(3); c=[c;hsv(n)]; 
      case(4); c=[c;bone(n)];
      case(5); c=[c;cool(n)];
    end
  end
  x=[x;linspace(cc(j),cc(j+1)-ctol,m)'];
end
x(end)=cc(end);
b=interp1(x,c,f,'linear');
colormap(gca,b);
set(gca,'clim',cc([1 end]));
