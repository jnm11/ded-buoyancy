function  [a fnmat]=ded_gc_fit_erf_profiles(nm,qqq,typ,force,display)

%a=ded_gc_fit_erf_profiles('gc/f7/mf/00200/25801',[],'final',[],1);
%for a in $(seq 1 1000) ; do rsync -vap --progress asahi:gc/f7/ma/00100/27100 ~/gc/f7/ma/00100 --exclude a --exclude ayz --exclude avar --exclude sev  --exclude final --exclude force --bwlimit 250 ; done
%for a in $(seq 1 1000) ; do rsync -vap --progress asahi:gc/f7/mb/01000/20986 ~/gc/f7/mb/01000 --exclude a --exclude ayz --exclude avar --exclude sev  --exclude final --exclude force  ; done
aa=[];

%nm='gc/f7/ma/00100/27100';
%nm='gc/f7/mb/01000/20986';
%qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f7/mb/01000/20986',q)
%qqq.trg=[10 inf];a=ded_gc_fit_erf_profiles('gc/f7/mf/00200/25801',[],1)
if nargin==0
  cd('~/');
  fns=cellstr_ls('gc/f6/g/*/*/*/[23]*/status',[],'dir');
  ded_gc_fit_erf_profiles(fns,[],'ay');
  return;
  qqq.trg=[  10 inf];a=ded_gc_fit_erf_profiles('gc/f6/f/21',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/22',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/23',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/24',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/25',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/30',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/31',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/32',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/33',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/34',q);
  qqq.trg=[-inf inf];a=ded_gc_fit_erf_profiles('gc/f6/f/35',q);
  return;
end


if nargin<2; q=struct; end;
if nargin<3; typ=[];   end;
if nargin<4; force=[]; end;
if nargin<5; display=[]; end;

if isempty(display)
  display=false;
end

if iscell(nm)
  a=struct;
  for j=1:length(nm)
    aa=ded_gc_fit_erf_profiles(nm{j},qqq,typ,force,display);
    if ~isempty(aa)
      aa=rmfieldifexist(aa,{'bd1',  'bd2',  'bdI',  'd1' ,  'd2' ,  'dI' ,  'ud1',  'ud2',  'udI',  'wd1',  'wd2',  'wdI'});
      a=struct_array_append(a,aa,nm{j},true);
    end
  end
  return
end
if ischar(nm)
  if any(nm=='*');
    cd(ded_dedalus_data_dir);
    fns=cellstr_ls([nm '/param.h5'],[],'dir');
    a=ded_gc_fit_erf_profiles(fns,qqq,typ,force,display);
    return;
  end
end

if isstruct(nm)
  a=nm;
  nm=a.nm;
else
  a=[];
end

if isempty(force); force=false; end

if isstruct(a)
  fnmat = [ded_dedalus_data_dir '/results/' a.fnmat '.mat'];
  if file_nt(fnmat,a.fn) 
    disp(sprintf('ded_gc_fit_erf_profiles: up to date "%s"',fnmat));
    load(fnmat);
    return;
  end
else
  fnmat = [ded_dedalus_data_dir '/results/' nm '/profile-' typ '.mat'];
  if typ(1)=='r'
    fna   = sprintf('%s/results/%s/%s.mat',ded_dedalus_data_dir,nm,typ(2:end));
  else
    switch(typ)
      case('final')
        fna  = ded_get_fn(nm,typ,[],'state');
      otherwise
        fna  = ded_get_fn(nm,typ);
    end
    if isempty(fna);
      a=[];
      return;
    end
    fna=fna{end};
  end
  if isempty(fna);
    a=[];
    return;
  end
  
  if ~isfile(fna) & isempty(a)
    disp(sprintf('ded_gc_fit_erf_profiles: no file %s',fna));
    return;
  end
  if file_nt(fnmat,fna) & isempty(a) & ~force 
    disp(sprintf('ded_gc_fit_erf_profiles: up to date "%s"',fnmat));
    load(fnmat);
    return;
  end
end
if isempty(a)
  [dd fn ext]=fileparts(fna);
  switch(ext)
    case('.mat')
      load(fna);
    case('.hdf5')
      a=ded_read_hdf(fna);
  end
end

nn=nm;nn(nn=='/')='-';
nn=['~/pics/mdns/' nn];

if ~isfield(a,'b')
  disp(sprintf('ded_gc_fit_erf_profiles: %s no data %s',nm,typ));
  return;
end
disp(sprintf('ded_gc_fit_erf_profiles: making "%s"',fnmat));

f=fieldnames(a);
for j=1:length(f)
  if ndims(a.(f{j}))==3
    a.(f{j})=squeeze(mean(a.(f{j}),2));
  end
end

p=ded_read_param(nm);

if isfield(a,'x')
  x=a.x;
  z=a.z;
  Nx=length(x);
  Nz=length(z);
else
  c=ded_coord(nm);
  z=c.Jz(:);
  x=c.Jx(:);
  Nx=c.NJx;
  Nz=c.NJz;
end

a=rmfieldifexist(a,{'d','bd','sd','ud','vd','wd'});

if typ(1)~='a';
  if ~isfield(a,'bb') a.bb=a.b.*a.b;end
  if ~isfield(a,'bu') a.bu=a.b.*a.u;end
  if ~isfield(a,'bw') a.bw=a.b.*a.w;end
  if ~isfield(a,'ww') a.ww=a.w.*a.w;end
  if ~isfield(a,'uu') a.uu=a.u.*a.u;end
  if ~isfield(a,'uw') a.uw=a.u.*a.w;end
end
  
H=p.H;
V=p.V;

nz=size(a.b,1);
w=ichebintw(nz)*H/2;
[w1 w2]=ichebendsw(nz);

ws=sqrt(w(:));
b0=w*a.b;
bz=(z(:)'.*w)*a.b;

if isfield(a,'dbdx')
  dbs=w*(a.dbdx.^2+a.dbdz.^2);
  dus=w*(a.dudx.^2+a.dwdz.^2 + (a.dudz+a.dwdx).^2/2);
else
  a.c='zx';
  parity=ded_read_parity(nm);
  ad=ded_diff(a,p,parity,c,{'b','u','w'});
  dbs=w*(ad.dbdx.^2+ad.dbdz.^2);
  dus=w*(ad.dudx.^2+ad.dwdz.^2 + (ad.dudz+ad.dwdx).^2/2);
end

uu2=w*a.u.^2;
b1=1./(w*(a.b.*(1-a.b)))/sqrt(2*pi);
b1(b0<1e-2)=0;




opt=optimset('display','none');
q=[1 10];

figs=isfield(qqq,'sz');
if figs;
  clf;preprint(qqq.sz,qqq.fnt);colour_lines
end

f1=max(find(b0>0.05*max(b0)));
x1=x(f1);x2=x(f1+1);
ba=b0(f1);bb=b0(f1+1);
X=((x2-x1)*0.05+x1*bb-x2*ba)/(bb-ba);

r=struct;
f = @(p,z) p(3)*serf(z,p(1),p(2))/serf(0,p(1),p(2));

for j=1:Nx
  
  if max(a.b(:,j))<1e-2
    P=[0 inf 0];
  else;
    ff = @(p) (f(p,z)-a.b(:,j)).*ws;
    [P p.rb(j)]=lsqnonlin(ff,[1.1*b0(j),b0(j)/10,1],[b0(j) 0 0],[H H/2 1],opt);
  end
  r.bA(j)=0;
  r.bB(j)=P(3)/serf(0,P(1),P(2));
  r.bh(j)=P(1);
  r.bw(j)=P(2);
  r.b1(j)=f(P,0);
  r.b2(j)=f(P,H);
  if mod(j,10)==0 & (figs | display)
    plot(a.b(:,j),z,f(P,z),z,'--');
    axis([-0.10 1.1 0 H]);
    text(1,H*0.95,sprintf('x=%5.2f, h=%5.3f w=%5.3f b1=%5.3f ',X-x(j),r.bh(j),r.bw(j),r.b1(j)),'horizontalalignment','right','verticalalignment','top');
    set(gca,'ygrid','on','xgrid','on','box','on');
    drawnow;
    %    fnm=sprintf('%s-b-%05u.eps',nn,j);
    %print('-depsc2',fnm);
  end
  
end


minu=min(a.u(:));
maxu=max(a.u(:));
rgu=[minu maxu]+(maxu-minu)/10*[-1 1];

for j=1:Nx
  u=a.u(:,j);
  u0=w*u/H;
  u1=w1*u;
  u2=w2*u;
  if all(u==mean(u))
    r.uh(j)=0;
    r.uw(j)=0;
    r.u1(j)=mean(u);
    r.u2(j)=mean(u);
    continue
  end
  
  switch(p.lbc)
    case 'slip'
      %f= @(p,z) w*u/H+p(3)*(serf(z,p(1),p(2))-1/H);
      f= @(p,z) w*u/H+p(3)*serf1(z,p(1),p(1)*p(2),H);
      ff = @(p) (f(p,z)-u).*ws;
      lh=max(0,min(H/3,(u0-u2)/(u1-u2)));
      [P r.ur(j)]=lsqnonlin(ff,[1 1 1],[lh 0 0],[H 1 inf],opt);
      r.u1(j)=f(P,0);
      P(2)=P(1)*P(2);
    case 'noslip'
      urg=max(u)-min(u);
      %f= @(p,z) w*u/H+p(3)*(serf(z,p(1),p(2))-1/H)-p(5)*(exp(-p(4)*z)-1/(p(4)*H));
      f= @(p,z) w*u/H+p(3)*serf1(z,p(1),p(2),H)-p(5)*(exp(-p(4)*z)-1/(p(4)*H));
      ff = @(p) (f(p,z)-u).*ws;
      [P r.ru(j)]=lsqnonlin(ff,[H/2 H/10 0 100/H 0],[0 0 0 10 0],[H H urg inf urg],opt);
      r.u1(j)=w*u/H+P(3)*(serf(0,P(1),P(2))-1/H)-P(5)/(P(4)*H);
  end
  r.uA(j)=w*u/H;
  r.uB(j)=P(3);
  r.uh(j)=P(1);
  r.uw(j)=P(2);
  r.u2(j)=f(P,H);
  if mod(j,10)==0 & (figs | display)
    plot(u,z,f(P,z),z,'--',rgu,r.uh([j j]),rgu,r.uh([j j])+r.uw(j),rgu,r.uh([j j])-r.uw(j));
    axis([rgu 0 H]);
    text(rgu(2)*0.95+rgu(1)*0.05,H*0.95,sprintf('x=%5.2f, h=%5.3f w=%5.3f u1=%5.3f ',X-x(j),r.uh(j),r.uw(j),r.u1(j)),'horizontalalignment','right','verticalalignment','top');
    set(gca,'ygrid','on','xgrid','on','box','on');
    drawnow;
    %fnm=sprintf('%s-u-%05u.eps',nn,j);
    %print('-depsc2',fnm);
  end
end

r.Tuu=w*(a.u.*a.u);
r.Tuw=w*(a.u.*a.w);
r.Tww=w*(a.w.*a.w);
r.Tbu=w*(a.b.*a.u);
r.Tbw=w*(a.b.*a.w);
r.Tbb=w*(a.b.*a.b);

if isfield(a,'uu')
  r.Fuu=w*(a.uu-a.u.*a.u);
  r.Fuw=w*(a.uw-a.u.*a.w);
  r.Fww=w*(a.ww-a.w.*a.w);
  r.Fbu=w*(a.bu-a.b.*a.u);
  r.Fbw=w*(a.bw-a.b.*a.w);
  r.Fbb=w*(a.bb-a.b.*a.b);
end

f=fieldnames(a);
z=z(:)';
for j=1:length(f)   % 'c' - collocation   'm' - moment
  if size(a.(f{j}),1)==nz
    r.([f{j} 'c1']) =  w1       * a.(f{j});
    r.([f{j} 'c2']) =  w2       * a.(f{j});
    r.([f{j} 'm0']) = (w.*z.^0) * a.(f{j});
    r.([f{j} 'm1']) = (w.*z.^1) * a.(f{j});
    r.([f{j} 'm2']) = (w.*z.^2) * a.(f{j});
    r.([f{j} 'm3']) = (w.*z.^3) * a.(f{j});
    r.([f{j} 'm4']) = (w.*z.^4) * a.(f{j});
  end
end
r.x=x(:)';
r.nm=nm;
r.X=X;

r.bz=bz;
r.dbs=dbs;
r.dus=dus;
r.uu2=uu2;

a=r;

[dd fn]=fileparts(fnmat);
if ~isdir(dd); mkdir(dd); end

hn=(a.uh-a.bh)./sqrt(a.uw.^2+a.bw.^2);
hp=(a.uh+a.bh)./sqrt(a.uw.^2+a.bw.^2);

a.q  = a.bB.*a.uA-2*a.bB.*a.uB.*(sqrt(pi).*(erf(hp).*hp-hn.*erf(hn))+exp(-hp.^2)-exp(-hn.^2))./(sqrt(a.bw.^2+a.uw.^2).*sqrt(pi).*(hn.^2-hp.^2));
a.q1 = a.bB.*a.uA-2*a.bB.*a.uB.*(1/sqrt(pi)./hp-1)./(hp.*sqrt(a.bw.^2+a.uw.^2));
save(fnmat,'a');

%plot(a.x,a.q,a.x,a.buI,a.x,a.q1);axis([0 p.L -1e-2 1e-2]);

return;

plot(x,r.uh,x,r.bh);axis([0 p.L 0 p.H]);



fb=find(rb<2e-3 & b0>1e-2);
fb=min(fb):max(fb);
fb=rg;

figure(1);clf;
subplot(3,1,1);
plot(c.Jx(fb),bh(fb));
axis([-inf inf 0 max(bh).*1.05]);
ylabel('bh');

subplot(3,1,2);
plot(c.Jx(fb),bw(fb));
axis([-inf inf 0 inf]);
ylabel('bw');

subplot(3,1,3);
plot(c.Jx(fb),rb(fb));
axis([-inf inf 0 inf]);
ylabel('err');

m0=w.*a.u./p.H;
m1=((c.Jz.^1)'.*w).*(a.u-m0)./p.H;
m2=((c.Jz.^2)'.*w).*(a.u-m0)./p.H;
m3=((c.Jz.^3)'.*w).*(a.u-m0)./p.H;
dU=b.u(end,:)-b.u(1,:);
uh=3./2.*m2./m1-p.H;
uw=2.*m1.^2./max(1e-2,8.*m3.*m1-9.*m2.^2+6.*m2.*m1.*p.H-4.*p.H.^2.*m1.^2);
ru=0.*uh;
uh=zeros(c.Nx,1);
uw=zeros(c.Nx,1);
clf;preprint([1 2],8);colour_lines

for j=1:c.Nx
  f = @(q) ded_uerf(b.z(:),q(1),q(2),q(3),m0(j),p.H); % ded_uerf(y,h,w,dU,U,H,n)
  ff = @(q) f(q)-b.u(:,j); % ded_uerf(y,h,w,dU,U,H,n)
  q=[uh(j) uw(j) b.u(end,j)-b.u(1,j)];
  q=lsqnonlin(ff,q,[-1 0 -2.*p.V.*p.H],[2 100 0],opt);
  ru(j)=sqrt(mean(ff(q).^2));
  plot(f(q),b.z,b.u(:,j),b.z);
  uh(j)=q(1);
  uw(j)=1./q(2);
  dU(j)=q(3);
  axis([-2.1 .25 0 2]);
  drawnow;
  nm=sprintf('~/pics/sdns/%s-u-%05u.eps',nn,j);
  print('-depsc2',nm);
end



fu=find((ru<p.V./10) & (-dU>p.V./5) & (uh>p.H./20));
fu=min(fu):max(fu);
fu=rg;

figure(2);clf;
subplot(4,1,1);
plot(c.Jx(fu),uh(fu));
axis([-inf inf 0 max(uh).*1.05]);
ylabel('uh');

subplot(4,1,2);
plot(c.Jx(fu),uw(fu));
axis([-inf inf 0 inf]);
ylabel('uw');

subplot(4,1,3);
plot(c.Jx(fu),-dU(fu));
axis([-inf inf 0 inf]);
ylabel('dU');

subplot(4,1,4);
plot(c.Jx(fu),ru(fu));
axis([-inf inf 0 inf]);
ylabel('err');

figure(3);clf;
subplot(3,1,1);
h=plot(c.Jx(rg),uh(rg),c.Jx(rg),bh(rg));
axis([-inf inf 0 1.05]);
legend(h,{'u','b'});
ylabel('h');

subplot(3,1,2);
h=plot(c.Jx(rg),uw(rg),c.Jx(rg),bw(rg));
axis([-inf inf 0 max(uw(rg)).*1.05]);
legend(h,{'u','b'});
ylabel('w');

subplot(3,1,3);
h=plot(c.Jx(rg),(bh(rg)-uh(rg))./uw(rg),c.Jx(rg),(bh(rg)-uh(rg))./bw(rg));
axis([-inf inf 0 inf]);
legend(h,{'u','b'});
ylabel('bh-uh');



subplot(2,1,1);
wf=find(c.Jx>1&c.Jx<8);
f=@(p)p(1).*sqrt(p(2)-c.Jx(wf))-uw(wf);
p=lsqnonlin(f,[max(uw(wf)),max(c.Jx(wf))]);
pp=polyfit(c.Jx(wf),uw(wf),2);
plot(c.Jx(wf),uw(wf),c.Jx(wf),p(1).*sqrt(max(0,p(2)-c.Jx(wf))),c.Jx(wf),polyval(pp,c.Jx(wf)));


subplot(2,1,2);
wf=find(c.Jx>1&c.Jx<8);
f=@(p)p(1).*sqrt(p(2)-c.Jx(wf))-bw(wf);
p=lsqnonlin(f,[max(bw(wf)),max(c.Jx(wf))]);
pp=polyfit(c.Jx(wf),bw(wf),1);
plot(c.Jx(wf),bw(wf),c.Jx(wf),0.*p(1).*sqrt(max(0,p(2)-c.Jx(wf))),c.Jx(wf),polyval(pp,c.Jx(wf)));

           

parity=ded_read_hdf([ded_dedalus_data_dir './' nm './' 'parity.hdf5']);
a=ded_read_javrg(nm,'a',trg,'combine');
nms=intersect(fieldnames(a),nms);
da=ded_diff(a,p,parity,c,nms);
Ia=ded_int(a,p,parity,c,nms);




fns=cellstr_ls('~/data/dedalus/results/gc/f6/g/*/*/*/[23]*/profile.mat');

for j=1:length(fns)
  load(fns{j})
  if isfield(a,'u0')
    disp([fns{j} ': u0']);
  end
end
rm -rf data/dedalus/results/gc/f6/g/1000/0400/2304/29400 data/dedalus/results/gc/f6/g/1000/1600/2304/23100 data/dedalus/results/gc/f6/g/2000/0400/2304/28600 data/dedalus/results/gc/f6/g/4000/2000/4608/21500 data/dedalus/results/gc/f6/g/8000/2000/4608/21300 






if false
P=zeros(Nx,3);
for k=1:3
  switch(k)
    case 1
      uw=a.uw-a.u.*a.w;
      mm='uw';
    case 2
      uw=a.uu-a.u.*a.u;
      mm='uu';
    case 3
      uw=a.ww-a.w.*a.w;
      mm='ww';
  end
  minu=min(uw(:));
  maxu=max(uw(:));
  rgu=[minu maxu]+(maxu-minu)/10*[-1 1];
  P=repmat(NaN,Nx,3);
  for j=1:Nx
    u=uw(:,j);
    p0=[r.uh(j) 1/r.uw(j) u(findmax(abs(u)))];
    if any(~isfinite(p0))
      continue;
    end
    ff = @(p) (f(p)-u).*ws;
    P(j,:)=lsqnonlin(ff,p0,[],[],opt);
    rb(j)=sqrt(mean(ff(P(j,:)).^2));
    if mod(j,10)==0 & figs
      plot(uw(:,j),z,f(P(j,:)),z,'--');
      axis([rgu 0 p.H]);
      text(rgu(2)*0.95+rgu(1)*0.05,p.H*0.95,sprintf('x=%5.2f',X-x(j)),'horizontalalignment','right','verticalalignment','top');
      set(gca,'ygrid','on','xgrid','on','box','on');
      drawnow;
      fnm=sprintf('%s-%s-%05u.eps',nn,mm,j);
      print('-depsc2',fnm);
    end
  end
  P(x>=X,:)=NaN;
  
  switch(k)
    case 1
      r.Fuw=w*uw;
      r.uwh=P(:,1)';
      r.uww=1./P(:,2)';
      r.ruw=rb(:)';
      r.Auw=P(:,3)';
    case 2
      r.Fuu=w*uw;
      r.uuh=P(:,1)';
      r.uuw=1./P(:,2)';
      r.ruu=rb(:)';
      r.Auu=P(:,3)';
    case 3
      r.Fww=w*uw;
      r.wwh=P(:,1)';
      r.www=1./P(:,2)';
      r.rww=rb(:)';
      r.Aww=P(:,3)';
  end
end
end



cd ~/data/dedalus/results/gc/f6/g/8000/2000/4608/21400


load('ay.mat');

nm='gc/f6/g/8000/2000/4608/21400';
ded_gc_fit_erf_profiles(nm,[],'ray',true);



% $$$ else
% $$$   z=a.z(:);
% $$$   x=a.x(:);
% $$$   Nx=length(x);
% $$$   Nz=length(z);
% $$$   c.dd=[mean(diff(x)) 0 mean(diff(z))];
% $$$   c.c= 'xyz';
% $$$   c.dim=3;
% $$$ end
% $$$ if isfield(a,'ub') & ~isfield(a,'bu')
% $$$   a.bu=a.ub;
% $$$   a=rmfield(a,'ub');
% $$$ end
% $$$ if isfield(a,'vb') & ~isfield(a,'bw')
% $$$   a.bw=a.vb;
% $$$   a=rmfield(a,'vb');
% $$$ end
% $$$ if isfield(a,'wb') & ~isfield(a,'bw')
% $$$   a.bw=a.wb;
% $$$   a=rmfield(a,'wb');
% $$$ end
% $$$ mfld={'bb'};
% $$$ for j=1:length(mfld)
% $$$   if ~isfield(a,mfld{j});
% $$$     a.(mfld{j})=NaN*a.b;
% $$$   end
% $$$ end


%parity=ded_parity_rmy(parity);
%x=c.Jx(:)';z=c.Jz(:);k=7*pi/p.L;o=randn(1);
%a.u=sin(k*x)+z.^2.*sin(2*k*x);ad=ded_diff(a,p,parity,c,{'b','u','w'});plot(x,ad.dudx-k*cos(k*x)-2*k*z.^2.*cos(2*k*x))
%a.w=cos(k*x)+z.^2.*cos(2*k*x);ad=ded_diff(a,p,parity,c,{'b','u','w'});plot(x,ad.dwdx+k*sin(k*x)+2*k*z.^2.*sin(2*k*x))
%plot(x,ad.dudz-2*z.*sin(2*k*x),x,ad.dwdz-2*z.*cos(2*k*x))



nm='gc/f7/mf/00200/25801';

nm='gc/f7/mg/4000/0250/4608/36591';
nm='gc/f7/mg/4000/0250/4608/36592'
b=ded_gc_fit_erf_profiles(nm,[],'final');
fnf=ded_get_fn(nm,'final',[],'state');
a=ded_read_hdf(fnf{end});
p=ded_read_param(nm);
c=ded_coord(nm);
w=ichebintw(size(a.w,1));
plot(c.x,w*a.w.^2);

dt=1;
%[a.psi a.omega a.div]=ded_psi(a.u/(p.W*dt),a.w/(p.W*dt),a.u_parity,a.w_parity,p.H,c.dJx);
%d=ded_zgrid(a,p.Nz,[],[],[],[],p.H,1);
%contour(c.Jx,d.z,d.psi,d.psi(5:10:end,10));



ded_gc_2db(nm,{'b','psi'},'ay');

 
mpiexec -n 30 ded_gc.py --pfn f7/mf/00200/25801 --rfn f7/mf/00200/25801 --reset --wu 0.4 gc/f7/mf/00200/25802
mpiexec -n 30 ded_gc.py --pfn f7/mf/00200/25801 --rfn f7/mf/00200/25801 --reset --wu 0.5 gc/f7/mf/00200/25803
mpiexec -n 32 ded_gc.py --pfn f7/mg/4000/0250/4608/36591 --rfn  f7/mg/4000/0250/4608/36591 --reset --wu 0.5  --U1 0.4 gc/f7/mg/4000/0250/4608/36592
mpiexec -n 32 ded_gc.py --pfn f7/mg/4000/0250/4608/36591 --rfn  f7/mg/4000/0250/4608/36591 --reset --ddiv True --dU1 True --dtjavar 0.01 --dtjsev 0.01 gc/f7/mg/4000/0250/4608/36593



nms={'gc/f7/mg/4000/0250/4608/36591','gc/f7/mg/4000/0250/4608/36592'};
for j=1:length(nms)
 fnf=ded_get_fn(nms{j},'final',[],'state');
 a=ded_read_hdf(fnf{end});
 ww(j,:)=w*a.w.^2;
end
plot(ww(:,1:20)')

 
 b=ded_gc_fit_erf_profiles(nms{j},[],'final');
fnf=ded_get_fn(nm,'final',[],'state');
a=ded_read_hdf(fnf{end});
p=ded_read_param(nm);
c=ded_coord(nm);
w=ichebintw(size(a.w,1));
plot(c.x,w*a.w.^2);


nm='gc/f7/mg/4000/0250/4608/36593';
avar=ded_read_g(nm,'avar');
sev=ded_read_g(nm,'sev');

