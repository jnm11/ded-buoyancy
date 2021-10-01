function r=ded_gc_profiles_slice(nm,X,xrg,typ,display)
%r=ded_gc_profiles_slice(nm,X,xrg,display) Calculate profiles at a fixed position as a they vary in time
% Used for lock exchange simulation
% Works from b and u files
% Also calculates front position

if nargin< 5
  display=[];
end
if nargin< 4
  typ=[];
end
if isempty(display)
  display=false;
end
if isempty(typ)
  typ='bu';
end

if iscell(nm)
  for j=1:length(nm)
    ded_gc_profiles_slice(nm{j},X,xrg,typ,display);
  end
  return;
end

fnmat=[ded_dedalus_data_dir '/results/' nm '/slice.mat'];
p=ded_read_param(nm);
dt=1;
avrg=false;
switch(typ)
  case 'bu'
    fnu=ded_get_fn(nm,'u');
    fnb=ded_get_fn(nm,'b');
  case 'ay'
    fnu=ded_get_fn(nm,'ay');
    fnb=fnu;
    avrg=true;
  case 'y'
    fnu=ded_get_fn(nm,'y');
    fnb=fnu;
   otherwise
    error(sprintf('ded_gc_profiles_slice: unknown type %s',typ));
end

if file_nt(fnmat,fnu{end})
  disp(sprintf('ded_gc_profiles_slice: %s up to date',nm));
  if nargout>0
    r=load(fnmat);
    r=r.a;
  end
  return;
end
if nargin==1 & isfile(fnmat)
  r=load(fnmat);
  r=r.a;
  return;
end

disp(sprintf('ded_gc_profiles_slice: %s making',nm));

p=ded_read_param(nm);
c=ded_coord(nm);
s=ded_read_stats(nm);
trg=[s.t(min(find(s.X>=xrg(1)))) s.t(min(find([s.X;inf]>xrg(2)))-1)];
if isempty(p);
  disp(sprintf('ded_gc_profiles_slice: %s no data'));
  return;
end
W=p.W;


x=c.Jx;
z=c.Jz;
nu=length(fnu);
nb=length(fnb);
n=min(nu,nb);
nz=length(z);
nx=length(x);
u=zeros(nz,n);
t=zeros(1,n);

b  = zeros(nz,n);
bI = zeros(nx,n);
bt = zeros(1,n);
h  = zeros(1,n);

dz=min(diff(z));
dx=mean(diff(x));

ws=ichebintw(nz)*p.H/2;
for j=1:n
  disp(fnu{j});
  a=ded_read_hdf(fnu{j});
  if c.dim==3; a.u=mean(a.u,2); end
  if avrg; dt=a(j).dt;end
  u(:,j)=interp1(x,a.u',X,'linear')/(dt*W);
  if isfield(a,'t')
    t(j)=a.t;
  else
    t(j)=(a.t1+a.t2)/2;
  end
end

for j=1:n
  disp(fnb{j});
  a=ded_read_hdf(fnb{j});
  if avrg; dt=a(j).dt;end
  b(:,j)=interp1(x,a.b',X,'linear')/(dt*W);
  bI(:,j)=ws*a.b;
end
f=find(t>=trg(1)&t<=trg(2));
t=t(f);
b=b(:,f);
bI=bI(:,f);
u=u(:,f);
n=length(t);


Ix=dx*(x>X);
Ix(x==X)=dx/2;Ix=Ix(:)';
M=Ix*bI;  % Total buoyancy in right hand side


b0=ws*b/p.H;
u0=ws*u/p.H;

ws=ws(:);
clf;
opt=optimset('display','none');
fb = @(P,z) (1-erf((z-P(1))/P(2)))/2*p.B;
for j=1:n
  ff = @(p) (fb(p,z)-b(:,j)).*ws;
  [P r.br(j)]=lsqnonlin(ff,[p.H/2,p.H/2],[0 dz],[p.H 10*p.H],opt);
  r.bh(j)=P(1);
  r.bw(j)=P(2);
  r.b1(j)=fb(P,0);
  r.b2(j)=fb(P,p.H);
  if display
    plot(z,fb(P,z),z,b(:,j));
    axis([0 p.H 0 p.B]);
    title(nm);
    xlabel('z');
    ylabel('b');
    drawnow;
  end
end

fu = @(P,z) -P(3)*erf((z-P(1))/P(2));
minu=min(u(:));
maxu=max(u(:));
if minu==maxu;maxu=minu+eps*(1+abs(minu));end;
for j=1:n
  ff = @(P) (fu(P,z)-u(:,j)).*ws;
  [P r.ur(j)]=lsqnonlin(ff,[p.H/2,p.H/2,u(1,j)],[0 dz -inf ],[p.H 10*p.H inf],opt);
  r.uh(j)=P(1);
  r.uw(j)=P(2);
  r.u1(j)=fu(P,0);
  r.u2(j)=fu(P,p.H);
  
  if display
    plot(z,fu(P,z),z,u(:,j));
    axis([0 p.H minu maxu]);
    title(nm);
    xlabel('z');
    ylabel('u');
    drawnow;
  end
end

P=[1e-2 p.H/2];
display=true;
xmin=p.L/2;
for j=1:n
% $$$   aa=ded_gc_find_head(x,bI(:,j),xmin,display)
% $$$   r.Xf(j) =       aa.Xf;  
% $$$   r.hf(j) =       aa.hf;  
% $$$   r.Xh(j) =       aa.Xh;  
% $$$   r.hh(j) =       aa.hh;  
% $$$   r.Xt(j) =       aa.Xt;  
% $$$   r.ht(j) =       aa.ht;  
% $$$   r.alpha(j) =    aa.alpha;  
% $$$   r.X0(j) =       aa.X0;  
% $$$   r.h0(j) =       aa.h0;  
% $$$   r.alpha1(j) =   aa.alpha1;  
% $$$   r.alpha2(j) =   aa.alpha2;    

  ff = @(P) hx(P(1),P(2),x,p.H,p.L)-bI(:,j);
  [P r.hr(j)]=lsqnonlin(ff,P,[0 p.H/4],[p.L/2 p.H],opt);
  if display
    plot(x,bI(:,j),x,hx(P(1),P(2),x,p.H,p.L));
    axis([0 p.L 0 p.H]);
    title(nm);
    xlabel('z');
    ylabel('bI');
    drawnow;
  end
  r.Xf(j)=P(1);
  r.ht(j)=P(2);
end
if display
  subplot(5,1,1);
  h=plot(t,r.bw);
  axis([0 inf 0 3*median(r.bw)]);
  ylabel('bw');
  xlabel('t');
  
  subplot(5,1,2);
  h=plot(t,r.uw);
  axis([0 inf 0 2*median(r.uw)]);
  ylabel('uw');
  xlabel('t');
  
  subplot(5,1,3);
  h=plot(t,r.u1-r.u2);
  axis([0 inf 0 inf]);
  ylabel('u1-u2');
  xlabel('t');
  
  subplot(5,1,4);
  h=plot(t,r.Xf);
  axis([0 inf 0 inf]);
  ylabel('X');
  xlabel('t');
  
  subplot(5,1,5);
  h=plot(t,r.ht);
  axis([0 inf 0 inf]);
  ylabel('S');
  xlabel('t');
end

r.t=t;
r.bu=ws'*(b.*u);
r.M=M;
r.bI=bI;
r.b=b;
r.u=u;
r.x=x;
r.z=z;

if n>2
  p=polyfit(r.t,r.Xf,1);
  r.V=p(1);
else
  r.V=NaN;
end
r.Xf  = @(t) p(1)*t+p(2);
r.Uf  = @(t) p(1);
r.trg = [t(1) t(end)];

mkdirifnotexist(fnmat);
a=r;
save(fnmat,'a');

return;

function r=hx(X,h,x,H,L)
z=x-L/2;
s=sign(z);
a=(h-H/2)/X;
Z=max(X-abs(z),0)/h;
b=exp(-sqrt(3)*Z);
r=H/2*(1-s.*b+(1-b).*min(H/2,max(-H/2,a*z)));
return;

cd('~/');
fns=cellstr_ls('gc/le/a/*/*/*/status',[],'dir');
a=ded_gc_profiles_slice(nm,24,[2 24]);

