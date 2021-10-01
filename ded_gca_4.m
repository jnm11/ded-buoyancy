function [a p]=ded_gca_4(nm,trg)
% ded_gca_3 depth integrated analysis of gravity current
%nm='gc/f6/g/26';
%nm='gc/test/50';
if nargin<2
  trg=[];
end
if isempty(trg)
  trg=[-inf inf];
end
f=ded_read_g(nm,'force');
[a p]=ded_gc_read_javrg(nm,trg);
T=a.dt;
trg=[a.t1 a.t2];
PS=p.B*p.g*p.hb;
nu= 1/p.Re;
xdmax=max([p.x0 p.x1 p.x2 p.x3 p.x4 p.x5 p.x6  p.x7]);
xrg=[0 p.L];
prg=[0 p.L -0.5 0.5];
xdmax=4;
fx=find(a.x>xdmax);
x=a.x(fx);
xrg=[x(1) x(end)];
dimz=1;
dimx=2;
dx=a.x(2)-a.x(1);
nz=length(a.z);
nx=length(a.x);

w0=ichebIw(nz)*(p.H/2);
w1=w0.*a.z;
w2=w0.*a.z.^2;
sz=[nz nx];
a.p=a.p-w0'*a.p(:,end)/p.H;


fld=fieldnames(a);
for j=1:length(fld)
  b=fld{j};
  if all(size(a.(b))==sz)
    a.([b 'zz'])=w2'*a.(b);
    a.([b  'z'])=w1'*a.(b);
    a.( b      )=w0'*a.(b);
  end
end
a.b=max(a.b,0);
tol=1e-4;
bd   = max(a.b, tol);
a.bz  = max(a.bz,0);
a.bzz = max(a.bzz,0);
a.bb  = max(a.bb,1e-7);
bz  = a.bz./bd;
hzz = sqrt(12*max(0,a.bzz/bd-bz.^2));
hz=2*bz;
hbb=bd.^2./max(a.bb,tol);
f=find(a.b<tol);
hb=max(0,a.b);
bz(f)=0;
hz(f)=0;
hzz(f)=0;
hbb(f)=0;


figure(1);clf;
h=plot(a.x,a.b,a.x,hz,a.x,hbb);
legend(h,{'b','hz','bb'},'location','best');

uu=max(0,a.uu/(p.H*p.U^2)-1);
vv=max(0,a.vv/(p.H*p.U^2));
ww=max(0,a.ww/(p.H*p.U^2));

figure(2);clf;
h=plot(a.x,uu,a.x,vv,a.x,ww);
legend(h,{'uu','vv','ww'});

%figure(3);clf;
%plot(a.x,(a.uu+a.p)/(p.H*p.U^2)-1);

