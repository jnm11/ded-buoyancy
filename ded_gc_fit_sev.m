function a=ded_gc_fit_sev(nm,typs)

if nargin<2
  typs=[];
end
if isempty(typs)
  typs={'b1','b2','bh','bw','bt','be','u1','u2','uh','uw','ut'};
  typs={'uh','uw','bh','bw'};
end

if ~iscell(nm)
  if any(nm=='*')
    nm=cellstr_ls([nm '/param.h5'],[],'dir');
  end
end

if iscell(nm)
  for j=1:length(nm)
    figure(j);
    a(j)=ded_gc_fit_sev(nm{j},typs);
  end
  return;
end

s=ded_read_g(nm,'sev');

if isempty(s)
  a=[]
  return;
end

p=ded_read_param(nm);
c=ded_coord(nm);

nz = c.NAz;
nt = length(s.t);
w  = ichebintw(nz)*p.H/2;
db = isfield(s,'db');
du = isfield(s,'divz');

a.t=s.t;
if db
  a.b1=zeros(nt,1);
  a.b2=zeros(nt,1);
  a.bh=zeros(nt,1);
  a.bw=zeros(nt,1);
  for j=1:nt
    [bp bn br] = fit_serf(c.Az,s.db(:,j),w);
    a.b1(j) = bp(1);  
    a.b2(j) = bp(2);
    a.bh(j) = bp(3);  
    a.bw(j) = bp(4);
    a.br(j) = br;
  end
end

if du
  a.u1=zeros(nt,1);
  a.u2=zeros(nt,1);
  a.uh=zeros(nt,1);
  a.uw=zeros(nt,1);
  for j=1:nt
    [up un ur] = fit_serf(c.Az,s.divz(:,j),w);
    a.u1(j) = up(1);  
    a.u2(j) = up(2);
    a.uh(j) = up(3);  
    a.uw(j) = up(4);
    a.ur(j) = ur;
  end
end

typs=intersect(typs,fieldnames(a));
n=length(typs);
figure;clf;
for j=1:n
  m=typs{j};
  subplot(n,1,j);
  cla;
  plot(a.t,a.(m));
  ylabel(m);
  axis('tight');
end
subplot(n,1,1);
title(nm);
