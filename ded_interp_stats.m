function [ts r]=ded_interp_stats(s,p,t,nn)

if nargin<3
  t=[];
end
if nargin<4
  nn=[];
end
if isempty(nn)
  nn=0;
end
if ~isfield('U',p)
  p.U=0;
end
if ~isfield('g',p)
  p.g=0;
end
if p.Ny==1
  p.W=0;
  p.Ny=0;
end
p.dx=0;
p.dy=0;
p.dz=0;
dtol=1e4;
if p.Nx>0
  p.dx=round(dtol*p.L/p.Nx)/dtol;
end
if p.Ny>0
  p.dy=round(dtol*p.W/p.Ny)/dtol;
end
if p.Nz>0
  p.dz=round(dtol*p.H/p.Nz)/dtol;
end
if strcmp(p.Tz,'Cheb')
  p.dz=p.dz*3/2;
end
ss={'g','U','Nx','Ny','Nz','L','H','W','dx','dy','dz'};
s1=p.name;
for j=1:length(ss)
  if p.(ss{j}) ~=0; s1=sprintf('%s, %s=%g',s1,ss{j},p.(ss{j}));end
end
  
ts={};

if ~isempty(s) 
  if isempty(t)
    t=s.t;
  end
  if isempty(s.t)
    s=[];
  end
end

if isempty(s) 
  for j=1:length(t)
    ts{j}=addnl(sprintf('%s, t=%6.2f',s1,t(j)),nn);
  end
  return
end
if isfield(s,'X')
  mt=(s.t(2:end)+s.t(1:end-1))/2;
  if isempty(mt)
    V=0;
  else
    V=diff(s.X)./diff(s.t);
    if length(mt)==1
      s.V=repmat(V,size(s.t));
    elseif length(mt)>0
      s.V=interp1(mt,V,min(max(mt),max(min(mt),s.t)),'pchip');
    else
      s.V=[];
    end
    V=1;
  end
end
t=double(t);
t=max(min(s.t),min(max(s.t),t));
nms=setdiff(fieldnames(s),{'t','name'});
for j=1:length(nms)
  r.(nms{j})= interpx(s,nms{j},t);
end

if ~isfield(r,'g');  r.g=repmat(p.g,size(t));end;
if ~isfield(r,'U');  r.U=repmat(p.U,size(t));end;
if ~isfield(r,'X');  r.X=repmat(NaN,size(t));end;
if ~isfield(r,'NS'); r.NS=repmat(NaN,size(t));end;
if ~isfield(r,'Re'); r.Re=repmat(NaN,size(t));end;

G=any(diff(r.g)~=0);
U=any(diff(r.U)~=0);
X=any(diff(r.X)~=0);
NS=any(diff(r.NS)~=0);
Re=any(diff(r.Re)~=0);
UE=isfield(r,'UE');
FU=isfield(r,'FU');

if UE
    UE=any(diff(s.UE)~=0);
end
if FU
  FU=any(diff(s.FU)~=0);
end

if strcmp(p.sType,'pm')
  G=0;
  X=0;
  U=0;
end

r.t=t;
for j=1:length(t)
  s2=sprintf('t=%6.2f, dt=%7.5f',t(j),r.dt(j));
  if p.U>0
    s2=sprintf('%s, Re=%5.0f',s2,p.U*p.H*p.Re); 
  else
    if isfield(r,'V')
      s2=sprintf('%s, V=%6.4f',s2,r.V(j)); 
      s2=sprintf('%s, Re=%5.0f',s2,r.V(j)*p.H*p.Re); 
    end
  end
  if UE; s2=sprintf('%s, E=%7.5f',s2,r.UE(j));  end  
  if FU; s2=sprintf('%s, E=%7.5f',s2,r.FU(j));  end  
  if NS; s2=sprintf('%s, N=%5.3f',s2,r.NS(j));  end
  if G;  s2=sprintf('%s, g=%6.4f',s2,r.g(j)); end
  if U;  s2=sprintf('%s, U=%6.4f',s2,r.U(j)); end
  if X;  s2=sprintf('%s, X=%6.4f',s2,r.X(j)); end
  
  st=sprintf('%s, %s',s1,s2);
  st=addnl(st,nn);
  ts{j}=st;
end


function st=addnl(st,nn)
if nn>0
  f=find(st==',');
  g=linspace(1,length(st),nn+2);g([1 end])=[];
  for k=length(g):-1:1
    h=findmin(abs(g(k)-f));
    st(f(h)+1)=[];
    st(f(h))=10;
    f(h)=[];
  end
end


function x=interpx(s,nm,t)
x=repmat(NaN,size(t));
if isfield(s,nm)
  n=length(s.(nm));
  if n==1
    x=repmat(s.(nm),size(t));
  elseif length(s.t)==n
    x  = interp1(s.t,double(s.(nm)),t,'linear');
  end
end
