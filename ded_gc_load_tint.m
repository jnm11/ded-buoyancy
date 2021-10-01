function a=ded_gc_load_tint(num,T,xrg)

if nargin<2
  T=[];
end
if nargin<3
  xrg=[];
end

if isempty(T)
  p=ded_gc_stats(num);
  T=p.T;
end

if isempty(xrg)
  xrg=[-inf inf];
end
  
typ={'ab','ap','au','aw','auu','avv','aww','auw'};
a=ded_read_2d(['gc/' num],'tint',typ);

if length(a.t)<2
  disp(sprintf('ded_gc_tint: tint data must have at least two entries %s',num));
end

if ~isfinite(T)
  T=a.t(end-1);
end

fx=find(a.x>=xrg(1) & a.x<=xrg(2));

T=min(a.t(end),T);
m=min(find(a.t>T))-1;
dt=a.t(end)-a.t(m);
a.x=a.x(fx);

if isempty(dt)
  disp(sprintf('ded_gc_tint: No converged times %7.3f %7.3f',T,a.t(end)));
  return;
end

ntyp=length(typ);
for k=1:ntyp
  n=typ{k};
  a.(n)=(a.(n)(:,fx,end)-a.(n)(:,fx,m))/dt;
end
