function a=ded_fix_stats(nm,trg)
if nargin<2
  t0=[];
end
if isempty(trg)
  trg=[-inf inf];
end

p=ded_read_param(nm);
if isempty(p)
  disp(sprintf('ded_fix_stats: %s: no param file',nm));
  return;
end
if strcmp(p.status,'Running')
  disp(sprintf('ded_fix_stats: %s: Simulation is running',nm));
  return;
end  

[a fns]=ded_read_stats(nm,0); % 0 zero mean read without fixing
if isempty(a)
  disp(sprintf('ded_fix_stats: %s: no stats file',nm));
  return;
end


if isfield(p,'dtstat')
  dt=p.dt;
else
  dt=mode(round(diff(a.t)*1000)/1000);
end
n=round(a.t/dt);



minn=max(trg(1)/dt,min(n));
maxn=min(trg(2)/dt,max(n));

n(n<minn | n>maxn)=[];

if isempty(n)
  disp(sprintf('ded_fix_stats: %s: No remaining times',nm));
  unix(sprintf('/bin/rm -f %s',fns));
  return;
end

t=(min(n):max(n))*dt;

disp(sprintf('ded_fix_stats: %s: clipping time %7.3f %7.3f',nm,t(1),t(end)))

mm=length(t);
f=zeros(mm,1);
e=zeros(mm,1);
for j=1:mm
  [e(j) f(j)]=min(abs(t(j)-a.t));
end
g=find(e<0.2*dt);
nf=length(f)-length(g);
if nf==0 & length(ft)==0
  disp(sprintf('ded_fix_stats: %s: no fixes necessary',nm));
  return;
end
if nf>=0
  disp(sprintf('ded_fix_stats: %s: %u fixes',nm,nf));
end

f=f(g);
e=e(g);

fnm=fieldnames(a);
for j=1:length(fnm)
  c=fnm{j};
  x=a.(c)(f);
  n=sum(~isfinite(x));
  if n>0
    disp([sprintf('ded_fix_stats: %s: field %s: %u non finite t=',nm,c,n) sprintf('%7.3f ',a.t(f(~isfinite(x))))]);
  end
  b.(c)=a.(c)(f);
end
a=b;

if isempty(a.t)
  disp(sprintf('ded_fix_stats: %s: fixed stats file is empty',nm));
  unix(sprintf('/bin/rm -f %s',fns));
  return;
end

ded_write_stats(fns,a);
