function [sys,b,fnm,s]=ded_gc_calc_larx(nm,ld,mint,fnm,ord)
%[sys,b,fnm,s]=ded_gc_calc_larx(nm,ld) Calculate LARX model relating X and U
if nargin<2
  ld=[];
end
if nargin<3
  mint=[];
end
if nargin<4
  fnm=[];
end
if nargin<5
  ord=[];
end
if isempty(ld)
  ld=0;
end
if isempty(ord)
  ord=[2 2 1];
end


if isempty(mint)
  mint=-inf;
end

if isempty(fnm)
  fnm=sprintf('~/%s/larx.mat',nm);
end

if isfile(fnm)
  load(fnm);
  if ld==1
    return;
  end
else
  sys=[];
end
B=struct();
figure;clf;
aa= jsubplot([1 3],[0.14 0.1],[0.01 0.01],[0.01 0.05]);
a1=aa(1);
a2=aa(2);
a3=aa(3);


  T=ded_convergence_T(nm);
  %a=ded_read_g(nm,'yz');
  %aa=ded_read_g(nm,'xyz');
  %aa.t=aa.t-mint;
  %maxt=min([max(a.t),max(aa.t)]);
  %mint=max([mmint,min(a.t),min(aa.t)]);
  %u=-aa.u/(p.H*p.L*p.W);
% $$$   a.t=a.t-mint;
% $$$   nt=floor(maxt/dt);
% $$$   b=struct('t',(0:nt)'*dt);
% $$$   b.b=make_monotonic(a.t,a.b,b.t);
% $$$   b.x=a.x;% $$$   b.U=make_monotonic(aa.t,u,b.t);
% $$$   c=ded_gc_find_front(b,T);
% $$$   b.X=c.X(:);
%u=a.g;
  a=ded_read_stats(nm);
  p=ded_read_param(nm);
  
  b.t=linspace(a.t(1),a.t(end),min(1e3,length(a.t)))';
  b.t=b.t(b.t>mint);
  b.Y=interp1(a.t,a.g,b.t,'linear');
  b.X=interp1(a.t,a.X,b.t,'linear');
  dt=b.t(2)-b.t(1);  
  if any(~isfinite(b.X)) | any(~isfinite(b.Y))
    error('ded_gc_calc_larx: Not finite');
  end
  s = iddata(b.X,b.Y,dt);
  axes(a1);
  plot(b.t,b.X);
  hold('on');
  axes(a2);
  hold('on');
  plot(b.t,b.Y);
  drawnow;
  %B.U(j)=b.U(end);
  %B.X(j)=b.X(end);
  
  
  
  axes(a1);
  if isfield(p,'PIDST')
    line(tt-mint,p.PIDS1 + (p.PIDS2-p.PIDS1)*floor(mod(tt/p.PIDST,2)),'color',.7*[1 1 1]);
  end
  
title(nm);
ylabel('X');
set([a1 a2],'xticklabel',[]);
axy=gu_setalim;
axes(a2);
ylabel('g');
gu_setalim;
axes(a3);
if isfield(p,'PIDST')
  line(tt-mint,p.PIDS1 + (p.PIDS2-p.PIDS1)*floor(mod(tt/p.PIDST,2)),'color',.7*[1 1 1]);
end
xlabel('t');

if false
  sys = idnlarx(ord);
  sys.Nonlinearity = 'sigmoidnet';
  sys.NonlinearRegressors = 'search';
  %sys = nlarx(s,sys,'focus','sim');
  sys = nlarx(s,sys);

else
  sys = procest(s,'P3UD');%sys = procest(data,type,'InputDelay',InputDelay)
end
disp(sys);
axes(a3);
compare(s,sys);
h=legend;
set(h,'location','best');
delete(title(''));
set(a3,'ylim',axy(3:4))
xlabel('t');
xlabel([]);



save(fnm,'sys','b')

return;

function nx=make_monotonic(t0,x,t1)
g=1:length(t0);
while 1
  t=t0(g);
  dt=median(diff(t));
  f=find(diff(t)<dt*.5);
  if isempty(f)
    break;
  end
  g(f)=[];
end
if size(x,1)>1 & size(x,2)>1
  nx=interp1(t,x(:,g)',t1,'linear','extrap')';
else
  nx=interp1(t,x(g),t1,'linear','extrap');
end


return;

ded_gc_calc_larx('2');
ded_gc_calc_larx('3');


yf = forecast(sys,[c.X(1:2) c.U(1:2)],length(c.U)-2,c.U(3:end));
plot(yf);
p=pidTuner(lsys);

lsys = linearize(sys,b.X(end),repmat(b.U(end),4,1));
pidtune(lsys,'PID');
