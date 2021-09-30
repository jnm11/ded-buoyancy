function [a p]=ded_tune_PID(nm,T)
%ded_tune_PID('gc/f6/f/01')
if nargin<2
  T=[];
end
if isempty(T)
  T=0;
end

a=ded_read_stats(nm);
p=ded_read_param(nm);

if T<0
  T=max(a.t(abs(a.X-p.PIDX)>2));
  if isempty(T)
    T=0;
  end
  disp(sprintf('ded_tune_PID: T=%6.1f',T));
end

t=linspace(min(a.t),max(a.t),length(a.t))';
t(t<T)=[];
n=1e3;
if length(t)>n
  t=linspace(t(1),t(end),n);
end

dt=t(2)-t(1);

X=interp1(a.t,a.X,t,'linear');

if p.PIDG
  Y=interp1(a.t,a.g,t,'linear');
else
  Y=interp1(a.t,a.U,t,'linear');
end

X=X(:);
Y=Y(:);
s1 = iddata(X,Y,dt);

sys = idnlarx([2 2 1]);
%sys.Nonlinearity = 'sigmoidnet';
%sys.NonlinearRegressors = 'search';
sys = nlarx(s1,sys,'focus','sim');
clf;
compare(s1,sys);

%yf = forecast(sys,[X(1:2) Y(1:2)],length(Y)-2,Y(3:end));
%plot(yf);


lsys = linearize(sys,p.PIDX,repmat(Y(end),4,1));
opt = pidtuneOptions('PhaseMargin',70,'DesignFocus','reference-tracking');
[CPID,info]=pidtune(lsys,'PID',0.5,opt);
% 1 selects a fast response time
%TPID = feedback(CPID*lsys, 1);
%step(TPID);
%disp(CPID);
%disp(info);
disp(sprintf('X=%5.2f, --PIDP %8.5f --PIDI %8.5f --PIDD %8.5f',p.PIDX,CPID.Kp,CPID.Ki,-CPID.Kd));

clear a
a.PIDP=CPID.Kp;
a.PIDI=CPID.Ki;
a.PIDD=CPID.Kd;


