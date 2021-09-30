function [a p]=ded_model_est(nm,T)
%ded_tune_PID('gc/f6/f/01')

if nargin==0
  for j=1:31
    fns{j}=sprintf('gc/gc2d7m/%02u',j);
  end
  [a p]=ded_model_est(fns);
  return;
end
if nargin<2
  T=[];
end
if isempty(T)
  T=0;
end
if ~iscell(nm)
  nm={nm};
end
dt=0.1;
for j=1:length(nm)
  [z{j} p(j)]=ded_model_data(nm{j},T,dt);
end
z=merge(z{:});

V = arxstruc(z,z,struc(1:6,1:6,1:6)); % try values in the range 1:5 for na, nb, nk%

Order = selstruc(V,'aic');
%Order=[4 4 4 4];
s1=sprintf('ded_model_est: order [%u %u %u], ',Order(1),Order(2),Order(3));
% $$$ advice(z, 'nonlinearity')

%sys = armax(z,Order);
sys = arx(z,Order);

figure(1);clf;
compare(z, sys)

opt = pidtuneOptions('PhaseMargin',70,'DesignFocus','reference-tracking');
[CPID,info]=pidtune(sys,'PID',.5,opt);
% 1 selects a fast response time
%TPID = feedback(CPID*lsys, 1);
%step(TPID);
%disp(CPID);
%disp(info);
s2=sprintf('X=%5.2f, --PIDP %8.5f --PIDI %8.5f --PIDD %8.5f',p(1).PIDX,CPID.Kp,CPID.Ki,CPID.Kd);
disp([s1 s2]);
clear a
a.PIDP=CPID.Kp;
a.PIDI=CPID.Ki;
a.PIDD=CPID.Kd;


function [z,p] = ded_model_data(nm,T,dt)
p=ded_read_param(nm);
if isempty(p) 
  disp(sprintf('ded_model_est: %s empty param',nm));
  return;
end
a=ded_read_stats(nm);
if isempty(a) 
  disp(sprintf('ded_model_est: %s empty stats',nm));
  return;
end

if T<0
  T=max(a.t(abs(a.X-p.PIDX)>2));
  if isempty(T)
    T=0;
  end
  disp(sprintf('ded_tune_PID: T=%6.1f',T));
end



t=min(max(T,a.t)):dt:max(a.t);

% $$$ t=linspace(min(a.t),max(a.t),length(a.t))';
% $$$ t(t<T)=[];
% $$$ n=1e3;
% $$$ if length(t)>n
% $$$   t=linspace(t(1),t(end),n);
% $$$ end
% $$$ 
% $$$ dt=t(2)-t(1);
u=interp1(a.t,a.U,t,'linear');
g=interp1(a.t,a.g,t,'linear');
X=interp1(a.t,a.X,t,'linear');

if max(g)==min(g)
  Y=u;
else
  Y=g;
end

X=X(:);
Y=Y(:);
rg1=1:round(length(X)/2);
rg2=round(length(X)/2):length(X);
%z1 = iddata(X(rg1),Y(rg1),dt);
%z2 = iddata(X(rg2),Y(rg2),dt);
z  = iddata(X,Y,dt);
%
