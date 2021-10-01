function ded_fit_PID(nm,m);
%ded_fit_PID('gc/f6/f/33');
p=ded_read_param(nm);
s=ded_read_stats(nm);

x=s.g-mean(s.g);
f=find(x(1:end-1).*x(2:end)<=0);
TC=(s.t(f)+s.t(f+1))/2;
dt=diff(s.t(f));
P=mean(dt(dt>median(dt)/2));
TS=100;
p=fit_exp(s.t/TS,abs(s.g-mean(s.g)));

p=[p(1) pi/P*TS];
ff=find(s.t>0);
opt=optimset('display','none');
for k=1:m
  f=@(p) fit_fun_shm_exp(s.t(ff),s.g(ff),p/TS,k);
  p=lsqnonlin(f,p,[0 0],[],opt);
end
[r y c]=fit_fun_shm_exp(s.t(ff),s.g(ff),p/TS,m);
plot(s.t,s.g,s.t(ff),y,s.t([1 end]),c([1 1]));
T=TS/p(1);
P=TS*2*pi/p(2);
title(sprintf('g=%6.4f T=%5.1f P=%5.1f',c(1),T,P));

return;

