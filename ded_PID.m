nm='gc/f6/i/4000/0800/4608/xxxxc';
s=ded_read_stats(nm);
p=ded_read_param(mn);
a=ded_read_g(nm,'sev');

X=p.PIDX-s.X;
dt=diff(s.t);
dX=diff(X)./dt;
D(1)=a.PIDDD(1)/p.PIDD;
IX(1)=a.PIDIT(1)/p.PIDI;
for j=1:length(dX)
  dt=s.t(j+1)-s.t(j);
  e=exp(-dt);
  DX(j+1)=e*DX(j)+(1-e)*dX(j);
  IX(j+1)=IX(j)+dt*X(j);
end
plot(a.t,a.PIDDD,s.t,p.PIDD*DX);


Ix=cumsum(X.*dt);

Kp=p.PIDP*X(2:end);
Kd=p.PIDD*dX;
KI=a.PIDIT(1)+p.PIDI*Ix(2:end);

t=s.t(2:end);
h=plot(t,Kp+Kd+KI,s.t,s.g)
legend(h,'PID','g');

mg=mean(a.PIDIT);
clf;
plot(a.t,a.PIDIT-mg,a.t,p.PIDD*a.PIDDD,s.t,p.PIDP*X);

plot(a.t,a.PIDDD,p.PIDD*dX);

