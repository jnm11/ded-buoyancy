function ded_gc_front_exp(x,b)
load('~/gc/f6/g/1000/1500/2304/23400/profile.mat');
f1=max(find(a.u0>0));
f2=min(find(a.u0(f1:end)<=.9*min(a.u0)));
rg=f1+(-f2:f2;
b=a.b0(rg);
u=a.u0(rg);
x=a.x(rg);
xdx=x(2)-x(1);

y=cumsum([0 u])*dx;
y=y-max(y);
f=@(p) dx*flip(cumsum(flip(exp(1000*p(2)*(y(1:end-1)-p(1))))));
ff=@(p) f(p)-b;
p=lsqnonlin(ff,[0 1]);
plot(x,f(p),x,b,'s');


plot(x,u,x,b);

