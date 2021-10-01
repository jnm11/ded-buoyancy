function [f u]=ded_b_front_fun(x,f0,P,U,L)
% P u f_x =f_xx % Solve advection diffusion equation
% x is coordinate
% f0=f(0)
% P is peclet number
% U is velocity
% L is length scale
s=2*P*U*x.^(3/2)/(3*sqrt(L));

f=f0*gammainc(s,2/3,'lower')/gamma(2/3);

u=-U*sqrt(x./L);
%f=f0*exp(-P*(int u dx));
return;

x=linspace(0,10,1e3);
[f u]=ded_b_front_fun(x,1,1,1,1);
subplot(2,1,1);plot(x,f);
subplot(2,1,2);plot(x,u);

