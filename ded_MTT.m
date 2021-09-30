function [B Q M F C r]=ded_MTT(x,F0,g,a,b)
%ded_MTT(x,F0,g,a) Morton-Turner-Taylor Plume
% F is buoyancy flux
% g is gravity 
% a is entrainment coefficient alpha
% b is shape factor for density relative to momentum

x=max(x,eps);

r =  6/5       *a      *x;
if 0
  rho = 5/6*(5*b^2*(b^2+1)^2*(1/36)*F0^2/(g*a^4))^(1/3);
  U   = sqrt(6)*sqrt(g*rho)/(2*b);
  B = 36/25*rho^1*a^2/b^2*x.^(+1/3);
  C = 36/25*rho^2*a^2/b^2*x.^(-4/3);
  Q = 36/25*U^1  *a^2    *x.^(+5/3);
  M = 36/25*U^2  *a^2    *x.^(+4/3);
  F = 72/25*U*rho*a^2/(b^2+1);
else
 B = (  12/25 *(b^2+1)^2*F0^2/g  *a^2/b^4.*x   ).^(1/3);
 C = (  25/324*(b^2+1)^4*F0^4/g^2/a^2/b^2./x.^4).^(1/3);
 M = (  81/400*(b^2+1)^2*F0^2*g^2*a^2/b^4.*x.^4).^(1/3);
 Q = ( 486/625*(b^2+1)  *F0  *g  *a^4/b^2.*x.^5).^(1/3);
 F = F0;
end

  
if abs(F-F0)>1e-10
  error('F does not match');
end
F=repmat(F0,size(x));

if nargout==1;
  B=[B(:);Q(:);M(:);F(:)];
end

return;

x=linspace(1,10,1e3);
F0=rand(1);
g=rand(1);
a=rand(1);
b=rand(1);
[B Q M F BB r]=ded_MTT(x,F0,g,a,b);

plot(x,B,x,Q,x,M,x,F,x,BB,x,r);






