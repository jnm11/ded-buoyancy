function u=ded_uerf(y,h,w,dU,U,H,C,D)
%u=ded_uerf(y,h,w,dU,U,H,n) Velocity profile centred on h, width 1/w ni region y [0,H]
% du/dy=0 when y=0,H
% int u du  = H*U
% u(H)-u(0) = U
erfy = erf(w*(y-h));    
erfh = h/H*erf(w*h);       
erfH = (H-h)/H*erf((H-h)*w);    
exph = exp(-w^2*h^2)/sqrt(pi);   
expH = exp(-w^2*(H-h)^2)/sqrt(pi); 

B=dU/((-w*H*exp(-w^2*(H-h)^2)-exp(-w^2*h^2)*w*H+(erf((H-h)*w)+erf(w*h))*sqrt(pi))/sqrt(pi));

u=U+B*(erfh+erfy-erfH)+B*(exph-expH)/H/w+w*B/(3*H)*(exph*(2*H^2-6*H*y+3*y.^2)+expH*(H^2-3*y.^2));

if nargin>6
  u=u+C.*(-(-3./4+(H.^2-3.*H.*h+3./2.*(h.^2)+3./2.*(y.^2)).*w.^2).*erfH./((3.*(H-h)).*w.^2)-(-3./4+(H.^2-3.*H.*y+3./2.*(h.^2)+3./2.*(y.^2)).*w.^2).*erfh./(3.*h.*w.^2)-erfy.*(-y+h)-2.*h.*exph.*(3./4+(H.^2-3.*H.*y+3./2.*(y.^2)).*w.^2)./(3.*w.*H)+(H-h).*expH.*(2.*H.^2.*w.^2-6.*w.^2.*y.^2-3)./(6.*w.*H));
end

if nargin>7
  u=u+D.*((1./3).*(2.*H.*h-h.^2-3.*y.^2).*erfH+(1./3).*(2.*H.^2-6.*H.*y+h.^2+3.*y.^2).*erfh+(-y+h).^2.*erfy+(2.*(1./2+h.^2.*(-3.*H.*y+3./2.*(y.^2)+H.^2).*w.^4+(1./2).*w.^2.*h.^2)).*exph./(3.*H.*w.^3)+(-1+(H.^2-3.*y.^2).*(H-h).^2.*w.^4-w.^2.*(H-h).^2).*expH./(3.*H.*w.^3));
end

return;

N=100;
H=2*rand(1);
h=rand(1)*H;
w=10/H*(1+rand(1));
dU=randn(1);
U=randn(1);
n=0;
y=linspace(0,H,N);
dy=1e-5;
u=ded_uerf(y,h,w,dU,U,H,n);
du=(ded_uerf(y+dy,h,w,dU,U,H,n)-ded_uerf(y-dy,h,w,dU,U,H,n))/(2*dy);
subplot(1,2,1);plot(u,y/H);
subplot(1,2,2);plot(du,y/H);line([0 0],[0 1]);

U1=sum(u(1:end-1)+u(2:end))/(2*(N-1));
disp(U1/U-1);disp((u(end)-u(1))/dU-1);
