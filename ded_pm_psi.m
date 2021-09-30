function [psiy psiz]=ded_pm_psi(d,y,z,h,W,H)

%[psiy psiz]=ded_pm_psi(0.6,0.2,0.3,1,3,4);disp(sprintf('%e ',[psiy psiz]))
%ju.periodic_gaussian(0.6,0.2,0.3,1,3,4)
if nargin==0
  test_ded_pm_psi();return;
end

Y=pi*2*y/W; % Scaled y coordinate;
Z=pi*2*z/H; % Scaled z coordinate;

d0=1/h;
Q0=h^2*pi;

w  = @(s) cos(s/2).^6.*(1+6*sin(s/2).^4+3*sin(s/2).^2);
WW = @(s) s/2+1/240*sin(s).*(9*sin(s).^4+20*sin(s).^2+120);

E   =  d.^2.*exp(-(H^2+W^2).* d.^2/2);
E0  = d0.^2.*exp(-(H^2+W^2).*d0.^2/2);
EY  = (WW(Y).*w(Z)-Y).*E;
EZ  = (WW(Z).*w(Y)-Z).*E;
HW0  = 3/8*pi^3*(H+W).*E0;
HW   = 3/8*pi^3*(H+W).*E;
dH0  = d0*H/sqrt(2)/pi;
dW0  = d0*W/sqrt(2)/pi;
dH   =  d*H/sqrt(2)/pi;
dW   =  d*W/sqrt(2)/pi;

GZ   = W/pi*d.^2.*ded_pm_psi_Fw(Y,dW).*w(Z).*exp(-2*z.^2.*d.^2);
GY   = H/pi*d.^2.*ded_pm_psi_Fw(Z,dH).*w(Y).*exp(-2*y.^2.*d.^2);

FF0 = H*W*d0.^2.*ded_pm_psi_Fw(pi,dH0).*ded_pm_psi_Fw(pi,dW0);
FF  = H*W* d.^2.*ded_pm_psi_Fw(pi,dH ).*ded_pm_psi_Fw(pi,dW );

C= Q0 * pi^2/4                 ./(HW0+FF0);
c= Q0/(H*W) * (HW0-HW + FF0-FF)./(HW0+FF0);

psiy = y.*c/2 +C.*(GZ-EY);
psiz = z.*c/2 +C.*(GY-EZ);

function test_ded_pm_psi()
Nx=500;
h=1;
dd=linspace(0,1/h,Nx);
W=20;
H=20;
yy=linspace(-W/2,W/2,201);
zz=linspace(-H/2,H/2,202);
e=1e-4;
for j=1:length(dd)
  d=dd(j);
  [y z]=ndgrid(yy,zz);
  [fy,fz]=ded_pm_psi(d,y,z,h,W,H);
  [fy1,fz1]=ded_pm_psi(d,y+e,z+0,h,W,H);
  [fy2,fz2]=ded_pm_psi(d,y-e,z-0,h,W,H);
  [fy3,fz3]=ded_pm_psi(d,y+0,z+e,h,W,H);
  [fy4,fz4]=ded_pm_psi(d,y-0,z-e,h,W,H);
  u=(fy1-fy2+fz3-fz4)/(2*e);
  %figure(j);set(gcf,'position',[48+5*j   511+5*j   600   1000]);
  subplot(3,1,1);mesh(zz,yy,u);title(sprintf('d=%5.3f int(u)=%5.3f',d,mean(u(:))*H*W));axis('tight');
  subplot(3,1,2);mesh(zz,yy,fy);title('fy');axis('tight');
  subplot(3,1,3);mesh(zz,yy,fz);title('fz');axis('tight');
  drawnow
end
  