%function [psi phi omega div]=ded_helmholtz_slip(a,H,dx)
%ded_helmholtz_slip(a) solve helmholtz equation with slip boundaries
%In this case the condition is that the vorticity omega vanishes on
%the boundaries

if 1
  nm='gc/f7/ma/00145/25125';
  c=ded_coord(nm);
  dx=c.dJx;
  p=ded_read_param(nm);
  H=p.H;
  fns=ded_get_fn(nm,'final',[],'state');
  a=ded_read_hdf(fns{2});
% $$$   dimz=1;
% $$$   nz=size(a.u,dimz);
% $$$   W=ichebintw(nz);
% $$$   Wc=ichebintwc(nz);
% $$$   uc=ichebf2c(a.u,dimz);
% $$$   plot(c.Jx,W*a.u,c.Jx,Wc*uc);
  
  x=c.Jx';
  z=c.Jz;
  % extend parity
  n=size(a.u,2);
  a.u = ded_parity(a.u,flip(a.u_parity));
  a.w = ded_parity(a.w,flip(a.w_parity));
  u=a.u;
  w=a.w;
end
dimx=2;
dimz=1;
e=[];
if 0
  clear('a');
  nz=10;
  nx=32;
  n=nx;
  H=1+2*rand(1);
  L=1+2*rand(1);
  H=5;
  p=8*pi/L;
  k=6*pi/L;
  z=ichebgrid(nz)*H/2;
  x=(0:nx-1)/nx*L;
  dx=L/nx;
  syms('X');
  syms('Z');
  fpsi = Z^3*(1-Z/H)^3*cos(p*X)+Z/H;
  fphi = Z^2*(1-Z/H)^2*cos(k*X);
  ps=randn(3,1);ph=randn(5,1);
  fpsi = (1-(2*Z/H)^2)*(ps(1)+ps(2)*Z/H+ps(3)*(Z/H)^2)*cos(p*X);
  fphi = 0*(1-(2*Z/H)^2)^2*(ph(1)+ph(2)*Z/H+ph(3)*(Z/H)^2+ph(4)*(Z/H)^3+ph(5)*(Z/H)^4)*cos(k*X);
  fdpsidz  = diff(fpsi,Z);fdpsidx  = diff(fpsi,X);
  fdphidz  = diff(fphi,Z);fdphidx  = diff(fphi,X);
  fdpsidzz = diff(fpsi,Z,Z);fdpsidxx  = diff(fpsi,X,X);fdpsidxz  = diff(fpsi,X,Z);
  fdphidzz = diff(fphi,Z,Z);fdphidxx  = diff(fphi,X,X);fdphidxz  = diff(fphi,X,Z);
  e.u   = double(subs(subs(fdphidx+fdpsidz,X,x),Z,z));
  e.w   = double(subs(subs(fdphidz-fdpsidx,X,x),Z,z));
  e.psi = double(subs(subs(fpsi,X,x),Z,z));
  e.phi = double(subs(subs(fphi,X,x),Z,z));
  e.dudx = double(subs(subs(fdphidxx+fdpsidxz,X,x),Z,z));
  e.dwdx = double(subs(subs(fdphidxz-fdpsidxx,X,x),Z,z));
  e.dudz = double(subs(subs(fdphidxz+fdpsidzz,X,x),Z,z));
  e.dwdz = double(subs(subs(fdphidzz-fdpsidxz,X,x),Z,z));
  e.div  = e.dudx+e.dwdz;
  
  a.u = +pr_diff(e.phi,dx,dimx) + ichebdifff(e.psi,dimz,H);
  a.w = -pr_diff(e.psi,dx,dimx) + ichebdifff(e.phi,dimz,H);
  a.dudx = pr_diff(a.u,dx,dimx);
  a.dwdx = pr_diff(a.w,dx,dimx);
  a.dudz = ichebdifff(a.u,dimz,H);
  a.dwdz = ichebdifff(a.w,dimz,H);
  a.div  = a.dudx+a.dwdz;
  a.uc   = fft(ichebf2c(a.u,dimz),[],dimx);
  a.wc   = fft(ichebf2c(a.u,dimz),[],dimx);
  
  E0=(-1).^(0:nz-1);
  E1=ones(1,nz);
  r = [max(max(abs(e.u-a.u))) max(max(abs(e.w-a.w)))];
  %r(3) = max(abs(E0*a.wc));
  %r(4) = max(abs(E1*a.wc));
  %disp(r);
  u=a.u;
  w=a.w;
end



nx = size(u,dimx);
nz = size(u,dimz);

k  = fft_modes(nx)/dx;
kk = k.^2;

% Chebychev transform;
u = ichebf2c(u,dimz);
w = ichebf2c(w,dimz);

W=ichebintw(nz)*H/2;
Wc=ichebintwc(nz)*H/2;
q=Wc*u; 
x=[x x+p.L];
plot(x,q,x,W*a.u);
keyboard;


% Fourier transform;
u = fft(u,[],dimx);
w = fft(w,[],dimx);


ux=k.*u;
wx=k.*w;
uz=ichebdiffc(u,dimz,H,1,true);
wz=ichebdiffc(w,dimz,H,1,true);




div   = ux+wz;
omega = wx-uz;

% $$$ E0=(-1).^(0:nz-1);
% $$$ E1=ones(1,nz);
% $$$ 

% $$$ 
% $$$ phic=zeros(size(w));
% $$$ psic=zeros(size(u));
% $$$ phic(:,1)=ichebintc(w(:,1),dimz,H,[],true);
% $$$ psic(:,1)=ichebintc(u(:,1),dimz,H,[],true);
% $$$ 
% $$$ f=find(kk~=0);
% $$$ phic(:,f)=icheb_helmholtz0(kk(f)/(2/H)^(2),  +div(:,f)/(2/H)^(2));
% $$$ psic(:,f)=icheb_helmholtz0(kk(f)/(2/H)^(2),-omega(:,f)/(2/H)^(2));

Iz  = sparse(ichebintc(eye(nz),dimz,H,1,true));
Dz  = ichebD(nz,H);
DD  = Dz*Dz;
II  = Iz*Iz;
N=0:nz-1;
E0  = (-1).^N;
E1  = ones(1,nz);
D0  = E0.*N.^2;  % First derivative at z=0
D1  = E1.*N.^2;  % First derivative at z=H
DD0 = E0.*(N-1).*N.^2.*(N+1)/3; % Second derivative at z=0
DD1 = E1.*(N-1).*N.^2.*(N+1)/3; % Second derivative at z=H

psic(:,1)=Iz*u(:,1);
phic(:,1)=Iz*w(:,1);

f=find(kk~=0);
psic(:,f)=icheb_helmholtz0(kk(f)/(2/H)^(2),-omega(:,f)/(2/H)^(2));
I=eye(nz); % Must use dense matrix solve ! sparse solve break;
Dz=full(Dz);
%Constant and linear terms are very close to degenerate
for j=2:length(k)
  %psic(:,j) = [DD+kk(j)*eye(nz);DD0+E0*kk(j);DD1+E1*kk(j)]\[Dz*u(:,j)-k(j)*w(:,j);0;0];
  %psic(:,j) = [Dz+kk(j)*Iz;DD0+E0*kk(j);DD1+E1*kk(j)]\[u(:,j)-k(j)*Iz*w(:,j);0;0];
  phic(:,j) = (u(:,j)-Dz*psic(:,j))/k(j);
end
A=[DD+kk(j)*eye(nz);DD0+E0*kk(j);DD1+E1*kk(j)];
Z=null(A);
disp(Z(1:5,:));

n=nx;

% $$$ Iuz = ichebintc(u,dimz,H,1,true);
% $$$ Iwz = ichebintc(w,dimz,H,1,true);
% $$$ phic=0*Iuz;
% $$$ psic = Iuz - ichebintc(k.*phic,dimz,H,1,true);
% $$$ phic = Iwz + ichebintc(k.*psic,dimz,H,1,true);

[phic psic ee]=icheb_phi_psi(k*H/2,u*H/2,w*H/2);disp(ee);

%[phi2 psi2]=icheb_phi_psi(k*H/2,u,w);
%keyboard


ddiv=ichebdiffc(phic,dimz,H,2,true)+kk.*phic;

uu = +k.*phic + ichebdiffc(psic,dimz,H,1,true);
ww = -k.*psic + ichebdiffc(phic,dimz,H,1,true);

uuz=ichebdiffc(uu,dimz,H,1,true);
%disp([max(abs(E0*ww)) max(abs(E1*ww)) max(abs(E0*uuz)) max(abs(E1*uuz))]);

%keyboard

ddiv  = ifft( ddiv, [], dimx, 'symmetric');
div   = ifft(  div, [], dimx, 'symmetric');
omega = ifft(omega, [], dimx, 'symmetric');
phi   = ifft(  phic, [], dimx, 'symmetric');
psi   = ifft(  psic, [], dimx, 'symmetric');
uu    = ifft(   uu, [], dimx, 'symmetric');
ww    = ifft(   ww, [], dimx, 'symmetric');


ddiv  = ichebc2f( ddiv(:,1:n), dimz, [], H);
div   = ichebc2f(  div(:,1:n), dimz, [], H);
omega = ichebc2f(omega(:,1:n), dimz, [], H);
phi   = ichebc2f(  phi(:,1:n), dimz, [], H);
psi   = ichebc2f(  psi(:,1:n), dimz, [], H);  
uu    = ichebc2f(   uu(:,1:n), dimz, [], H);
ww    = ichebc2f(   ww(:,1:n), dimz, [], H);  

if false
  W=ichebintw(nz)/2;
  p=polyfit(x,W*phi,1);
  phi=phi-polyval(p,x);
  psi=psi+p(1)*z;
end
if false
  W=ichebintw(nz)/2;
  mp=W*phi;mp=mp-mean(mp);
  phi=phi-mean(phi(:));
  cc=cos(x*pi/p.L);
  phi=phi-cc*(sum(cc.*mp)/sum(cc.^2));
  p=polyfit(x,W*phi,1);
  phi=phi-polyval(p,x);
  psi=psi+p(1)*z;
end

figure(1);clf;
subplot(3,2,1);mesh(x,z,uu-a.u);  title('du');  xlabel('x');
subplot(3,2,2);mesh(x,z,ww-a.w);  title('dw');  xlabel('x');
subplot(3,2,3);mesh(x,z,phi);title('phi');      xlabel('x');
subplot(3,2,4);mesh(x,z,psi);title('psi');      xlabel('x');
subplot(3,2,5);mesh(x,z,div);title('div');      xlabel('x');
subplot(3,2,6);mesh(x,z,omega);title('omega');  xlabel('x');
