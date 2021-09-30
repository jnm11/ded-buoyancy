function [f g w d z]=ded_helmholtz2(u,v,L,H,Px,nz)
% Return Helmholtz deomposition
% f is streamfunction
% g is potential
% w is vorticity
% d is divergence
%
%  u = g_x+f_y
%  v = g_y-f_z
%
% if k is wavenumver
% boundary conditions are 
% 
%  kf+ig' = 0 
% ikg +f' = 0 no-slip
% ikg'+f''= 0 slip
%
if nargin==0
  for j=1:50
    ded_helmholtz2_test;
  end
  return;
end

if nargin<3
  nz=[];
end


dimz=1;
dimx=2;

sz=size(u);
dx=L/sz(dimx);


M=prod(sz(2:end));
N=size(u,dimz);

fg=zeros(N,M);
ff=zeros(N,M);

wx = fft_modes(M,dx);


HS=(H/2)^2;
k=real(wx.^2);

% Fourier and Chebychev transform velocity
uc = fft(ichebf2c(reshape(u,N,M),dimz),[],dimx);
wc = fft(ichebf2c(reshape(v,N,M),dimz),[],dimx);


dudz=zeros(N,M);
dwdz=zeros(N,M);

% Calculate first derivatives
dudz(1:N-1,:) = ichebdiffc(uc,dimz,H);
dwdz(1:N-1,:) = ichebdiffc(wc,dimz,H);
dudx          = bsxfun(@times,wx,uc);
dwdx          = bsxfun(@times,wx,wc);

fd = dudx+dwdz;
fw = dwdx-dudz; 

for j=1:M
  disp(j);
  rw=imag(wx(j));
  a=zeros(4,6,2);
  a(1,:,1)=[rw 0 0 0 0 0];a(1,:,2)=[0 i   0 0   0 0]; % rwf+ig' = 0  z=-1
  a(2,:,1)=[0 0 0 rw 0 0];a(2,:,2)=[0 0   0 0   i 0]; % rwf+ig' = 0  z=1
  if 0  % no slip
    a(3,:,1)=[0 1 0 0 0 0];a(1,:,2)=[i*rw 0 0   0 0 0]; % rwf+ig' = 0  z=-1
    a(4,:,1)=[0 0 0 0 1 0];a(2,:,2)=[0   0 0 i*rw 0 0]; % rwf+ig' = 0  z=1
  else %no slip
    a(3,:,1)=[0 0 1 0 0 0];a(1,:,2)=[i*rw 0 0   0 0 0]; % rwf+ig' = 0  z=-1
    a(4,:,1)=[0 0 0 0 0 1];a(2,:,2)=[0   0 0 i*rw 0 0]; % rwf+ig' = 0  z=1
  end
  if j==1
    a(end+1,1,1)=1;
  end
  
  xx = icheb_helmholtz_double(HS*k(j), [fw(:,j) fd(:,j)],a,zeros(1,size(a,1)));
  ff(:,j) = xx(:,1);
  fg(:,j) = xx(:,2);
end

if isempty(nz)
  z=[];
else
  z=linspace(-1,1,nz);
  sz(1)=nz;
end



g=reshape(ichebc2f(ifft(fg,[],dimx,'symmetric'),dimz,z),sz);
f=reshape(ichebc2f(ifft(ff,[],dimx,'symmetric'),dimz,z),sz);
w=reshape(ichebc2f(ifft(fw,[],dimx,'symmetric'),dimz,z),sz);
d=reshape(ichebc2f(ifft(fd,[],dimx,'symmetric'),dimz,z),sz);

z=(1+z)*H/2;


return
function ded_helmholtz2_test

M=round(7+5*rand(1));
N=round(7+5*rand(1));

L=1+rand(1);
H=1+rand(1);

k=floor(M/2*rand(1))*2*pi/L;

x=(0:M-1)/M*L;
z=(1+ichebgrid(N))/2*H;

p=[randn(N-1,1);0];
%p=[p;0]-[0;p];

dp   = poly_diff(p,1);
ddp  = poly_diff(dp,1);

C=randn(1)+i*randn(1);

fx   = real(     C*exp(i*k*x));
dfx  = real( i*k*C*exp(i*k*x));
ddfx = real(-k^2*C*exp(i*k*x));

fz   = polyval(  p,z);
dfz  = polyval( dp,z);
ddfz = polyval(ddp,z);

a.psi   =  fz*fx;
a.phi   =  0;
u     =  dfz*fx;
w     = -fz*dfx;
a.omega = -fz*ddfx - ddfz*fx;
a.div   = 0; 
H=H;
L=L;

[f g w d]=ded_helmholtz2(u,w,L,H,0,[]);
e1=sqrt(sum((a.psi(:)-f(:)).^2));
e2=sqrt(sum((a.phi(:)-g(:)).^2));
e3=sqrt(sum((a.omega(:)-w(:)).^2));
disp(sprintf('%2u %2u %8.2e %8.2e %8.2e',M,N,e1,e2,e3));



return;
      
sum(sum(f.*a.psi))/sqrt(sum(sum(a.psi.^2))*sum(sum(f.^2)))
sum(sum(w.*a.omega))/sqrt(sum(sum(a.omega.^2))*sum(sum(w.^2)))

%plot(sum(abs(ff).^2,1),'s');

fx=fx/sqrt(sum(fx.^2));
fz2=f*fx';

cz=(fz'*fz2)/sqrt((fz'*fz)*(fz2'*fz2))
cz=(fz'*fz2)/(fz'*fz)

plot(z,fz2,z,fz);

zz=linspace(0,H,100);
plot(zz,polyval(p,zz),z,fz,'s',z,f);

