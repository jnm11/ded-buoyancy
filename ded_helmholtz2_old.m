function [f g w d z]=ded_helmholtz2(a,nz)
% Return Helmholtz deomposition
% f is streamfunction
% g is potential
% w is vorticity
% d is divergence

if nargin==0
  for j=1:50
    ded_helmholtz2_test;
  end
  return;
end

if nargin<2
  nz=[];
end
sz=size(a.u);
dx=a.x(2)-a.x(1);
dz=a.z(2)-a.z(1);


dimz=1;
dimx=2;

M=prod(sz(2:end));
N=size(a.u,dimz);

fg=zeros(N,M);
ff=zeros(N,M);

wx = fft_modes(M,dx);

HS=(a.H/2)^2;
k=real(wx.^2);

% Fourier and Chebychev transform velocity
uc = fft(ichebf2c(reshape(a.u,N,M),dimz),[],dimx);
wc = fft(ichebf2c(reshape(a.w,N,M),dimz),[],dimx);

dudz=zeros(N,M);
dwdz=zeros(N,M);

% Calculate first derivatives
dudz(1:N-1,:) = ichebdiffc(uc,dimz,a.H);
dwdz(1:N-1,:) = ichebdiffc(wc,dimz,a.H);
dudx          = bsxfun(@times,wx,uc);
dwdx          = bsxfun(@times,wx,wc);

Iuc=ichebintc(uc,dimz,a.H);
Iwc=ichebintc(wc,dimz,a.H);

mu = diff(ichebc2f(Iuc,dimz,[-1 1]),1,dimz);
mv = diff(ichebc2f(Iwc,dimz,[-1 1]),1,dimz);

b1=[1 0 0];
b2=[1 0 0];

fd = dudx+dwdz;
fw = dwdx-dudz; 

mv(2:end)=0;
mu(2:end)=0;

for j=1:M
  fg(:,j) = icheb_helmholtz(HS*k(j), HS*fd(:,j),b1,[1 0 mv(j)]);
  ff(:,j) = icheb_helmholtz(HS*k(j),-HS*fw(:,j),b1,[1 0 mu(j)]);
end

if isempty(nz)
  z=[];
else
  z=linspace(-1,1,nz);
end

sz(1)=nz;

g=reshape(ichebc2f(ifft(fg,[],dimx,'symmetric'),dimz,z),sz);
f=reshape(ichebc2f(ifft(ff,[],dimx,'symmetric'),dimz,z),sz);
w=reshape(ichebc2f(ifft(fw,[],dimx,'symmetric'),dimz,z),sz);
d=reshape(ichebc2f(ifft(fd,[],dimx,'symmetric'),dimz,z),sz);

z=(1+z)*a.H/2;


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
a.u     =  dfz*fx;
a.w     = -fz*dfx;
a.omega = -fz*ddfx - ddfz*fx;
a.div   = 0; 
a.H=H;
a.L=L;
a.x=x;
a.z=z;

[f g w d]=ded_helmholtz2(a);
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

