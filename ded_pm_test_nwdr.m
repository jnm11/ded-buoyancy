function ded_pm_test_nwdr(fn)
if nargin==0
  ded_pm_test_nwdr('pm/test/03','nwdr');
  pause;
  ded_pm_test_nwdr('pm/test/03','inwdr');
  ded_pm_test_nwdr('pm/S2/01-16','nwdr');%fn='pm/S2/01-16';nm='nwdr';
  return;
end

  
a=ded_read_hdf([ded_dedalus_data_dir '/' fn '/force/fyz.hdf5']);
p=ded_read_param(fn);
c=ded_coord(fn);
y=c.Ay;
z=c.Az;
dy=c.dAy;
dz=c.dAz;
dA=dy*dz;

f=a.nwdr;

%Take account of offset grid by replicating
IA=1;
if strcmp(p.Ty,'Fourier')
  y(end+1)=y(1)+p.W;
  f(end+1,:)=f(1,:);
  f([1 end],:)=f([1 end],:)/2;
else
  IA=IA*2;
  y=[-flip(y);y];
  f=cat(1,flip(f,1),f);
end
if strcmp(p.Tz,'Fourier')
  IA=IA*2;
  z(end+1)=z(1)+p.H;
  f(:,end+1)=f(:,1);
  f(:,[1 end])=f(:,[1 end])/2;
else
  z=[-flip(z);z];
  f=cat(2,flip(f,2),f);
end

% Midpoint integration weight
%w=repmat(dA,size(f));
%([1 end],:)=w([1 end],:)/2;
%w(:,[1 end])=w(:,[1 end])/2;
%sum(w(:))-p.H*p.W  



N=20;
A=zeros(N,1);
[yy zz]=ndgrid(y,z);
rr=sqrt(yy.^2+zz.^2);
q = atan2(yy,zz);

A=zeros(N+1,1);
for n=0:N
  A(n+1) = sum(sum(exp(n*i*q).*f*dA));
end
A=A/IA;
A(1)=A(1)-1;

dc;
figure(1)
mesh(y,z,f);
title(sprintf('A=%f',A(1)));
[yy zz]=ndgrid(y,z);
figure(2)
plot(1:N,abs(A(2:end)).^2);
return;

p.W=20;
p.H=10;
r=2;
HW=p.W/2;
HH=p.H/2;
R=sqrt(HW^2+HH^2);
N=200;
M=round(N*HH/HW);
y=HW*(-N:2:N)/(N+1);
z=HH*(-M:2:M)/(M+1);
[yy zz]=ndgrid(y,z);
qq=atan2(yy,zz);
rr=sqrt(yy.^2+zz.^2);
dy=y(2)-y(1);
dz=z(2)-z(1);
dA=dy*dz;
disp(dA*length(y)*length(z));



dr=2*sqrt(dy^2+dz^2);
f=interp1([0 r-dr/2 r+dr/2 R],[0 0 1 1],rr,'pchip')./(min(HW./abs(sin(qq)),HH./abs(cos(qq)))-r);
I1=r*pi+2*HH*asinh(HW/HH)+2*HW*asinh(HH/HW);
I2=dA*sum(f(:));
disp([I1 I2]);


f=f/(dA*sum(f(:)));





A=reshape(A,Ny*Nz,2*n);
c=A\f(:);
g=reshape(f(:)-A*c,Ny,Nz);
figure(2)
mesh(g)
%f=f/(dA*sum(f(:)));
disp(reshape(c,1,2*n));
disp(reshape(A'*f(:),1,2*n));
disp(reshape(A'*g(:),1,2*n));


% $$$ D=zeros(n,n);
% $$$ for j=1:n
% $$$   for k=1:n
% $$$     D(j,k)=dA*sum(sum((sin(j*q).*cos(k*q)).^2));
% $$$   end
% $$$ end
% $$$ 
% $$$ 
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ N=zeros(n,n,n,n);
% $$$ for j1=1:n
% $$$   for k1=1:n
% $$$     for j2=1:n
% $$$       for k2=1:n
% $$$         N(j1,k1,j2,k2)=dA*sum(sum(sin(j1*q).*cos(k1*q).*sin(j2*q).*cos(k2*q)));
% $$$       end
% $$$     end
% $$$   end
% $$$ end
% $$$ N=reshape(N,n*n,n*n);
% $$$ A=A(:);
% $$$ c=reshape(A\N,n,n);
% $$$ for j=1:n
% $$$   for k=1:n
% $$$     f=f-c(j,k)*sin(j*q).*cos(k*q);
% $$$   end
% $$$ end
% $$$ mesh(f);
% $$$ 
% $$$ 
