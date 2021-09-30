
Ny=31;
Nz=52;
H=5;
W=6;
y=linspace(-W/2,W/2,Ny);
z=linspace(-H/2,H/2,Nz);
f=randn(Ny,Nz);

dA=W/Ny*H/Nz;

[yy zz]=ndgrid(y,z);
q=atan2(yy,zz);
n=8;
A=zeros(Ny,Nz,n,n);
for j=1:n
  for k=1:n
    A(:,:,j,k)=sin(j*q).*cos(k*q);
  end
end
A=reshape(A,Ny*Nz,n*n);
f=f(:);
c=A\f;
g=f-A*c;
g=reshape(g,Ny,Nz);
mesh(g)
f=f/(dA*sum(f(:)));
%/Users/jnm/smop/smop-master/smop/main.py /Users/jnm/matlab/dedalus/ded_pm_remove_q.m
