function ded_calc_size(x,n,q)

N=n.*[1 1 2/3];

dx=x./N;
dx(~isfinite(dx))=[];
dx=exp(mean(log(dx)));

f=~isfinite(n);

nn=round(x/dx.*[1 1 3/2]/2)*2;

n(f)=nn(f);

M=8*36;


j=0:5;
k=(0:2)';
P=3.^k*2.^j;
P=-sort(-P(:));

m1=NaN;
m2=NaN;
for k=1:length(P)
  m1=P(k);
  if floor(n(2)/m1)*m1==n(2)
    break;
  end
end

m2=M/m1;

n(3)=m2*round(n(3)/m2);


disp(sprintf('--Nx %4i --Ny %3i --Nz %4i -W %6.4f --m1 %3i --m2 %3i # %6.0f dx=%8.6f',n(1),n(2),n(3),x(2),m1,m2,prod(n)/1e6,dx));
