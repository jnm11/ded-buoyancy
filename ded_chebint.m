m=12;
n=64;
p=randn(1,m);
p(end+1:n)=0;

clear all
n=64;
T=@(x,k) cos(k*acos(x));
z=linspace(0,1,1e3);
for k=1:n-1
  an=zeros(n,1);
  an(k+1)=1;
  f1 = cheval('regular',an,z);
  f2 = T(z,k);
  plot(z,f1,z,f2);
  e(k)=sqrt(mean((f1-f2).^2));
end
disp(max(e));


t=2*(((1-n)/2):((n-1)/2))'/n;
z=sin(pi*t/2);


tt=[t-1;t+1];
zz=[z(end:-1:1);z];
disp(max(abs(zz-sin(pi*tt/2))));

p1=randn(n,1);
f1=cheval('regular',p1,z);
p2 = mdct(f1).*(-1).^(0:n-1)';
disp(max(abs(p1-p2)));

return;


T=@(x,k) cos(k*acos(x));
k=round(20*rand(1));
dx=1e-5;
x=linspace(-1+dx/2,1-dx/2,1e3);
x=x(2:end-1);
f=T(x,k);
F=@(x)(T(x,k+1)/(1+k)-T(x,k-1)/(k-1))/2;
dF=(F(x+dx/2)-F(x-dx/2))/dx;
clf
plot(x,f,x,dF)
disp(max(abs(f-dF)));
title(max(abs(f-dF)));








da=polyval(p(1:n-1).*(n-1:-1:1)',z);
Ia=polyval([p./(n:-1:1)';0],z);

p=[randn(m,1);zeros(n-m,1)];
a=0;
for j=0:length(p)-1
  a=a+T(z,j)*p(j+1);
end





if 0
  s=tt*pi/2;
  ss=linspace(-pi,pi,1e4);
  for k=0:5
    a=T(z,k);
    aa=[a(end:-1:1);a];
    bb=(-1)^k*cos(k*ss);
    clf;
    subplot(2,1,1);
    plot(s,aa,s,(-1)^k*cos(k*s),ss,bb);
    subplot(2,1,2);
    plot(s,aa-(-1)^k*cos(k*s));
    axis('tight');
    drawnow;
    pause;
  end
end



f=fft(aa(1:end-1));disp(max(abs(imag(f))));
f=fft(aa(2:end));disp(max(abs(imag(f))));
f=fft(aa,2*n+1);disp(max(abs(imag(f))));


a=T(z,4);
aa=[a(end:-1:1);a];
f=real(fft(aa,2*n-1)));
c=real(ifft(f));
c=c(n:-1:1);
f=f(1:n)/n;
f(1)=f(1)/2;
b=0;
for j=0:n-1
  b=b+T(z,j)*f(j+1);
end
plot(z,a,z,c,z,b);



figure(1);
subplot(2,1,1);
plot(z,a,z,b);
subplot(2,1,2);
plot(z,a-b);
axis('tight');


for j=0:n-1
  for k=0:n-1
    A(j+1,k+1)=2*T(z,j)*T(z,k)'/n;
  end
end








