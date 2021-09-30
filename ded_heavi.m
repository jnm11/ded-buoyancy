


n=1024;
x=linspace(0,2*pi,n);
y=sin(x)+1/3*sin(3*x);
plot(x,y);


L=16;
n=512;
x=(0:n-1)/n*L;
dx=x(2)-x(1);

x1=4;
x2=12;

D=n/L;
DD=n/L*logspace(-3,1,100);
DD=logspace(0,1,1000);
for j=1:length(DD)
  D=DD(j);
  y = (erf((x-x1)*D)-erf((x-x2)*D))/2;
  f=fft(y);
  w=fftfreq(n);
  wc=2/3*pi;
  f(abs(w)>=wc)=0;
  ny=ifft(f);
  %  plot(x,y,x,ny);
  e1(j)=max(abs(y-ny));
  e2(j)=sqrt(mean((y-ny).^2));
end
loglog(DD,e1,DD,e2);
axis('tight');





y=sign(abs(x)-1/2);
plot(x,y);

f=fft(y);
w=fftfreq(n);
wc=pi/5;
f(abs(w)>=wc)=0;
ny=ifft(f);

2y3=y;
for j=1:1000
  f=fft(y3);
  f(abs(w)>=wc)=0;
  y4=ifft(f);
  y3=min(1,max(-1,[-sort(-y4(x<0)) sort(y4(x>=0))]));
  plot(x,y,x,y3,x,y4);
  title(j)
  drawnow;
  if max(abs(y4-y3))<1e-3
    break;
  end
end
