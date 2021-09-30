function ded_tanh(typ)


% The function all interpolate between 0 and 1 and have a gradient of  1 at 0
% pchip 25
% erf  2.6--2.7
% tanh 2.7--3.1
% p7   4.5
% p11  2.9
% p11 looks like a good choice similar spectral bandwidth but more compact
% if we choose 3 dx for the multiplier will be very narrow indeed
w=1;

L=1;
L1=L/4;
L2=3*L/4;
p7=[-5 0 21 0 -35 0 35 16]/32;c7=32/35;
p11=[-63 0 385 0 -990 0 1386 0 -1155 0 693 256]/512;c11=512/693;

switch(typ)
  case 'pchip'
    ff = @(x) interp1([-1e6 -3/4 3/4 1e6],[0 0 1 1],x,'pchip'); % Gradient is 1
  case 'erf'
    ff = @(x) (1+erf(x*sqrt(pi)))/2;
  case 'tanh'
    ff = @(x) (1+tanh(2*x))/2;
  case '7'
    ff = @(x) polyval(p7,max(-1,min(1,c7*x)));
  case '11'
    ff = @(x) polyval(p11,max(-1,min(1,c11*x)));
end

nn=round(2.^(-1+linspace(4,12,9+32)))*2;
for k=1:length(nn)
  a.n=nn(k);
  a.m=a.n/2;
  a.N=2^15;
  a.y=zeros(1,a.N);
  a.x=linspace(0,L,a.n+1);a.x(end)=[];
  a.X=linspace(0,L,a.N+1);a.X(end)=[];
  dx=a.x(2)-a.x(1);
  a.f = @(ww) ff(ww/dx*(a.x-L1))-ff(ww/dx*(a.x-L2));
  
  w1=0.1;
  w2=10;
  
  
  for j=1:6
    e1=ef(w1,a);
    if e1<0
      break;
    end
    w1=w1/10;
  end
  for j=1:6
    e2=ef(w2,a);
    if e2>0
      break;
    end
    w2=w2*10;
  end
  
  if e1*e2>0
    error('Starting points do not change sign')
  end
  
  for j=1:1000
    dw=w2-w1;
    w=(e2*w1-e1*w2)/(e2-e1);
    w=max(w1+1e-2*dw,w);
    w=min(w2-1e-2*dw,w);
    e=ef(w,a);
    disp(sprintf('w=[%8.1e %5.3f %8.1e], e=[%8.1e %8.1e %8.1e]',w1-w,w,w2-w,e1,e,e2));
    if e<=0
      w1=w;
      e1=e;
    else
      w2=w;
      e2=e;
    end
    if e1*e2>0
      disp('Starting points do not change sign')
      keyboard
    end
    if abs(w1-w2)<1e-7
      break;
    end
  end
  W(k)=w;
  disp(1/w);
end
subplot(2,1,1)
semilogx(nn,W);
axis('tight');
subplot(2,1,2)
semilogx(nn,1./W);
axis('tight');

function [e a]=ef(w,a)
a.z=a.f(w);
a.y([1:a.m end-a.m+1:end])=fft(a.z);
a.Z=ifft(a.y,'symmetric')*a.N/a.n;
e=max(-min(a.Z),max(a.Z-1))-1e-4;




return;





p7=[-5 0 21 0 -35 0 35 16]/32;c7=32/35;
p11=[-63 0 385 0 -990 0 1386 0 -1155 0 693 256]/512;c11=512/693;

cc=c7;
pp=p7;
d1p=c*poly_diff(pp);
d2p=c*poly_diff(d1p);
d3p=c*poly_diff(d2p);
d4p=c*poly_diff(d3p);

f1 = @(x) interp1([-1e6 -3/4 3/4 1e6],[0 0 1 1],x,'pchip'); % Gradient is 1
f2 = @(x) (1+erf(x*sqrt(pi)))/2;
f3 = @(x) (1+tanh(2*x))/2;
f4 = @(x) polyval(pp,max(-1,min(1,c*x)));
d1f = @(x) polyval(d1p,max(-1,min(1,c*x)));
d2f = @(x) polyval(d2p,max(-1,min(1,c*x)));
d3f = @(x) polyval(d3p,max(-1,min(1,c*x)));
d4f = @(x) polyval(d4p,max(-1,min(1,c*x)));

x=linspace(-2,2,1e3);
mx=filter_midpoint(x);
mmx=x(2:end-1);
mmmx=mx(2:end-1);
dx=x(2)-x(1);
subplot(4,1,1);
plot(x,f1(x),x,f2(x),x,f3(x),x,f4(x));
subplot(4,1,2);
plot(mx,diff(f1(x))/dx,mx,diff(f2(x))/dx,mx,diff(f3(x))/dx,x,d1f(x))
subplot(4,1,3);
plot(mmx,diff(f1(x),2)/dx^2,mmx,diff(f2(x),2)/dx^2,mmx,diff(f3(x),2)/dx^2,x,d2f(x));
subplot(4,1,4);
plot(mmmx,diff(f1(x),3)/dx^3,mmmx,diff(f2(x),3)/dx^3,mmmx,diff(f3(x),3)/dx^3,x,d3f(x));
axis([-2 2 -10 10]);
