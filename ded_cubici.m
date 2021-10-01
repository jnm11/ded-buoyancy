function f=ded_cubici(y,x1,x2)  % cubic spline interpolating between 0 and 1 between x1 and x2
if x1==x2
  f = (1-sign(y-x1))/2;
else
  x=(2*y-x1-x2)/(x2-x1);
  z=max(0,(-1+2*abs(x)));
  f=0.5+sign(x).*min(3, 2*abs(x).*(-2*x.*x+3)+z.*z.*z)/6;
end
return

x1=rand(1);
x2=x1+rand(1);
y= linspace(x1-1,x2+1,1e3);
f=ded_cubici(y,x1,x2);
plot(y,f);

L=4;
y=linspace(0,L,1e3);
x1=3;x2=3.5;x3=4;x4=4.5;
f=ded_cubici(y,x1,x2)-ded_cubici(y+L,x3,x4);
plot([y y+L],[f f]);

x1=2.5;x2=3.5;x3=3.5;x4=4.5;
f=ded_cubici(y,x1,x2)-ded_cubici(y,x3,x4)-ded_cubici(y+L,x3,x4);
plot([y y+L],[f f]);

x1=23;x2=24;x3=24;x4=25;L=24;
x1=-1;x2=0;x3=1;x4=2;L=24;
x=linspace(0,L,1e3);xx=[x-L x x+L];

f = ded_cubici(x,x1,x2)-ded_cubici(x,x3,x4)
if x4>L
  f = f + ded_cubici(x+L,x4,x3)
end
  if x1<0  
  f = f - ded_cubici(x-L,x2,x1)
end

figure(1);
subplot(2,1,1);
plot(xx,[f f f]);
subplot(2,1,2);
plot(diff([f f f]));

