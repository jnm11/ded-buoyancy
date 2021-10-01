

fb = @(a,b) beta((a+1)/2,(b+1)/2)/2;
m=1e3;

x=linspace(0,pi/2,m)';
dx=x(2)-x(1);
n=40;
f1=cos(x).^(n)/fb(n,0);


f2=2*((n+1)*cos(x).^2-(n+1/2)).*cos(x).^n*sqrt((n+2)/(n+1)/beta((2*n+3)/2,1/2));
f2=cos(x).^(n+1)-fb(2*n+1,0)/f(2*n,0)*cos(x).^(n);
plot(x,f1,x,f2);
disp(sum(f1)*dx-f1(1)/2*dx);


f=cos(x).^n.*[cos(x).^0 cos(x).^1 cos(x).^2];
q=orth(f);

f1=f1/sqrt(f1'*f1);
f2=f2-f1*(f1'*f2);f2=f2/sqrt(f2'*f2);
plot(x,f1,x,f2)






d5  :=   3*(  30*cos(x)+   5*cos(3*x)- 3*cos(5*x))/416
d7  :=  15*( 105*cos(x)+  35*cos(3*x)- 7*cos(5*x)- 5*cos(7*x))/9728
d9  :=   7*( 378*cos(x)+ 168*cos(3*x)               -27*cos(7*x)- 7*cos(9*x))/20480
d11 := 315*(1386*cos(x)+ 726*cos(3*x)+99*cos(5*x)-99*cos(7*x)-55*cos(9*x)-9*cos(11*x))/4063232
u5  :=   (450*sin(x)+25*sin(3*x)-9*sin(5*x))/416
u7  :=    (11025*sin(x)+1225*sin(3*x)-147*sin(5*x)-75*sin(7*x))/9728
u9  :=   (23814*sin(x)+3528*sin(3*x)-243*sin(7*x)-49*sin(9*x))/20480
u11 :=     (4802490*sin(x)+838530*sin(3*x)+68607*sin(5*x)-49005*sin(7*x)-21175*sin(9*x)-2835*sin(11*x))/4063232



        