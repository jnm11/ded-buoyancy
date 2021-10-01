
N=1e3;
j=3;
m=2;
n = 1+2*j*m;
L=j;
x=linspace(0,L,n+1);x(end)=[];
X=linspace(0,L,N+1);X(end)=[];
t = sin(x*pi/2/L).^2;
T = sin(X*pi/2/L).^2;
a = j+1/2;
f = betainc(t,a,a); 
F = betainc(T,a,a); 
dF = sin(X*pi/L).^(2*j);
df = sin(x*pi/L).^(2*j);
subplot(2,1,1);plot(x,f,'s',X,F);
subplot(2,1,2);plot(X,dF);
FF=fft(dF);
ff=fft(df);
disp(abs(FF(1:2+ceil((1+n)/2))))
disp(abs(ff(1:2+ceil((1+n)/2))))



        
