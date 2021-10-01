



% constrain a function between 0 and 1 for ffts of length n and n/2

N=256;
n1=N;
n2=N*3/2;

% Try to fit H(z-x)
% equality constraint of zero must be satsified for x>z+dz

z=0.1;
dz=0.2;
x1=(0:n1-1)/n1;
x2=(0:n2-1)/n2;
X=linspace(0,1,1e4);


g1=double(abs(x1-0.5)<z); % This is the step function we wish to fit
g2=double(abs(x1-0.5)<z); % This is the step function we wish to fit

g1=exp(-20*(x1-0.5).^2/z);
f1=g1(x<=0.5-z-dz);
f2=g1(x>=0.5+z+dz);
f0=g1(abs(x1-0.5)<z+dz);

ga=f0;

gg = @(f) sum((f-ga).^2);

c = @(f) ded_fftstep_con2([f1 f f2],n1,n2);

opt=optimset('MaxFunEval',1e5,'tolfun',1e-10,'tolx',1e-10);
f0=fmincon(gg,f0,[],[],[],[],zeros(size(f0)),[],c,opt);
f=[f1 f0 f2];

plot(X,abs(X-0.5)<z,x1,f,'s-',x1,ifft(fft(f),n1,'symmetric'),x2,ifft(fft(f),n2,'symmetric')*n2/n1);
