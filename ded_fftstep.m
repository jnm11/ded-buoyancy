



% constrain a function between 0 and 1 for ffts of length n and n/2

N=256;
n1=N;
n2=N*3/2;




x=(0:n1-1)/n1;
x2=(0:n2-1)/n2;
X=linspace(0,1,1e4);
y1=0.25;
y2=0.75;
g=interp1([0 y1-eps y1+eps y2-eps y2+eps 1],[0 0 1 1 0 0],x,'pchip');
G=interp1([0 y1-eps y1+eps y2-eps y2+eps 1],[0 0 1 1 0 0],X,'pchip');

f=repmat(.5,1,n1);

x1=x;
s1=zeros(1,n1);
s2=zeros(1,n2);
s1(x1<y1)=-1;s1(x1>y1 & x1<y2)=1;s1(x1>y2)=-1;
s2(x2<y1)=-1;s2(x2>y1 & x2<y2)=1;s2(x2>y2)=-1;



c = @(f) ded_fftstep_con(f,n1,n2,s1,s2);

gg = @(f) sum((f-g).^2);

opt=optimset('MaxFunEval',1e5,'tolfun',1e-10,'tolx',1e-10);
f=fmincon(gg,f,[],[],[],[],zeros(1,n1),ones(1,n1),c,opt);

plot(X,G,x,f,'s-',x,ifft(fft(f),n1,'symmetric'),x2,ifft(fft(f),n2,'symmetric')*n2/n1);
