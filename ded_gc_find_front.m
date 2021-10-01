function b=ded_gc_find_front(n,T)

if nargin<2
  T=[];
end
b.X=[];
b.t=[];
b.Xm=NaN;
if isempty(n)
  return;
end

if isstruct(n)
  a=n;
else
  a=ded_read_flux(nm);
end
if isempty(a)
  return;
end

if isnumeric(T)
  b.T=T;
else
  b.T=ded_convergence_T(a.nm);
end

if isfinite(b.T)
  T=b.T;
elseif isfield(a,'t')
  T=a.t(end);
else
  T=NaN;
end
       
nt=length(a.t);
for j=1:nt
  b.X(j)=find_X(a.b(:,j),a.x);
end
b.t=a.t;
b.T=T;
f=find(a.t>=T);
b.Xm=find_X(mean(a.b(:,f),2),a.x);
return

function X=find_X(b,x)
if 0
  n=round(50+50*rand(1));
  x=linspace(0,1,n+1);x(end)=[];
  b=sin(2*pi*10*x)+cos(2*pi*11*x);
  b=randn(1,n);
  xx=x;bb=b;
  
  f=fft(b);
  m=floor((n-1)/2);

  ff=zeros(1,2*n+1);
  ff(1:m-1)=f(1:m-1);
  ff(end-n+m:end)=f(m:n);
  b=2*ifft(ff,[],'symmetric');
  dx=(xx(2)-xx(1))/2;
  x=(0:2*n)*dx;
  subplot(2,1,1)
  plot(xx,bb,'-s',x,b,'^-');
  subplot(2,1,2)
  plot(xx,bb-b(1:2:end));
end
% $$$ dx=x(2)-x(1);
% $$$ xx=(0:2*length(x))*dx/2;
% $$$ b=interp1([x;x(end)+dx],[b;b(1)],xx,'pchip');
% $$$ x=xx;

b=b-max(b,[],1)/20;
f=max(find(b(1:end-1).*b(2:end)<0));
if isempty(f)
  X=NaN;
else
  b1=b(f);
  b2=b(f+1);
  X = (b1*x(f)-b2*x(f-1))/(b1-b2);
end
