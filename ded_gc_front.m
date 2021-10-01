function xr=ded_gc_front(x,f,display)

if nargin<3
  display=[];
end
if isempty(display)
  display=[];
end

% $$$ f=b.b;
% $$$ f=a.b(1,:);
% $$$ x=c.Jx(:);
x=x(:);
f=f(:);

tol=max(f)/100;
f2=min(find(f<tol));
ff=f(f2+(-10:10));
xx=x(f2+(-10:10));
xx=xx(:);
ff=ff(:);
x0=xx(1);
xx=xx-x0;

n=3;
g   = @(q,x) exp(-polyval(q,x));
dg  = @(q,x) -polyval((n:-1:1).*q(1:n),x).*exp(-polyval(q,x));
ddg = @(q,x) polyval((n:-1:1).*q(1:n),x).^2.*exp(-polyval(q,x))-polyval((n-1:-1:1).*(n:-1:2).*q(1:n-1),x).*exp(-polyval(q,x));;
gg=@(q) g(q,xx)-ff;
q=zeros(1,n+1);
q(end-1)=5/(xx(end)-xx(1));
q(end)=-log(ff(1));
q=lsqnonlin(gg,q);
qq=q;
qq(end)=qq(end)+log(tol);
xr=roots(qq);
xr=xr(imag(xr)==0 & real(xr)>=0 & real(xr)<=xx(end));
subplot(3,1,1);plot(xx,ff,xx,g(q,xx),xr,g(q,xr),'s');
subplot(3,1,2);plot(xx, dg(q,xx),xr, dg(q,xr),'s');
subplot(3,1,3);plot(xx,ddg(q,xx),xr,ddg(q,xr),'s');

xr=xr+x0;




