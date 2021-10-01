function c=ded_gc_nose(nm,display)

%ded_gc_nose('gc/gc2d7h/20');
%s=ded_read_stats(nm);
%p=ded_read_param(nm);
%T=ded_convergence_T(nm);
%a=ded_read_g(nm,'tavg');

if nargin<2
  display=nargout==0;
end

a=ded_read_g(nm,'left');
if isempty(a)
  c=[];
  return;
end
maxb=max(a.b(:));
opt=optimset('display','none','tolx',1e-8,'tolfun',1e-8);
%figure(1);clf;
c.dx=a.x(2)-a.x(1);
a.nt=min(a.nt,size(a.b,2));
if display
  figure(1);
  clf;
end
for j=1:a.nt
  f1=max(find(a.b(:,j)>3*maxb/4));
  f2=max(find(a.b(:,j)>1*maxb/4));
  k1=(f2+f1)/2;
  k2=2+(f2-f1)/2;
  rg{j} = max(1,k1-3*k2):min(a.nx,k1+3*k2);
  y=a.x(rg{j});
  u=a.u(rg{j},j);
  my=mean(y);
  y=y-my;
  b=a.b(rg{j},j);
  %ff = @(p,z) (1-erf(p(1)*(z-p(2))))/2 - p(3)*gsl_si(p(4)*(z-p(2)));
  %p=10/(max(y)-min(y))*[1 0 0.01 1];p(3)=0.01;
  fb = @(p,z) (1-erf(p(1)*(z-p(2))))/2;
  pb=10/(max(y)-min(y))*[1 0];
  ffb=@(p) fb(p,y)-b;
  pb=lsqnonlin(ffb,pb,[0 -inf],[],opt);
  
  fu = @(p,z) p(1)-p(2)*sqrt(max(0,z-p(3))) + p(4)*max(0,p(3)-z);
  pu=[0.05 1 0 0.1];
  ffu=@(p) fu(p,y)-u;
  pu=lsqnonlin(ffu,pu,[0 0 -inf 0],[],opt);
  if display
    figure(1);
    subplot(2,1,1);
    plot(y,b,'s',y,fb(pb,y));
    hold('on');
    subplot(2,1,2);
    plot(y,u,'-s',y,fu(pu,y));
    hold('on');
  end
  c.w(j)=2*erfinv(0.9)/pb(1);
  c.X(j)=my+pb(2);
  c.e(j)=sqrt(mean(ffb(pb).^2));
  c.up=pu;
  c.ux=pu(3)+my;
  c.u0=pu(1);
  c.u1=pu(2);
  c.u2=pu(4);
end 
c.bp=abs(fft(a.b,[],1)).^2;
c.up=abs(fft(a.u,[],1)).^2;
c.bp([1 a.nx/2+1:end],:)=[];
c.up([1 a.nx/2+1:end],:)=[];
if a.nt>2
  c.bp(:,1:2)=[];
  c.up(:,1:2)=[];
end


c.t=a.t;
c.n=c.w/c.dx;
%c.kn=c.k/c.dx;

if display
  figure(2);clf;
  subplot(3,1,1);
  plot(c.t,c.n);
  ylabel('width n');
  subplot(3,1,2);
  plot(c.t,c.X);
  ylabel('X');
  subplot(3,1,3);
  plot(c.t,c.e);
  ylabel('rms error');
% $$$ subplot(4,1,3);
% $$$ plot(c.t,c.A);
% $$$ subplot(4,1,4);
% $$$ plot(c.t,c.kn);
  
% $$$ x=linspace(-20,20,1e4);
% $$$ y=gsl_si(x);
% $$$ plot(x,y);
  
  figure(3);
  mbp=mean(c.bp,2);
  mup=mean(c.up,2);
  subplot(2,1,1);loglog(mbp);axis([100 inf -inf max(mbp(100:end))]);ylabel('b power');
  subplot(2,1,2);loglog(mup);axis([100 inf -inf max(mup(100:end))]);ylabel('u power');
end

return;
clear('all');
fns=ded_dir_nms('gc/gc2d7[aeh]/*');
%fns=ded_dir_nms('gc/gc2d7a/*');
j=0;
for k=1:length(fns)
  try
    cc=ded_gc_nose(fns{k});
    pp=ded_read_param(fns{k});
  catch
    cc=[];
    pp=[];
  end
  if isempty(cc) | isempty(pp)
    disp(sprintf('failed: %s',fns{k}));
    continue
  end
  j=j+1;
  c(j)=cc;
  p(j)=pp;
end

clear X N w Re sr
for j=1:length(c)
  X(j)=mean(c(j).X);
  N(j)=mean(c(j).n);
  w(j)=mean(c(j).w);
  Re(j)=p(j).Re;
  sr{j}=p(j).series;
end

clf;
pp.lt='-';
[h lh uc]=groupplot(Re,w,sr,pp);
legend(lh,uc)
fR=find(Re<8e3);

f=@(p) max(w)*p(1).*(Re(fR)/max(Re(fR))).^(-p(2))-w(fR);
p=lsqnonlin(f,[1 1],[0 0]);
x=linspace(Re(1)*.9,Re(end)*1.1,1e4);
hold('on');
f = @(R) max(w)*p(1).*(R/max(Re(fR))).^(-p(2));

plot(x, f(x));
set(gca,'xscale','log','yscale','log');
axis('tight');
pa=ded_read_param('gc/gc2d7a/01');dxa=pa.L/pa.Nx;
ph=ded_read_param('gc/gc2d7h/01');dxh=ph.L/ph.Nx;
f(7e3)/dxa
f(2e3)*2

for j=1:20
  R=j*1e3;
  disp(sprintf('Re: %5.0f dx: %7.5f 2dx: %7.5f',R,f(R),2*f(R)));
end
