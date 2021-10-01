function ded_gc_front_b(nm)
%ded_gc_front_b('gc/f7/l/05000/01');
%ded_gc_front_b('gc/f7/l/05000/02');
%nm='gc/f7/l/05000/02'
fns=ded_get_fn(nm,'b');
c=ded_coord(nm);
nt=length(fns);
nz=c.NJz;
nx=c.NJx;
bl=zeros(nz,nt);
br=zeros(nz,nt);
bt=zeros(nx,nt);
bb=zeros(nx,nt);
t=zeros(1,nt);
for j=1:nt
  a=ded_read_hdf(fns{j});
  dim=ndims(a.b);
  if dim==3
    a.b=squeeze(mean(a.b,2));
  end
  bl(:,j)=a.b(:,1);
  br(:,j)=a.b(:,end);
  bb(:,j)=a.b(1,:);
  bt(:,j)=a.b(end,:);
  t(j)=a.t;
end

p=ded_read_param(nm);
w=ichebintw(c.NJz);
Ib=w*a.b;
b1=max(a.b,[],1);
b1=2*b1/max(b1)-1;

f=max(find(b1(1:end-1).*b1(2:end)<0));
x1=c.Jx(f);x2=c.Jx(f+1);
a1=b1(f);a2=b1(f+1);
X=(x1*a2-x2*a1)/(a2-a1);
dX=(x1-x2)/(a2-a1);
X1=X-3*dX;
X2=X+3*dX;
f1=find(c.Jx>=X1&c.Jx<=X2);
f2=min(find(c.Jx>X2))+(0:10);
f3=max(find(c.Jx<X1))-(0:10);
f4=find(c.Jx>(p.L-0.5+X)/2);
f5=1:20;2
clf;

% $$$ fp=get(gcf,'position');
% $$$ fp(3:4)=[768 768];
% $$$ set(gcf,'position',fp);

xf=c.Jx(f1);
bf=a.b(1,f1);
fb=@(p,xx) p(3)*( 1-erf((xx-p(2))/p(1)*100) )/2 + exp(-((xx-p(2))/p(1)*100).^2 ).*(p(4)+p(5)*((xx-p(2))/p(1)*100));
pp = [5 X 1 0 0];
Pe=p.Re*p.Scb; 

% $$$ fb=@(q,xx) q(4)-ded_b_front_fun(max(0,xx-q(2)),q(5),q(1)*Pe,p.U,q(3));
% $$$ pp = [1 X 1 1 1 1];
% $$$ fb(pp,xf);

ff = @(p) fb(p,xf(:))-bf(:);

pp=lsqnonlin(ff,pp);disp(pp);
xx=linspace(0,p.L,1e5);
xx1=xx(xx>=X1&xx<=X2);
%plot(xf,bf,'s',xx1,fb(pp,xx1));

figure(1);clf;
ah=jsubplot([1 4],[0.1 0.02],[0.02 0.05],[0.01 0.1]);
axes(ah(1));
plot(c.Jx(f1),a.b(1,f1),'-s',xx1,fb(pp,xx1),'-');axis('tight');
title(sprintf('%s front and fit',nm));
axes(ah(2));
plot(c.Jx(f1),a.b(1,f1)-fb(pp,c.Jx(f1))','-s');axis('tight');
title('front minus fit');
axes(ah(3));
plot(c.Jx(f4),a.b(1,f4),'s-');axis('tight');
title('right forcing');
axes(ah(4));
plot(c.Jx(f5),1-a.b(1,f5),'s-');axis('tight');
set(ah,'xtick',[]);
title('left forcing');

% int |db/dx| dx as  afuntion of time
pbt=sum(abs(diff(bt)));
pbb=sum(abs(diff(bb)));
pbl=sum(abs(diff(bl)));
pbr=sum(abs(diff(br)));

figure(2);clf;
subplot(4,1,1);plot(t,pbr);ylabel('right');
title('Total b deviation');
subplot(4,1,2);plot(t,pbl-1);ylabel('left');
subplot(4,1,3);plot(t,pbt);ylabel('top');
subplot(4,1,4);plot(t,pbb-1);ylabel('botton');
xlabel('t');

