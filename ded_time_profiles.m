function ded_time_profiles(nm,X,T)
%ded_time_profile('gc/gc2d7n/73',1:10);
if nargin<2
  X=[];
end
if nargin<3
  T=[];
end

p=ded_read_param(nm);
c=ded_coord(nm);
a=ded_read_javrg(nm,'a',[0 inf],'separate');

if isempty(X)
  X=linspace(2,p.L-2);
end
if isfield(a,'t');
  t=[a.t];
else
  t=([a.t1]+[a.t2])/2;
end
if isempty(T)
  T=1:length(a);
else
  for j=1:length(T)
    T(j)=findmin(abs(t-T(j)));
  end
end
T=unique(T);

a=a(T);

x=c.Jx;
z=c.Jz;


u=cat(3,a.u);
b=cat(3,a.b);

n=round(1/c.dJx);

figure(1);clf;
f=X*n;
h=jsubplot([length(f) 1],[0.02 0.02],[0.02 0.02],[0.02 0.02]);
minu=min(u(:));
maxu=max(u(:));

for j=1:length(f);
  axes(h(j));
  plot(squeeze(u(:,f(j),:)),z);
  title(sprintf('%5.2f',x(f(j))));
  axis([minu maxu 0 p.H]);
end



figure(2);clf;
h=jsubplot([length(f) 1],[0.02 0.02],[0.02 0.02],[0.02 0.02]);
for j=1:length(f);
  axes(h(j));
  plot(squeeze(b(:,f(j),:)),z);
  title(sprintf('%5.2f',x(f(j))));
  axis([0 1 0 p.H]);
end

