function  b=ded_mavrg(nm,trg,fx,fu,typ,LL)
%ded_mavrg(nm,x,v) Moving average in frame given with position x and velocity v

p=ded_read_param(nm);
fns=ded_get_fn(nm,typ);
[a atyps tnm stnm]=ded_read_scales(fns{1});
if length(a.sz)==3
  a.sz(2)=[];
end
for k=1:length(fns)
  try
    t{k}=h5read(fns{k},stnm);
  catch
    fns={fns{1:k-1}};
    break;
  end
end
a.t=cat(1,t{:});
minx=fx(min(a.t));
maxx=fx(max(a.t));
nx=length(a.x);
dx=(a.x(end)-a.x(1))/(nx-1);

b.x=(floor(LL(1)/dx):ceil(LL(2)/dx))*dx;  % New grid 
b.z=a.z;
b.sz=[length(b.x) length(b.z)];
for k=1:length(atyps)
  b.(atyps{k})=zeros(b.sz);
end
b.T=zeros(b.sz);

for k=1:length(fns)
  disp(k)
  fn=fns{k};
  if t{k}(end)<trg(1) | t{k}(1)>=trg(2)
    continue
  end
  nt=length(t{k});
  sz=[flip(a.sz) nt];
  n=length(sz);
  sz(end+1:2)=1;
  for j=1:length(atyps)  % Need to reorder this so as to load all fields and adjust for velocity
    m=atyps{j};
    try
      x=h5read(fn,[tnm '/'  m]);
      %if prod(size(x))==prod(sz)
      c.(m) = permute(reshape(x,sz),[2 1 3])/p.W;
      %end
    catch
      disp(sprintf('ded_mavrg: Could not read %s from %s',[tnm '/'  m],fn));
    end
  end
  for i=1:length(t{k})
    u=fu(t{k}(i));
    c.u(:,:,i)  = c.u(:,:,i)-u;
    c.uu(:,:,i) = c.uu(:,:,i)-u^2;
    c.uv(:,:,i) = c.uv(:,:,i)-u*c.v(:,:,i);
    c.uw(:,:,i) = c.uw(:,:,i)-u*c.w(:,:,i);
  end
  for i=1:length(t{k})
    if t{k}<trg(1) | t{k}>=trg(2)
      coninue
    end
    X=b.x+fx(t{k}(i));
    b.T=b.T+interp1(a.x,ones(sz([2 1])),X,'spline',0);
    for j=1:length(atyps)
      m=atyps{j};
      if isfield(m,c)
        if i<=size(c.(m),3)
          b.(m)=b.(m)+interp1(a.x,c.(m)(:,:,i),X,'spline',0);
        end
      end
    end
  end
end

for j=1:length(atyps)
  m=atyps{j};
  b.(m)=(b.(m)./b.T)';
end
b.t1=trg(1);
b.t2=trg(2);

return;


nm='gc/ccle/022';
a=ded_read_stats(nm);
p=ded_read_param(nm);
k=4;
n=1;
tol1=1e-4;
tol2=0.05;
display=0;
t=a.t;
X=a.X;
v=(X(k+1:end)-X(1:end-k))./(t(k+1:end)-t(1:end-k));
tv=(t(k+1:end)+t(1:end-k))/2;
[c t1 t2 z]=fit_poly_clip(tv,v,n,tol1,tol2,display);

t1=18;
t2=42;
f=find(t>=t1 & t<t2);
px=polyfit(a.t(f),a.X(f),2);
pv=poly_diff(px,1);

fx = @(t) polyval(px,t);
fu = @(t) polyval(pv,t);
subplot(2,1,1);plot(t,X,t,fx(t));
subplot(2,1,2);plot(tv,v,t,fu(t));





clf;imagesc(b.x,b.z,b.b);set(gca,'ydir','normal');
axis([-4 1 0 1]);
hold('on');
contour(b.x,b.z,b.b,0.05:0.1:0.95);

dc
contour(b.x,b.z,b.b,0.98*[1 1]);   
axis([-6 1 0 1]);




%H=10 cm, Re = 3139
%The nose of the currents takes 23.85 sec to get to 75 cm and the ...
%total time the simulation should run is 23.85+33=56.85 sec. In that time the current propagates 178.8 cm, i.e. 17.8H in 41.2 T






