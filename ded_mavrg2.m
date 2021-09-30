function  b=ded_mavrg2(nm,trg,fx,fu,typ,minx)
%ded_mavrg(nm,x,v) Moving average in frame given with position x and velocity v

p=ded_read_param(nm);
parity=ded_read_parity(nm);
cc=ded_coord(nm);
W=p.W;
fns=ded_get_fn(nm,typ);
[tt nmb]=ded_get_times(fns);
[tt f]=sort(tt);
f=f(find(tt>=trg(1) & tt<=trg(2)));
if isempty(f)
  b=[];
  return;
end

nmb={nmb{f}};
fns=unique(nmb);


a=ded_read_hdf(fns{1});
sz=[length(a.z) length(a.x)];
if isfield(a,'t1')
  mint=a.t1;
elseif isfield(a,'sim_time')
  mint=min(a.sim_time);
else
  b=[];
  disp(sprintf('ded_mavrg2: %s no time %s',nm,fns{1}));
  keyboard
  return;
end

while 1
  a=ded_read_hdf(fns{end});
  if ~isempty(a)
    break;
  end
  fns={fns{1:end-1}};
end

if isfield(a,'t1')
  maxt=a.t1;
elseif isfield(a,'sim_time')
  maxt=max(a.sim_time);
end

nx=length(a.x);
dx=(a.x(end)-a.x(1))/(nx-1);

b.x=a.x-p.L+4;
b.z=a.z;
b.sz=[length(b.x) length(b.z)];

atyps=fieldnames(a);
f=[];
a.sz=[length(a.z) length(a.x)];
for k=1:length(atyps)
  m=atyps{k};
  a.(m)=squeeze(a.(m));
  asz=size(a.(m));
  if length(asz)==length(sz)+1
    asz(end)=[];
  end
  if all(asz==sz)
    f(end+1)=k;
    %disp(m);
  else
    %disp(sprintf('%s %u %u %u %u',m,size(a.(m),1),size(a.(m),2),sz(1),sz(2)));
  end
end

atyps={atyps{f}};

b.T=0;
for j=1:length(atyps)
  b.(atyps{j})=0;
end

%disp(['atyps ' cellsprintf('%s ',atyps)]);

sz2=sz([2 1]);
for k=1:length(fns)
  fn=fns{k};
  [dd bfn]=fileparts(fn);
  try
    tt=ded_get_times(fn);
  catch
    disp(sprintf('ded_mavrg2: %s cannot read time',fn));
    continue;
  end
  if all(tt<trg(1)) 
    continue
  end
  if all(tt>=trg(2))
    continue;
  end
  try
    c=ded_read_hdf(fn);
  catch
    disp(sprintf('ded_mavrg2: %s cannot read',fn));
    continue;
  end
  for j=1:length(atyps)
    m=atyps{j};
    if isfield(c,m)
      c.(m)=squeeze(c.(m))/W;
      if any(~isfinite(c.(m)(:)))
         disp(sprintf('ded_mavrg2: %s %s not finite',fn,m));
         continue;
      end
    end
  end
  for kk=1:length(tt)
    t=tt(kk);
    if t<trg(1) | t>=trg(2)
      continue;
    end
    disp(sprintf('ded_mavrg2: %s %s, t: [%6.3f %6.3f %6.3f] X: %6.3f, V: %6.3f',nm,bfn,trg(1),t,trg(2),fx(t),fu(t)));
    u=fu(t);
    if ~isfinite(u)
      error(sprintf('ded_mavrg2: fu(t) is not finite'));
    end
    if isfield(c,'u')  ;  c.u(:,:,kk) =  c.u(:,:,kk)-u;             end
    if isfield(c,'bu') ; c.bu(:,:,kk) = c.bu(:,:,kk)-u*c.b(:,:,kk); end
    if isfield(c,'bs') ; c.bs(:,:,kk) = c.bs(:,:,kk)-u*c.s(:,:,kk); end
    if isfield(c,'uu') ; c.uu(:,:,kk) = c.uu(:,:,kk)-u^2;           end
    if isfield(c,'uv') ; c.uv(:,:,kk) = c.uv(:,:,kk)-u*c.v(:,:,kk); end
    if isfield(c,'uw') ; c.uw(:,:,kk) = c.uw(:,:,kk)-u*c.w(:,:,kk); end
    
    X=fx(t);
    fX=find(b.x+X>minx);
    for j=1:length(atyps)
      m=atyps{j};
      if isfield(c,m)
        %disp(size(b.(m)));disp(size(c.x));disp(size(c.(m)(:,:,kk)'));disp(size(X))
        b.(m)=b.(m)(:,fX)+interp1(c.x-X,c.(m)(:,:,kk)',b.x(fX),'spline',0);
      end
    end
    b.T(fX)=b.T(fX)+interp1(c.x,ones(sz2),X,'spline',0);
  end
end

b.T=max(eps,b.T');
for j=1:length(atyps)
  m=atyps{j};
  b.(m)=b.(m)'./(b.T);
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
t=a.t1;
X=a.X;
v=(X(k+1:end)-X(1:end-k))./(t(k+1:end)-t(1:end-k));
tv=(t(k+1:end)+t(1:end-k))/2;
[c t1 t2 z]=fit_poly_clip(tv,v,n,tol1,tol2,display);

t1=18;
t2=42;
f=find(t>=t1 & t<t2);
px=polyfit(a.t1(f),a.X(f),2);
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






