function a=ded_gc_traj(nm,xrg,force)
if nargin <2
  xrg=[];
end
if nargin <3
  force=[];
end
if iscell(nm)
  for j=1:length(nm)
    a(j)=ded_gc_traj(nm{j},xrg,force);
  end
  return;
end

if isempty(force)
  force=false;
end
d=ded_dedalus_data_dir;
dd=[d '/results/' nm];
fntraj   = [dd '/traj.mat'];
h5stats  = [d '/' nm '/stats.hdf5'];

if file_nt(h5stats,fntraj) | ~isfile(fntraj) | force
  s=ded_read_stats(nm);
  if isempty(s)
    disp(sprintf('ded_gc_traj: %s stats or param missing',nm));
    a=[];
    return;
  end
  disp(sprintf('ded_gc_traj: %s making',fntraj));

  dt=median(diff(sort(s.t)));
  [w X t]=jgrid(s.t',s.X',dt,'cubic');
  w([1 end])=[];
  X([1 end])=[];
  t([1 end])=[];
  
  f=find(w>0 & isfinite(X));
  X=interp1(t(f),X(f),t,'linear');
  
  nv=round(0.1/median(diff(unique(t))));
  
  f= find(X>=xrg(1) & X<=xrg(2));
  if isempty(f)
    a.p=[NaN NaN NaN];
    a.t1=NaN;
    a.t2=NaN;
  else
    a.p=polyfit(t(f),X(f),2);
    a.t1=t(f(1));
    a.t2=t(f(end));
  end
  a.dp=a.p(1:2).*[2 1];
  
  Xf= @(t) polyval(a.p,t);
  Uf= @(t) polyval(a.dp,t);
  
  mt = (t(1:end+1-nv)+t(nv:end))/2;
  dt = (t(1:end+1-nv)-t(nv:end));
  dX = (X(1:end+1-nv)-X(nv:end));
  V=dX./dt;
  a.t=t;
  a.X=X;
  a.Xf=Xf;
  a.mt=mt;
  a.Uf=Uf;
  a.V=V;
  a.x1=xrg(1);
  a.x2=xrg(2);
  if ~isempty(dd)
    mkdirifnotexist(fntraj);
    save(fntraj,'a');
  end
  minv=min(a.V);maxv=max(a.V);rgv=maxv-minv;minv=minv-0.05*rgv;maxv=maxv+0.05*rgv;
  mint=min(a.t);maxt=max(a.t);rgt=maxt-mint;mint=mint-0.05*rgt;maxt=maxt+0.05*rgt;
  minX=min(xrg(1),min(a.X));maxX=max(xrg(2),max(a.X));
  
  figure;clf;
  subplot(3,1,1);plot(a.mt,a.V,a.t,Uf(a.t),repmat([a.t1 a.t2],2,1),repmat([minv;maxv],1,2));title(nm);axis([mint maxt minv maxv]);
  subplot(3,1,2);plot(a.t, a.X,a.t,Xf(a.t),repmat([mint;maxt],1,2),repmat(    xrg(:)',2,1));axis([mint maxt minX maxX]);
  subplot(3,1,3);plot(Xf(a.t),Uf(a.t));axis([minX maxX -inf inf]);
  drawnow;
end
if nargout>0 & ~isempty(dd)
  disp(sprintf('ded_gc_traj: %s loading',fntraj));
  load(fntraj);
end
