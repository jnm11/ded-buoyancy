function [c s fld]=ded_read_javrg(nm,typ,trg,sep,fns,scl)
%[a p]=ded_gc_read_javrg(nm,typ,trg) Read time averaged data for the time interval trg
% The data is normalised by the total time
% IN addition it can be scaled by scl useful for y averaged data where y-p.W
if nargin<6
  scl=[];
end
if nargin<5
  fns=[];
end
if nargin<4
  sep=[];
end
if nargin<3
  trg=[];
end
if isempty(trg)
  trg=[-inf inf];
end
if isempty(sep)
  sep=false;
end
if isempty(fns)
  fns=ded_get_fn(nm,typ);
end
if isempty(scl)
  scl=1;
end

c=struct();
s=repmat(struct(),length(fns));

fld=[];
if isempty(fns)
  c=[];
  return;
end
j=0;
m=0;
dim=[];
if trg(2)==-1
  fns={fns{end}};
  trg(2)=inf;
end
for jj=1:length(fns)
  b=ded_read_hdf(fns{jj});
  if isempty(b)
    continue;
  end
  if b.t1<trg(1) |  b.t2>trg(2)
    disp(sprintf('ded_read_javrg: read %s [%7.3f %7.3f] [%7.3f %7.3f] skipping',fns{jj},b.t1,b.t2,trg(1),trg(2)));
    continue
  end
  j=j+1;
  disp(sprintf('ded_read_javrg: read %s [%7.3f %7.3f] [%7.3f %7.3f] loading',fns{jj},b.t1,b.t2,trg(1),trg(2)));
  fld=setdiff(fieldnames(b),{'xx','yy','zz','x','y','z','dt','t1','t2','i1','i2','n','write_num','wn'});
  if j==1
    c.n=b.n;
    c.t1=double(b.t1);
    c.t2=double(b.t2);
    c.i1=b.i1;
    c.i2=b.i2;
  else
    c.n=c.n+b.n;
    c.t1=min(c.t1,b.t1);
    c.t2=max(c.t2,b.t2);
    c.i1=min(c.i1,b.i1);
    c.i2=max(c.i2,b.i2);
  end
  for k=1:length(fld)
    fnm=fld{k};
    if isfield(b,fnm)
      if ~isfield(c,fnm) 
        c.(fnm)    = squeeze(b.(fnm));
        c.dt.(fnm) = double(b.dt);
      else
        c.(fnm)    = c.(fnm)   +squeeze(b.(fnm));
        c.dt.(fnm) = c.dt.(fnm)+double(b.dt);
      end
    end
  end
  if sep
    s(j).n=b.n;
    s(j).t1=b.t1;
    s(j).t2=b.t2;
    s(j).i1=b.i1;
    s(j).i2=b.i2;
    s(j).scl=scl;
    for k=1:length(fld)
      fnm=fld{k};
      s(j).(fnm)    = scl*squeeze(b.(fnm))/double(b.dt);
      s(j).dt.(fnm) = double(b.dt);
    end
  end
end
if j==0
  s=[];
  c=[];
  return;
end
for k=1:length(fld)
  if isfield(c,fld{k}) & isfield(c.dt,fld{k})
    c.(fld{k})=scl*c.(fld{k})/c.dt.(fld{k});
  else
    c=rmfield(c,fld{k});
  end
end
if sep
  s=s(1:j);
end
c.scl=1;
return

a=ded_read_hdf('~/data/dedalus/gc/f7/g/1400/1000/a/a-00004.hdf5');
u=a.u/a.dt;
uu=a.uu/a.dt;
subplot(3,1,1);
plot(a.x,uu-u.^2)

s1=ded_read_hdf('~/data/dedalus/gc/f7/g/1400/1000/a/a-00003.hdf5');
s2=ded_read_hdf('~/data/dedalus/gc/f7/g/1400/1000/a/a-00004.hdf5');
u=(s1.u+s2.u)/(s1.dt +s2.dt);
uu=(s1.uu+s2.uu)/(s1.dt +s2.dt);
subplot(3,1,2);
plot(s1.x,uu-u.^2)

subplot(3,1,3);
[a p]=ded_gc_read_javrg('gc/f7/g/1400/1000','a',[55 176]);
plot(a.x,a.uu-a.u.^2)
keyboard