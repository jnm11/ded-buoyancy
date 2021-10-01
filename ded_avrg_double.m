function c=ded_avrg_double(nms1,nms2)
if isempty(nms1) | isempty(nms2)
  c=[];
  return;
end
n=min(length(nms1),length(nms2));
fnm={};
for j=1:n
  c1=ded_read_hdf(nms1{j});
  c2=ded_read_hdf(nms2{j});
  if c1.dt~=c2.dt
    error('ded_avrg_double: dt is not equal');
  end
  if c1.t~=c2.t
    error('ded_avrg_double: t is not equal');
  end
  
  if j==1
    sz=[];
    if isfield(c1,'z')
      sz(end+1)=length(c1.z);
      c.z=c1.z;
    end
    if isfield(c1,'y')
      sz(end+1)=length(c1.y);
      c.y=c1.y;
    end
    if isfield(c1,'x')
      sz(end+1)=length(c1.x);
      c.x=c1.x;
    end
    c.i=c1.i;
    c.t=c1.t;
    c.n=0;
    c.wn=[];
    c.dt=0;
    f=[];
    f1=fieldnames(c1);
    for k=1:length(f1)
      if prod(size(c1.(f1{k})))==prod(sz)
        f(end+1)=k;
      end
    end
    f1={f1{f}};
    f=[];
    f2=fieldnames(c2);
    for k=1:length(f2)
      if prod(size(c2.(f2{k})))==prod(sz)
        f(end+1)=k;
      end
    end
    f2={f2{f}};
  end
  for k1=1:length(f1)
    for k2=1:length(f2)
      nm=sort([f1{k1} f2{k2}]);
      if j==1
        c.(nm)=0;
        fnm{end+1}=nm;
      end
      c.(nm)=c.(nm)+c1.(f1{k1}).*c2.(f2{k2})*c1.dt;
    end
  end
  c.i=min(c.i,c1.i);
  c.t=min(c.t,c1.t);
  c.dt=c.dt+c1.dt;
  c.n=c.n+1;
  c.wn(end+1)=c1.wn;
end
for k=1:length(fnm)
  c.(fnm{k})=c.(fnm{k})/c.dt;
end
c.typs=fnm;

        