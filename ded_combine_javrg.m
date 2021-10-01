function c=ded_combine_javrg(nms)
if isempty(nms)
  c=[];
  return;
end

for j=1:length(nms)
  cc=ded_read_hdf(nms{j});
  if j==1
    c=cc;
    sz=[];
    if isfield(c,'z')
      sz(end+1)=length(c.z);
    end
    if isfield(c,'y')
      sz(end+1)=length(c.y);
    end
    if isfield(c,'x')
      sz(end+1)=length(c.x);
    end
    f=[];
    fnm=fieldnames(c);
    for k=1:length(fnm)
      if prod(size(cc.(fnm{k})))==prod(sz)
        f(end+1)=k;
      end
    end
    fnm={fnm{f}};
  else
    for k=1:length(fnm)
      c.(fnm{k})=c.(fnm{k})+cc.(fnm{k});
    end
    c.i1=min(c.i1,cc.i1);
    c.t1=min(c.t1,cc.t1);
    c.i2=max(c.i2,cc.i2);
    c.t2=max(c.t2,cc.t2);
    c.dt=c.dt+cc.dt;
    c.n=c.n+cc.n;
    c.wn(end+1)=cc.wn;
  end
end
if isfield(c,'dt')
  for k=1:length(fnm)
    c.(fnm{k})=c.(fnm{k})/c.dt;
  end
end
c.typs=fnm;

        