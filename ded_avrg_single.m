function c=ded_avrg_single(nms)
if isempty(nms)
  c=[];
  return;
end

for j=1:length(nms)
  cc=ded_read_hdf(nms{j});
  if j==1
    sz=[];
    if isfield(cc,'z')
      sz(end+1)=length(cc.z);
      c.z=cc.z;
    end
    if isfield(cc,'y')
      sz(end+1)=length(cc.y);
      c.y=cc.y;
    end
    if isfield(cc,'x')
      sz(end+1)=length(cc.x);
      c.x=cc.x;
    end
    c.i=cc.i;
    c.t=cc.t;
    c.n=0;
    c.wn=[];
    c.dt=0;
    f=[];
    fnm=fieldnames(cc);
    for k=1:length(fnm)
      if prod(size(cc.(fnm{k})))==prod(sz)
        f(end+1)=k;
        c.(fnm{k})=zeros(sz);
      end
    end
    fnm={fnm{f}};
  end
  for k=1:length(fnm)
    c.(fnm{k})=c.(fnm{k})+cc.(fnm{k})*cc.dt;
  end
  c.i=min(c.i,cc.i);
  c.t=min(c.t,cc.t);
  c.dt=c.dt+cc.dt;
  c.n=c.n+1;
  c.wn(end+1)=cc.wn;
end
for k=1:length(fnm)
  c.(fnm{k})=c.(fnm{k})/c.dt;
end
c.typs=fnm;

if isfield(c,'b')
  c.Eb=ded_mix_entropy(c.b);
  c.typs={c.typs{:},'Eb'};
end
if isfield(c,'s')
  c.Es=ded_mix_entropy(c.s);
  c.typs={c.typs{:},'Es'};
end


  
        