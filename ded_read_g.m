function a=ded_read_g(dd,nm,typ,trg,fnm)
%a=ded_read_g('pm/025','ayz',[],[5 inf])
if nargin < 3
  typ=[];
end
if nargin < 4
  trg=[];
end
if nargin < 5
  fnm=[];
end

if isempty(trg)
  trg=[-inf inf];
end

fns=ded_get_fn(dd,nm,[],fnm);

if isempty(fns)
  a=[];
  return;
end
if trg(2)==-1
  fns={fns{end}};
end


if nargin<3
  typ=[];
end
if ~isempty(typ)
  if ~iscell(typ)
    typ={typ};
  end
end
try
  aa=ded_read_hdf(fns{1});
catch
  a=[];
  return;
end

if isfield(aa,'ddiv');
  aa=renameStructField(aa,'ddiv','dd');
end

fnms=fieldnames(aa);
if isfield(aa,'t1') & isfield(aa,'t2')
  scl=0;
  sz=[];
  if isfield(aa,'b')
    sz=size(aa.b);
  else
    if isfield(aa,'x') & ~any(nm=='x')
      a.x=aa.x;
      sz(end+1)=length(aa.x);
    end
    if isfield(aa,'y') & ~any(nm=='y')
      a.y=aa.y;
      sz(end+1)=length(aa.y);
    end
    if isfield(aa,'z') & ~any(nm=='z')
      a.z=aa.z;
      sz(end+1)=length(aa.z);
    end
    sz(end+1:2)=1;
  end
  if isempty(typ)
    for j=1:length(fnms)
      if any(intersect(size(aa.(fnms{j})),sz))
        typ{end+1}=fnms{j};
      end
    end
    typ=setdiff(typ,{'x','y','z','i1','i2','n','dt','t1','t2','kx','ky','kz'});
  end
elseif isfield(aa,'t') & isfield(aa,'i')
  scl=0;
  typ=setdiff(fieldnames(aa),{'x','y','z','i1','i2','n','dt','t1','t2','kx','ky','kz','dt','i','t','wn'});
else
  scl=1;
  [a atyps]=ded_read_scales(fns{1});
end

t=ded_get_times(fns);

if isempty(t)
  f=length(fns);
else
  f=find( t>=trg(1) & t<=trg(2) );
end
fns={fns{f}};

for k=1:length(fns)
  try
    a=ded_read_hdf(fns{k});
  catch
    disp(sprintf('ded_read_g: failed to read file "%s"',fns{k}));
    continue;
  end
  if k==1
    aa=a;
  else
    aa=struct_array_append(aa,a);
  end
end

a=struct();
for j=1:length(typ)
  if isfield(aa,typ{j})
    n=1+max([0 find(size(aa(1).(typ{j}))>1)]);
    a.(typ{j})=cat(n,aa.(typ{j}));
  end
end
nmm=intersect({'t','t1','t2','dt','i','i1','i2'},fieldnames(aa));
for j=1:length(nmm)
  a.(nmm{j}) = cdouble(aa,nmm{j});
end

if isfield(a,'t1') & isfield(a,'t2') & ~isfield(a,'t')
  a.t=(a.t1+a.t2)/2;
end
a.typ=nm;
a.nm=dd;
return;

function x=cdouble(aa,nm)
x=zeros(1,length(aa));
for j=1:length(aa)
  x(j)=double(aa(j).(nm));
end


% $$$   if isfield(a,'dt')
% $$$     b.dt=sum(a.dt);
% $$$     for j=1:length(typ)
% $$$       b.(typ{j})=sum(a.(typ{j}),n)/b.dt;
% $$$     end
% $$$   end
if 1
else
  if isempty(typ)
    typ=setdiff(fieldnames(aa),{'Tx','Ty','Tz','constant','iteration','sim_time','timestep','wall_time','world_time','write_number','x','y','z','kx','ky','kz'});
  end
  for k=1:length(fns)
    fn=fns{k};
    t{k}=h5read(fn,'/scales/sim_time');
    nt=length(t{k});
    if nt>0
      sz=[flip(a.sz) nt];
    else
      sz=flip(a.sz);
    end
    n=length(sz);
    sz(end+1:2)=1;
    for j=1:length(typ)
      try
        ftyp=min(find(typ{j}=='_'));
        if ~isempty(ftyp)
          loc=['/tasks/' typ{j}(1:ftyp-1)];
          attr=typ{j}(ftyp+1:end);
          x=h5readatt(fn,loc,attr);
        else
          x=h5read(fn,['/tasks/' typ{j}]);
        end
        if iscell(x)
          b(k).(typ{j}) = x;
        elseif prod(size(x))==prod(sz)
          b(k).(typ{j}) = reshape(x,sz);
        end
      catch
        %        disp(sprintf('ded_read_g: failed to read /tasks/%s from %s',typ{j},fn));
      end
    end
  end
  
  typ=fieldnames(b);
  
  for j=1:length(typ)
    a.(typ{j})=cat(n,b.(typ{j}));
  end
  a.t=cat(1,t{:});
  a.nt=length(a.t);
  
  if ~isempty(trg)
    f=find(a.t>= trg(1) & a.t<=trg(2));
    a.t=a.t(f);
    a.nt=length(a.t);
    for j=1:length(typ)
      sz=size(a.(typ{j}));
      if length(sz)==2
        if sz==[a.nt 1]
          sz=a.nt;
        end
      end
      switch(length(sz))
        case 1
          a.(typ{j})=a.(typ{j})(f);
        case 2
          a.(typ{j})=a.(typ{j})(:,f);
        case 3
          a.(typ{j})=a.(typ{j})(:,:,f);
        case 4
          a.(typ{j})=a.(typ{j})(:,:,:,f);
        case 5
          a.(typ{j})=a.(typ{j})(:,:,:,:,f);
      end
    end
  end
end

