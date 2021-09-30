function a=ded_read(dd,nm,typ)

fns=ded_get_fn(dd,nm);

if isempty(fns)
  a=[];
  return;
end

if nargin<3
  typ=[];
end
if ~isempty(typ)
  if ~iscell(typ)
    typ={typ};
  end
end

[a atyps]=ded_read_scales(fns{1});
if isempty(typ)
  typ=atyps;
end

for k=1:length(fns)
  fn=fns{k};
  t{k}=h5read(fn,'/scales/sim_time');
  nt=length(t{k});
  sz=[flip(a.sz) nt];
  for j=1:length(typ)
    x=h5read(fn,['/tasks/' typ{j}]);
    if prod(size(x))==prod(sz)
      b(k).(typ{j}) = reshape(x,sz);
    end
  end
end

for j=1:length(typ)
  a.(typ{j})=cat(length(sz),b.(typ{j}));
end
a.t=cat(1,t{:});
a.nt=length(a.t);

