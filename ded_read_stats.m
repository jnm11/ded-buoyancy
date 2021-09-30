function [s fns]=ded_read_stats(nm,fix)
if nargin<2
  fix=1;
end
quiet=1;
if iscell(nm)
  s=[];
  f=[];
  k=0;
  for j=1:length(nm)
    [ss fns{j}]=ded_read_stats(nm{j},fix);
    if isempty(ss)
      continue;
    end
    k=k+1;
    f(k)=j;
    s=struct_array_append(s,ss,nm{j},quiet);
  end
  fns={nm{f}};
  return;
end

fns=ded_get_fn(nm,'stats');
if isempty(fns)
  s=[];
  return;
end
fns=fns{1};
s=ded_read_hdf(fns);
if isempty(s)
  return;
end
nm=fieldnames(s);
for j=1:length(nm)
  sz=size(s.(nm{j}));
  if sz(2)>1
    s.(nm{j})=permute(s.(nm{j}),[2 1]);
  end
end
if fix
  [t a]=unique(s.t);
  for j=1:length(nm)
    try
      s.(nm{j})=s.(nm{j})(a,:);
    end
  end
end

nm=fieldnames(s);
for j=1:length(nm)
  if isempty(s.(nm{j}))
    s=rmfield(s,nm{j});
  end
end
