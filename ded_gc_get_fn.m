function fns=ded_gc_get_fn(n,m)
j=1;
k=1;
while 1
  fn1=sprintf('~/data/dedalus/gc/%s/%s/%s_s%i.h5',n,m,m,k);
  fn2=sprintf('~/data/dedalus/gc/%s/%s/%s_s%i/%s_s%i_p0.h5',n,m,m,k,m,k);
  fn3=sprintf('~/data/dedalus/gc/%s/%s/%s_s%i.hdf5',n,m,m,k);
  fn4=sprintf('~/data/dedalus/gc/%s/%s/%s_s%i/%s_s%i_p0.hdf5',n,m,m,k,m,k);
  if isfile(fn1)
    fns{j}=fn1;
    j=j+1;
  elseif isfile(fn2)
    fns{j}=fn2;
    j=j+1;
  elseif isfile(fn3)
    fns{j}=fn3;
    j=j+1;
  elseif isfile(fn4)
    j=j+1;
    fns{j}=fn4;
  end
  k=k+1;
  if k>j+100
    break;
  end
end
if j==1
  fns={};
end
