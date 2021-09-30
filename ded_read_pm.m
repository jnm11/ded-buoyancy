function a=ded_read_pm(n)

fns=ded_get_fn(n,'pm');
if isempty(fns)
  return;
end

for k=1:length(fns)
  fn=fns{k};
  a(k)=ded_read_hdf(fn);
end

%Rename fields by copying and then deleting
nfn={'t','x','x'};
[rfn fa fb]=intersect({'sim_time','z','y'},fieldnames(a));
nfn={nfn{fa}};
for k=1:length(fns)
  for j=1:length(rfn)
    a(k).(nfn{j})=a(k).(rfn{j});
  end
end
a=rmfield(a,rfn);

rfn=intersect({'constant','kx','ky','kz','wall_time','world_time','write_number','iteration','timestep'},fieldnames(a));
a=rmfield(a,rfn);

nm=fieldnames(a);
for k=1:length(a)
  sz=[length(a(k).x) length(a(k).t)];
  n=prod(sz);
  for j=1:length(nm)
    x=a(k).(nm{j});
    if length(x(:))==n
      a(k).(nm{j})   = reshape(x,sz);
    end
  end
end
a=ded_merge(a);
