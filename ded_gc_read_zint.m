function a=ded_gc_read_zint(n)

fns=ded_gc_get_fn(n,'zint');

if isempty(fns)
  a=[];
  return;
end

for k=1:length(fns)
  fn=fns{k};
  a(k).x = h5read(fn,'/scales/x/1');
  a(k).y = h5read(fn,'/scales/y/1');
  a(k).t = h5read(fn,'/scales/sim_time');  
  sz=[length(a(k).y) length(a(k).x) length(a(k).t)];
  a(k).b  = reshape(h5read(fn,'/tasks/b'  ),sz);
end
a=ded_merge(a);
