function a=ded_read_yavg(n)
fns=ded_get_fn(n,'yavg');
if isempty(fns)
  fns=ded_get_fn(n,'y');
end
if isempty(fns)
  a=[];
  return;
end

%h5disp(fn);
%a=h5info(fn);


for k=1:length(fns)
  fn=fns{k};
  a(k).x=h5read(fn,'/scales/x/1');
  a(k).z=h5read(fn,'/scales/z/1');
  a(k).dt=h5read(fn,'/scales/timestep');
  a(k).t=h5read(fn,'/scales/sim_time');
  a(k).i=h5read(fn,'/scales/iteration');
  nt=length(a(k).t);
  nx=length(a(k).x);
  nz=length(a(k).z);
  sz=[nz nx nt];
  bb=h5read(fn,'/tasks/b'  );
  a(k).b = reshape(h5read(fn,'/tasks/b'  ),sz);
  a(k).p = reshape(h5read(fn,'/tasks/p'  ),sz);
  a(k).u = reshape(h5read(fn,'/tasks/u'  ),sz);
  a(k).w = reshape(h5read(fn,'/tasks/w'  ),sz);
end

a=ded_merge(a);
