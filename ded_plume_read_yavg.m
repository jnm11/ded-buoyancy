function b=ded_plume_read_yavg(n)
fns=ded_get_fn(['plume/' n],'yavg');
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
  nt(k)=length(a(k).t);
  nx=length(a(k).x);
  nz=length(a(k).z);
  sz=[nz nx nt(k)];
  a(k).T = reshape(h5read(fn,'/tasks/T'  ),sz);
  a(k).S = reshape(h5read(fn,'/tasks/S'  ),sz);
  a(k).C = reshape(h5read(fn,'/tasks/C'  ),sz);
  a(k).p = reshape(h5read(fn,'/tasks/p'  ),sz);
  a(k).u = reshape(h5read(fn,'/tasks/u'  ),sz);
  a(k).w = reshape(h5read(fn,'/tasks/w'  ),sz);
end

b.T  = cat(3,a.T);
b.S  = cat(3,a.S);
b.C  = cat(3,a.C);
b.p  = cat(3,a.p);
b.u  = cat(3,a.u);
b.w  = cat(3,a.w);
b.t  = cat(1,a.t);
b.dt = cat(1,a.dt);
b.i  = cat(1,a.i);
b.x=a(1).x;
b.z=a(1).z;


%a=ded_merge(a);
