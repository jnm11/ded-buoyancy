function a=ded_plume_read_param(n)
fn=sprintf('~/data/dedalus/gc/%s/param.h5',n)
%mm={'H','L','Nx','Ny','Nz','Re','Sc','T','U','U1','W','average','h','lbc','name','ubc','parallel'};
i=h5info(fn);
mm={i.Datasets.Name};
for j=1:length(mm)
  m=mm{j};
  a.(m)=h5read(fn,['/' m]);
  if iscell(a.(m)) & length(a.(m))==1
    a.(m)=a.(m){1};
  end
  if isempty(a.(m))
   end
end

return
