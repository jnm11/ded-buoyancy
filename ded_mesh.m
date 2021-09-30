function ded_mesh(n,m,mm)
if nargin<3
  mm=m;
end
fns=ded_get_fn(n,m);


c=ded_coord(n);
nn=length(fns);
for j=1:nn
  a=ded_read_hdf(fns{j});
  if ~isfield(a,mm) ; continue ; end
    
  f=a.(mm);
  nz=size(f,1);
  ny=size(f,2);
  if nz==c.Nz;
    z=c.z;
    y=c.y;
    x=c.x;
  elseif nz==c.NAz;
    z=c.Az;
    y=c.Ay;
    x=c.Ax;
  elseif nz==c.NJz;
    z=c.Jz;
    y=c.Jy;
    x=c.Jx;
  end
  mesh(y,z,f);
  axis('tight');
  title(sprintf('%s %s %3i %6.3f [%6.3f %6.3f]',n,m,j,a.t,min(f(:)),max(f(:))));
  drawnow;
end
