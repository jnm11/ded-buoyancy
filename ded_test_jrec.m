function ded_test_jrec(nm,typ)
%ded_test_jrec('gc/test/01','b')

d=ded_dedalus_data_dir;
if nargin==0
  nms=cellstr_ls([d '/gc/jrec/*']);
  nms=cellstrremoveprefix(nms,[d '/']);
  for j=1:length(nms)
    ded_test_jrec(nms{j},'b');
  end
  return;
end
if nargin<2
  typ='b';
end
pr=ded_read_param(nm);
fns=cellstr_ls(sprintf('%s/%s/j%s/*',d,nm,typ));
disp(sprintf('ded_test_jrec(''%s'',''%s''), n=%u',nm,typ,length(fns)));
keyboard
for j=1:length(fns)
  a=ded_read_hdf(fns{j});
  a=ded_zgrid(a,[],{typ});
  if ndims(a.(typ))==2
    f=a.(typ);
    imagesc(a.x,a.z,f);
    set(gca,'ydir','normal');
  else
    f=permute(a.(typ),[2 3 1]);
    cla;
    p = patch(isosurface(a.x,a.y,a.z,f,max(f(:))*0.03));%,'FaceColor','red','EdgeColor','none','clipping','on');
    isonormals(a.x,a.y,a.z,f,p);
    p.Clipping = 'on';
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    p.FaceLighting = 'gouraud';
    p.AmbientStrength = 0.3;
    p.DiffuseStrength = 0.8;
    p.SpecularStrength = 0.9;
    p.SpecularExponent = 25;
    p.BackFaceLighting = 'unlit';
    axis([0 pr.L 0 pr.W 0 pr.H]);
    daspect([1 1 1]);
    lighting('gouraud');
    daspect([1 1 1]);
    set(gca,'ambientlightcolor',[1 1 1]);
    camup([0 0 1]); 
    campos(a.x(end)*[2 -1 1])
    camlight;
  end
  title(sprintf('%s, t=%7.3f, %s:[%7.3f,%7.3f],',nm,a.t,typ,min(f(:)),max(f(:))));
  drawnow;
end
