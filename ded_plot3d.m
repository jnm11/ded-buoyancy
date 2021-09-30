function [a b]=ded_plot3d(nm,typ,filmp)

if nargin<2
  typ=[];
end
if nargin<3
  filmp=[];
end
if nargin<4
  fnt=[];
end

nm='008';


dd='plume/008';
nm='scalar';
typ='S';


a=ded_read_3d(dd,nm,typ);


nx=length(a.x);
nz=round(2*p.Nz);
z=filter_midpoint(linspace(0,p.L,nz+1));


figure(1);
clf;
axes;
p=ded_read_param(dd);
set(gca,'xlim',[0 p.L],'ylim',[0 p.W],'zlim',[0 p.H],'box','on','DataAspectRatio',[1 1 1]);
camlight; 
lighting('phong');
view([30 60]);
colormap(prism(28))
camup([1 0 0 ]); 
campos([25 -55 5]) 
  
xlabel('x');
ylabel('y');
zlabel('z');

mins=max(0,min(a.s(:)));
maxs=min(1,max(a.s(:)));
a.s=min(1,max(0,a.s));
sc=0.99;

sc=0.9;

for j=1:a.nt
  %  set(htt,'string',sprintf('%s, t=%7.2f',ts,a.t(j)));
  X=a.s(:,:,:,j);
  if ishandle(p)
    delete(p);
  end
  p = patch(isosurface(a.x, a.y, a.z, X));
  isonormals(a.x,a.y,a.z,X,p)
  %  p.FaceColor = 'interp';
  %p.EdgeColor = 'none';
  title(sprintf('t=%7.2f',a.t(j)));
  drawnow;
  pause;
end
