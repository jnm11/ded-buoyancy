function a=ded_plot_3d(dd,nm,typ,filmp)

%a=ded_plot_3d('plume/009','scalar','T')
if nargin<2
  typ=[];
end
if nargin<3
  filmp=[];
end
if nargin<4
  fnt=[];
end

a=ded_read_3d(dd,nm,typ);

figure(1);
clf;

p=ded_read_param(dd);
lim=[0 p.H 0 p.W 0 p.L];
  
xlabel('x');
ylabel('y');
zlabel('z');


mins=max(0,min(a.T(:)));
maxs=min(1,max(a.T(:)));
a.T=min(1,max(0,a.T));
sc=0.99;

sc=0.9;

ff=linspace(0.01,0.99,8);
nf=length(ff);
p.col=jet(nf);

sz=[a.nz a.ny a.nz];
xx=repmat(reshape(a.x,[1 1 a.nx]),[a.nz a.ny    1]);
yy=repmat(reshape(a.y,[1 a.ny 1]),[a.nz 1    a.nx]);
zz=repmat(reshape(a.z,[a.nz 1 1]),[1    a.ny a.nx]);


for j=1:a.nt
  cla;
  %  set(htt,'string',sprintf('%s, t=%7.2f',ts,a.t(j)));
  X=permute(a.T(:,:,:,j),[2 3 1]);
  if ishandle(p)
    delete(p);
  end
  for i=1:nf
    hp(i) = patch(isosurface(a.x, a.y, a.z, X,ff(i)),'FaceColor',p.col(i,:),'EdgeColor','none','clipping','on');
    isonormals(a.x,a.y,a.z,X,hp(i));
  end
  view([40 20]);
  camlight(80,40);
  camlight('left');
  lighting('phong');
  title(sprintf('t=%7.2f',a.t(j)));
  daspect([1 1 1]);
  axis(lim);
  xlabel('x');ylabel('y');zlabel('z');
  drawnow;
  pause(0.01);

% $$$ camlight; 
% $$$ lighting('phong');
% $$$ view([70 -30]);
% $$$ colormap(prism(28))
% $$$ camup([1 0 0 ]); 
% $$$ campos([25 -55 5]) 

end

