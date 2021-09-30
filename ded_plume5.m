function ded_plume5(n)
if nargin<1
  fns=sort(cellstr_ls('~/data/dedalus/plume/*/*.h5'));
  fn=fns{end};
else
  fn=sprintf('~/data/dedalus//plume/plume_s%i/plume_s%i_p0.h5',n,n);
end

if ~isfile(fn)
  disp(fn);
  return;
end
disp(fn);
%h5disp(fn);
a.x=h5read(fn,'/scales/x/1');
a.z=h5read(fn,'/scales/z/1');
a.dt=h5read(fn,'/scales/timestep');
a.t=h5read(fn,'/scales/sim_time');
a.i=h5read(fn,'/scales/iteration');
a.b=h5read(fn,'/tasks/b');
a.T=h5read(fn,'/tasks/T');
a.u=h5read(fn,'/tasks/u');
a.w=h5read(fn,'/tasks/w');
L=round(a.x(end));
H=round(a.z(end));
nx=length(a.x);
nz=round(H/L*nx);
z=linspace(a.z(1),a.z(end),nz);

maxu=max(a.u(:));minu=min(a.u(:));
maxw=max(a.w(:));minw=min(a.w(:));
maxb=max(1,max(a.b(:)));minb=min(0,min(a.b(:)));
maxT=max(1,max(a.T(:)));minT=min(0,min(a.T(:)));

if minu==maxu
  minu=-1;
  maxu=1;
end
if minw==maxw
  minw=-1;
  maxw=1;
end

figure(1);
clf;
a1=axes('position',[0.01 0.76 0.9 0.23]);
a2=axes('position',[0.01 0.51 0.9 0.23]);
a3=axes('position',[0.01 0.26 0.9 0.23]); 
a4=axes('position',[0.01 0.01 0.9 0.23]); 
axes(a1);hb=imagesc(a.x,z,NaN*a.b(:,:,1),[minb maxb]);c1=colorbar('position',[0.93 0.76 0.02 0.23]);
axes(a2);hT=imagesc(a.x,z,NaN*a.T(:,:,1),[minT maxT]);c1=colorbar('position',[0.93 0.51 0.02 0.23]);
axes(a3);hu=imagesc(a.x,z,NaN*a.u(:,:,1),[minu maxu]);c2=colorbar('position',[0.93 0.26 0.02 0.23]);
axes(a4);hw=imagesc(a.x,z,NaN*a.w(:,:,1),[minw maxw]);c3=colorbar('position',[0.93 0.01 0.02 0.23]);
set([a1 a2 a3],'xticklabel',[]);
set([a1 a2 a3 a4],'xlim',[0 L],'ylim',[0 H],'DataAspectRatio',[1 1 1],'box','on','ydir','normal');
set(a3,'FontSize',8);
axes(a1);ht=text(L/200,H,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);


S=512/L;
set(gcf,'units','pixels');
pf=get(gcf,'position');
set(gcf,'position',[pf(1:2),20+round(S*L),20+4*round(S*H)]);

for j=1:length(a.t);
  b=interp1(a.z,a.b(:,:,j),z,'linear');
  T=interp1(a.z,a.T(:,:,j),z,'linear');
  u=interp1(a.z,a.u(:,:,j),z,'linear');
  w=interp1(a.z,a.w(:,:,j),z,'linear');
  set(ht,'string',sprintf('t=%5.4f, i=%4u, u=%7.5f, w=%7.5f',a.t(j),a.i(j),mean(u(:)),mean(w(:))));
  set(hb,'CData',b');
  set(hT,'CData',T');
  set(hu,'CData',u');
  set(hw,'CData',w');
  drawnow;
  pause(0.1);
end
