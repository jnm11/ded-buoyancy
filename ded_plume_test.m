function a=ded_plume_test(n)

j=1;
for
  fn=sprintf('%s/pm/test/yavg/yavg_s%u.hdf5',ded_dedalus_data_dir,j);
  if ~isfile(fn);
    break;
  end
  j=j+1;
  a=ded_read_hdf(fn);
end
imagesc(a)



  fn=fns{end};
else
  fn=sprintf('%s/pm/plume_s%i/plume_s%i_p0.h5',ded_dedalus_data_dirn,n);
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
a.c=h5read(fn,'/tasks/c');
a.s=h5read(fn,'/tasks/s');
a.T=h5read(fn,'/tasks/T');
a.u=h5read(fn,'/tasks/u');
a.w=h5read(fn,'/tasks/w');
try
  a.fT=h5read(fn,'/tasks/fT');
  a.fs=h5read(fn,'/tasks/fs');
  a.fc=h5read(fn,'/tasks/fc');
  a.wg=h5read(fn,'/tasks/wg');
  a.wr=h5read(fn,'/tasks/wr');
end

%a.fT=h5read(fn,'/tasks/fT');
%a.wg=h5read(fn,'/tasks/wg');
%a.IS1=h5read(fn,'/tasks/IS1');
%a.IS2=h5read(fn,'/tasks/IS2');
%a.fT=a.fT(:,:,1);
%a.wg=a.wg(:,:,1);
%a.IS1=a.IS1(:,:,1);
%a.IS2=a.IS2(:,:,1);

% $$$ a.Q =h5read(fn,'/tasks/Q');

L=round(a.z(end));
H=round(a.x(end));
nx=length(a.x);
nz=round(L/H*nx);
z=linspace(a.z(1),a.z(end),nz);

ylim=[1 9]*H/10;
xlim=[0 L];
f2=find(a.x>=ylim(1) & a.x<=ylim(2));
f1=find(a.z>=xlim(1) & a.z<=xlim(2));

minc=min(min(min(a.c(f1,f2,:))));maxc=max(max(max(max(a.c(f1,f2,:)))),minc+1e-3);
mins=min(min(min(a.s(f1,f2,:))));maxs=max(max(max(max(a.s(f1,f2,:)))),mins+1e-3);
minT=min(min(min(a.T(f1,f2,:))));maxT=max(max(max(max(a.T(f1,f2,:)))),minT+1e-3);
minu=min(min(min(a.u(f1,f2,:))));maxu=max(max(max(max(a.u(f1,f2,:)))),minu+1e-3);
minw=min(min(min(a.w(f1,f2,:))));maxw=max(max(max(max(a.w(f1,f2,:)))),minw+1e-3);

figure(1);
clf;

if L>2*H
  asz=[1 5];
else
  asz=[2 3];
end

aa=jsubplot(asz,[0.05 0.03],[0.02 0.02],[0.01 0.01]);

axes(aa(1));hc=imagesc(z,a.x,NaN*a.c(:,:,1),[minc maxc]);c1=colorbar;%('position');%,[0.93 0.81 0.02 0.18]);
axes(aa(2));hs=imagesc(z,a.x,NaN*a.s(:,:,1),[mins maxs]);c1=colorbar;%('position');%,[0.93 0.61 0.02 0.18]);
axes(aa(3));hT=imagesc(z,a.x,NaN*a.T(:,:,1),[minT maxT]);c2=colorbar;%('position');%,[0.93 0.41 0.02 0.18]);
axes(aa(4));hu=imagesc(z,a.x,NaN*a.u(:,:,1),[minu maxu]);c3=colorbar;%('position');%,[0.93 0.21 0.02 0.18]);
axes(aa(5));hw=imagesc(z,a.x,NaN*a.w(:,:,1),[minw maxw]);c3=colorbar;%('position');%,[0.93 0.01 0.02 0.18]);
set(aa(:,1:end-1),'xticklabel',[]);
set(aa(2:end,:),'yticklabel',[]);
set(aa,'xtick',[0 L]);
set(aa,'ytick',[0 H]);
set(aa,'xlim',[0 L],'ylim',ylim,'DataAspectRatio',[1 1 1],'box','on','ydir','normal');
set(aa(3),'FontSize',8);
axes(aa(1));htc=text(L/200,H,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
axes(aa(2));hts=text(L/200,H,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
axes(aa(3));htT=text(L/200,H,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
axes(aa(4));htu=text(L/200,H,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
axes(aa(5));htw=text(L/200,H,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);


set([htc hts htT htu htw],'BackgroundColor',[0 0 0 .3]);
if length(aa(:))==6
  delete(aa(6));
end

S= 600/sqrt(prod(asz)*H*L);
FX=S*asz(1)*L;
FY=S*asz(2)*H ;
set(gcf,'units','pixels');
pf=get(gcf,'position');
set(gcf,'position',[pf(1:2),20+round(FX),20+round(FY)]);

for j=1:length(a.t);
  c=interp1(a.z,a.c(:,:,j),z,'linear');
  s=interp1(a.z,a.s(:,:,j),z,'linear');
  T=interp1(a.z,a.T(:,:,j),z,'linear');
  u=interp1(a.z,a.u(:,:,j),z,'linear');
  w=interp1(a.z,a.w(:,:,j),z,'linear');
  set(htc,'string',sprintf('t=%5.4f, i=%4u, c=%6.3f',a.t(j),a.i(j),mean(c(:))));
  set(hts,'string',sprintf('s=%6.3f',mean(s(:))));
  set(htT,'string',sprintf('T=%6.3f',mean(T(:))));
  set(htu,'string',sprintf('u=%6.3f',mean(u(:))));
  set(htw,'string',sprintf('w=%6.3f',mean(w(:))));
  set(hc,'CData',c');
  set(hs,'CData',s');
  set(hT,'CData',T');
  set(hu,'CData',u');
  set(hw,'CData',w');
  drawnow;
  pause(0.1);
end

a.aa=aa;

return;

figure(3);
nt=size(a.s,3);
mw=round(max(a.wg(:)));
mw=100;
for j=1:nt
  subplot(2,1,1);plot(a.x,mean(a.s(:,:,j),1),'b',a.x,mean(a.fs(:,:,j),1),'b--',a.x,mean(a.wg(:,:,j),1)/mw,'r');
  subplot(2,1,2);plot(a.x,mean(a.T(:,:,j),1),'b',a.x,mean(a.fT(:,:,j),1),'b--',a.x,mean(a.wg(:,:,j),1)/mw,'r');
  pause(0.1);
end
return;
dc;
aa=jsubplot([1 5],[0.05 0.03],[0.02 0.02],[0.01 0.01]);
axes(aa(1));imagesc(a.z,a.x,mean(a.fT,3)');set(gca,'ydir','normal','dataaspect',[1 1 1]);colorbar;
axes(aa(2));imagesc(a.z,a.x,mean(a.fs,3)');set(gca,'ydir','normal','dataaspect',[1 1 1]);colorbar;
axes(aa(3));imagesc(a.z,a.x,mean(a.fc,3)');set(gca,'ydir','normal','dataaspect',[1 1 1]);colorbar;
axes(aa(4));imagesc(a.z,a.x,mean(a.wg,3)');set(gca,'ydir','normal','dataaspect',[1 1 1]);colorbar;
axes(aa(5));imagesc(a.z,a.x,mean(a.wr,3)');set(gca,'ydir','normal','dataaspect',[1 1 1]);colorbar;
