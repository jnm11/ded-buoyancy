function b=ded_pm_f7_e(nm,b,qq)

fn=ded_get_fn(nm,'final',[],'state');
if isempty(b)
  b=ded_read_hdf(fn{end});
end
c=ded_coord(nm);
m=findmax(squeeze(max(max(b.s))));
[f g w d ss rr q]=helmholtz_fft(b.u(:,:,m),b.v(:,:,m),c.dy,c.dz);
C=contourc(c.y,c.z,f,10);

figure(1);
clf;
preprint([4 4],8);
ah=jsubplot([2 2],[0.07 0.07],[0.05 0.05],[0.02 0.05]);

axes(ah(1,1));
imagesc(c.y,c.z,max(0,b.s(:,:,m)));
set(gca,'clim',[0 3.5]);
%hold('on');h=plot_contours(C);set(h,'color',[1 1 1],'linewidth',1.1);
%colorbar;
axis([-2 2 -2 2]);
title('a) Sediment','horizontalalignment','right')

axes(ah(2,1));
imagesc(c.y,c.z,max(0,b.b(:,:,m)));
axis([-2 2 -2 2]);
title('b) Saline','horizontalalignment','right')

axes(ah(1,2));
imagesc(c.y,c.z,q);
axis([-2 2 -2 2]);
title('c) Okubo-Weiss','horizontalalignment','right')

axes(ah(2,2));
h=plot_contours(C);
axis([-2 2 -2 2]);
set(h,'color',[0 0 1],'linewidth',1.1);
title('d) Streamlines','horizontalalignment','right')

set(ah(:,1),'xticklabel',[]);
set(ah(2,:),'yticklabel',[]);
set(ah,'box','on');
axes(ah(2,2));xlabel('$x$','interpreter','latex');
axes(ah(1,2));xlabel('$x$','interpreter','latex');
axes(ah(1,1));ylabel('$y$','interpreter','latex');
axes(ah(1,2));ylabel('$y$','interpreter','latex');
nm(nm=='/')='-';
print('-depsc2',sprintf('~/Dropbox/NSF-OCE-Sediments/2020NSF/ofigs/%s-xy.eps'));


figure(2);
rg=[c.y(1) c.y(end) c.x(1) c.x(end)]; 
rg=[qq.rgy qq.rgx]; 
clf;
preprint([4 4],8);
ah=jsubplot([3 1],[0.07 0.07],[0.05 0.05],[0.02 0.05]);
m=round(size(b.b,1)/2);

rgx=find(c.x>=rg(3)&c.x<=rg(4));

ss=squeeze(max(0,b.s(m,:,rgx)))';
bb=squeeze(max(0,b.b(m,:,rgx)))';
uu=squeeze(      b.u(m,:,rgx))';
% $$$ ss=ss./mean(ss,2);
% $$$ bb=ss./mean(bb,2);
% $$$ uu=uu./std(ss,1,2);
% $$$ ss=ss/max(ss(:));
% $$$ bb=bb/max(bb(:));
% $$$ uu=uu/max(uu(:));

axes(ah(1));
imagesc(c.y,c.x,ss);
set(gca,'clim',qq.sc,'ydir','normal');
axis(rg);
title('a)','horizontalalignment','right')

axes(ah(2));
imagesc(c.y,c.x,bb);
set(gca,'clim',qq.bc,'ydir','normal');
axis(rg);
title('b)','horizontalalignment','right')

axes(ah(3));
imagesc(c.y,c.x,uu);
set(gca,'clim',qq.uc,'ydir','normal');
axis(rg);
title('c)','horizontalalignment','right')


set(ah(:,1),'xticklabel',[]);
set(ah(2,:),'yticklabel',[]);
set(ah,'box','on');
axes(ah(1));xlabel('$x$','interpreter','latex');
axes(ah(2));xlabel('$x$','interpreter','latex');
axes(ah(3));xlabel('$x$','interpreter','latex');
axes(ah(1));ylabel('$z$','interpreter','latex');
set(ah,'dataaspectratio',[1 1 1]);
print('-depsc2',sprintf('~/Dropbox/NSF-OCE-Sediments/2020NSF/ofigs/%s-xz.eps'));


return;
nms={'pm/f7/e/27','pm/f7/e/28','pm/f7/e/26','pm/f7/e/25','pm/f7/e/24'};
for j=1:length(nms)
  a=ded_s_stats(nms{j},'s');
  b=ded_s_stats(nms{j},'b');
end

nm='pm/f7/e/27';
b=[];
q.sc=[0 4];
q.bc=[0 1];
q.uc=[-0.2 2];
q.rgx=[ 4 14];
q.rgy=[-2  2];

b=ded_pm_f7_e(nm,b,q)

% $$$ nms={'pm/f7/e/27','pm/f7/e/28','pm/f7/e/26','pm/f7/e/25','pm/f7/e/24'};
% $$$ for j=1:length(nms)
% $$$   p(j)=ded_read_param(nms{j});
% $$$   print
% $$$ end

%       name     Skg     Ski     Skp   num      t 
% pm/f7/e/24 0.00774 0.00774 0.00458  e/24   39.7 
% pm/f7/e/25 0.02000 0.02000 0.01000  e/25   56.7 
% pm/f7/e/26 0.04000 0.04000 0.02000  e/26   41.3 
% pm/f7/e/27 0.10000 0.10000 0.10000  e/27   44.3 



% $$$ nms={'pm/f7/e/27','pm/f7/e/28','pm/f7/e/26','pm/f7/e/25','pm/f7/e/24'};
% $$$ N=zeros(1000,length(nms));
% $$$ M={};
% $$$ for k=1:length(nms)
% $$$   S=ded_get_fn(nms{k},'s');
% $$$   h=linspace(0,10,1000);
% $$$   M{k}=zeros(MM,length(S));
% $$$   for j=1:length(S)
% $$$     disp(fns{j});
% $$$     a=ded_read_hdf(S{j});
% $$$     s=sort(a.s(:));
% $$$     n=histc(s,h);
% $$$     N(:,k)=N(:,k)+n;
% $$$   end
% $$$   plot(h(2:end),log(N(2:end)));
% $$$ end
% $$$  
% $$$ rsync -vap poseidon:pm/f7/e/27/b/b-00120.hdf5 ~/pm/f7/e/27/b
% $$$ rsync -vap poseidon:pm/f7/e/27/s/s-00120.hdf5 ~/pm/f7/e/27/s
% $$$ rsync -vap poseidon:pm/f7/e/27/b/b-00120.hdf5 ~/pm/f7/e/27/b
% $$$ rsync -vap poseidon:pm/f7/e/27/u/u-00120.hdf5 ~/pm/f7/e/27/u
% $$$ rsync -vap poseidon:pm/f7/e/27/v/v-00120.hdf5 ~/pm/f7/e/27/v
% $$$ rsync -vap poseidon:pm/f7/e/27/w/w-00120.hdf5 ~/pm/f7/e/27/w
% $$$ a=ded_read_hdf('pm/f7/e/27/s/s-00120.hdf5');
% $$$ 
% $$$ a=ded_read_hdf('pm/f7/e/27/s/s-00120.hdf5');
% $$$ 
