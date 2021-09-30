function [a p]=ded_pm_plot_forcing(nm)
fnb=[ded_dedalus_data_dir '/' nm];
fnf=[fnb '/force/force_s1.hdf5'];
fnp=[fnb '/param.h5'];
a=ded_read_hdf(fnf);
if isempty(a)
  disp(sprintf('ded_pm_plot_forcing: No forcing data for %s',nm));
  return;
end
p=ded_read_param(fnb);
n={'fd','fb','fu','fv','fw','wb','wd','wu','wv','ww','psi','fdx','fdy','fdz'};
n=intersect(n,fieldnames(a));
m=length(n);

for j=1:m
  sz=size(a.(n{j}));
  ky=floor((1+sz(2))/2):ceil((1+sz(2))/2);
  kz=floor((1+sz(1))/2):ceil((1+sz(1))/2);
  ny{j}=['y' n{j}];a.(ny{j})=reshape(mean(a.(n{j})(kz,:,:),1),sz([2 3]));
  nz{j}=['z' n{j}];a.(nz{j})=reshape(mean(a.(n{j})(:,ky,:),2),sz([1 3]));
end

LL=[0 p.L];
HH=p.H/2*[-1 1];
WW=p.W/2*[-1 1];

dc;
figure;clf;
ah=jsubplot([2 m],[0.01 0.01],[0.02 0.04],[0.01 0.04]);
for j=1:m
  axes(ah(1,j));imagesc(a.x,a.y,a.(ny{j}));title(ny{j});xlabel('x');ylabel('y');
  axes(ah(2,j));imagesc(a.x,a.z,a.(nz{j}));title(nz{j});xlabel('x');ylabel('z');
end
set(ah,'ydir','normal','dataaspect',[1 1 1],'ytick',[],'xtick',[]);
set(ah,'xlim',LL,'ylim',WW);


figure;clf;
ah=jsubplot([1 m],[0.05 0.05],[0.02 0.02],[0.01 0.04]);
for j=1:m
  y=squeeze(a.(ny{j}));
  axes(ah(j));plot(a.x,mean(y,1),a.x,max(y,[],1),a.x,min(y,[],1));title(ny{j});
  xlabel('x');
end
set(ah,'xtick',0:p.L,'xticklabel',[],'xlim',LL,'ylim',[-inf inf],'xgrid','on','ygrid','on');
set(ah(end),'xticklabel',cellsprintf('%i',0:p.L));

figure;clf;
ah=jsubplot([1 m],[0.05 0.05],[0.02 0.02],[0.01 0.04]);
for j=1:m
  y=squeeze(a.(nz{j}));
  axes(ah(j));plot(a.x,mean(y,1),a.x,max(y,[],1),a.x,min(y,[],1));title(nz{j});
  xlabel('x');
end
set(ah,'xtick',0:p.L,'xticklabel',[],'xlim',LL,'ylim',[-inf inf],'xgrid','on','ygrid','on');
set(ah(end),'xticklabel',cellsprintf('%i',0:p.L));


figure;clf;
ah=jsubplot([m 1],[0.01 0.01],[0.02 0.04],[0.01 0.04]);
for j=1:m
  y=squeeze(a.(ny{j}));
  axes(ah(j));plot(mean(y,2),a.y,max(y,[],2),a.y,min(y,[],2),a.y);title(ny{j});
  ylabel('y');
end
set(ah,'ytick',[],'xtick',[],'ylim',WW,'xlim',[-inf inf]);

figure;clf;
ah=jsubplot([m 1],[0.01 0.01],[0.02 0.04],[0.01 0.04]);
for j=1:m
  y=squeeze(a.(nz{j}));
  axes(ah(j));plot(mean(y,2),a.z,max(y,[],2),a.z,min(y,[],2),a.z);title(nz{j});
  ylabel('z');
end
set(ah,'ytick',[],'xtick',[],'ylim',HH,'xlim',[-inf inf]);



mx=filter_midpoint(a.x);
figure;clf;
ah=jsubplot([2 m],[0.01 0.01],[0.02 0.04],[0.01 0.04]);
for j=1:m
  axes(ah(1,j));plot(mx,diff(mean(squeeze(a.(ny{j})),1)));title(ny{j});xlabel('x');
  axes(ah(2,j));plot(mx,diff(mean(squeeze(a.(nz{j})),1)));title(nz{j});xlabel('x');
end
set(ah,'ytick',[],'xtick',[],'xlim',LL,'ylim',[-inf inf]);

figure;clf;

n=intersect(n,{'fd','fu','fv','fw','fb'});
m=length(n);


k   =findmin(abs(a.x-(   0+p.x0)/2));
k(2)=findmin(abs(a.x-(p.x0+p.x1)/2));
k(3)=findmin(abs(a.x-(p.x1+p.x2)/2));
k(4)=findmin(abs(a.x-(p.x2+p.x3)/2));
ah=jsubplot([length(k),m],[0.01 0.01],[0.03 0.03],[0.03 0.03]);

for j=1:m
  x=a.(n{j});
  minx=floor(100*min(x(:)))/100;
  maxx= ceil(100*max(x(:)))/100;
  for kk=1:length(k)
    axes(ah(kk,j))
    y=x(:,:,k(kk));
    imagesc(a.y,a.z,y);
    if minx<maxx & isfinite(minx) & isfinite(maxx)
      set(ah(kk,j),'clim',[minx maxx]);
      colorbar('ytick',[minx maxx]);
    end
    title(sprintf('%s x=%4.2f',n{j},a.x(k(kk))));
  end
end
set(ah,'ydir','normal','dataaspect',[1 1 1],'ytick',[],'xtick',[]);
set(ah,'xlim',WW,'ylim',HH);


figure;
wb=squeeze(sum(sum(a.wb,1),2));
f=wb<=max(wb)/100;
de=squeeze(sum(sum(a.div.^2,1),2));
de(f)=0;
plot(a.x,squeeze(sum(sum(a.div.^2,1),2)));
title('fd-div(u)');
return;
