function a=ded_gc_plot_forcing(nm)
fnf=[ded_dedalus_data_dir '/' nm '/force/force_s1.hdf5'];
fnp=[ded_dedalus_data_dir '/' nm '/param.h5'];
a=ded_read_hdf(fnf);
if isempty(a)
  disp(sprintf('ded_gc_plot_forcing: No forcing data for %s',nm));
  return;
end
p=ded_read_param(nm);
n={'psi','div','fd','fb','fu','fw','wb','wu','wv','ww','wnoise'};
%n={'psi','fu','fw'};
n=intersect(n,fieldnames(a));
m=length(n);
figure(1);clf;
ah=jsubplot([1 m],[0.01 0.01],[0.02 0.04],[0.01 0.04]);
if isfield(a,'y')
  for j=1:length(n)
    a.(n{j})=squeeze(mean(a.(n{j}),2));
  end
end

for j=1:m
  if strcmp(n{j},'psi')
    axes(ah(j));contour(a.x,a.z,a.(n{j}),20);
  else
    axes(ah(j));imagesc(a.x,a.z,a.(n{j}));
  end
  y=a.(n{j});
  title(sprintf('%s min %8.3e max %8.3e', n{j},min(y(:)),max(y(:))));
end
set(ah,'ydir','normal','dataaspect',[1 1 1],'ytick',[],'xtick',[]);
set(ah,'xlim',[0 p.L],'ylim',[0 p.H]);

figure(2);clf;
ah=jsubplot([1 m],[0.05 0.05],[0.02 0.02],[0.01 0.04]);
for j=1:m
  y=squeeze(a.(n{j}));
  axes(ah(j));plot(a.x,mean(y,1),a.x,max(y,[],1),a.x,min(y,[],1));title(n{j});
  title(sprintf('%s min %8.3e max %8.3e', n{j},min(y(:)),max(y(:))));
end
set(ah,'xtick',0:p.L,'xticklabel',[],'xlim',[0 p.L],'ylim',[-inf inf],'xgrid','on','ygrid','on');
set(ah(end),'xticklabel',cellsprintf('%i',0:p.L));

figure(3);clf;
ah=jsubplot([m 1],[0.01 0.01],[0.02 0.04],[0.01 0.04]);
for j=1:m
  y=squeeze(a.(n{j}));
  axes(ah(j));plot(mean(y,2),a.z,max(y,[],2),a.z,min(y,[],2),a.z);
  title(n{j})
  %title(sprintf('%s min %8.3e max %8.3e', n{j},min(y(:)),max(y(:))));
end
set(ah,'ytick',[],'xtick',[],'ylim',[0 p.H],'xlim',[-inf inf]);


mx=filter_midpoint(a.x);
figure(4);clf;
ah=jsubplot([1 m],[0.01 0.01],[0.02 0.04],[0.01 0.04]);
for j=1:m
  X=diff(mean(squeeze(a.(n{j})),1));
  axes(ah(j));plot(mx,X);title(sprintf('d%s/dx min %8.3e max %8.3e', n{j},min(X),max(X)));
end
set(ah,'ytick',[],'xtick',[],'xlim',[0 p.L],'ylim',[-inf inf]);

if isfield(a,'fu') & isfield(a,'fw')
  figure(5);clf;
  ah=jsubplot([1 4],[0.01 0.03],[0.02 0.04],[0.01 0.04]);
  a=struct_rename_fields(a,{'fu','fw'},{'u','w'});
  a.t=[];
  %a=ded_zgrid(a,length(a.z),{'u','w'},{'u','w'},{},{},p.H);
  a.H=p.H;
  [psi phi omega div a.z]=ded_helmholtz2(a,length(a.z));
% $$$   dimx=2;
% $$$   dx=a.x(2)-a.x(1);
% $$$   div=pr_diff(a.u,dx,dimx)+a.wdz;
  axes(ah(1));
  contour(a.x,a.z,psi,20);
  title(sprintf('psi min %8.3e max %8.3e',min(psi(:)),max(psi(:))));
  axes(ah(2));
  contour(a.x,a.z,phi,20);
  title(sprintf('phi min %8.3e max %8.3e',min(phi(:)),max(phi(:))));
  axes(ah(3));
  contour(a.x,a.z,omega,20);
  title(sprintf('omega min %8.3e max %8.3e',min(omega(:)),max(omega(:))));
  axes(ah(4));
  contour(a.x,a.z,div,20);
  title(sprintf('div min %8.3e max %8.3e',min(div(:)),max(div(:))));
  set(ah(1:3),'xticklabel',{});
  set(ah,'xgrid','on','ygrid','on');
  if p.forcing==6
    xx=repmat([p.x0 p.x1 p.x2 p.x3],2,1);
  else
    xx=repmat([p.x0 p.x1 p.x2 p.x3 p.x4 p.x5 p.x6],2,1);
  end
  for j=1:length(ah)
    axes(ah(j));
    line(xx,repmat([0;p.H],1,size(xx,2)),'color',[0 0 0]);
  end
  
end



return;

a=ded_gc_plot_forcing('007');
u=squeeze(a.fu);
wu=squeeze(a.wu);
dx=a.x(3)-a.x(1);
du=(u(:,[2:end 1])-u(:,[end 1:end-1]))/dx;
d=squeeze(a.fd);
plot(a.x,du,a.x,d);
plot(a.x,wu.*(du-d)/max(wu(:)));

a=ded_gc_plot_forcing('gc/qgcf3');
nm={'fb','fd','fu'};
[b c]=ded_zgrid(a,128,nm,nm,nm,nm,2);
plot(b.fu,b.z);

nm='gc/test/11';
a=ded_read_hdf(['~/' nm '/force/force_s1.hdf5']);
p=ded_read_param(nm);
nn={'fb','fu','fw','psi'};
[b c]=ded_zgrid(a,128,nn,nn,nn,nn,p.H);
m=length(b.x);
n=length(b.z);
rgz=round(linspace(1,n,5))
subplot(4,1,1);
plot(b.x,b.fw([1 end],:));axis('tight');
xlabel('x');ylabel('w');
subplot(4,1,2);
plot(b.x,b.fu(rgz,:));axis('tight');
xlabel('x');ylabel('u');
subplot(4,1,3);
plot(b.x,[b.psi(1,:);2+b.psi(end,:)]);axis('tight');
xlabel('x');ylabel('psi');
subplot(4,1,4);
plot(b.z,b.fu(:,25));
xlabel('z');ylabel('u');



figure;
plot(b.z,b.fu(10,:));


a=ded_gc_plot_forcing(nm);

