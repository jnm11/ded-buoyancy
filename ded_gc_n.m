function ded_gc_n(xt)

if nargin <1
  xt=10;
end
fns={'gc/f7/n/005','gc/f7/n/001','gc/f7/n/002','gc/f7/n/003','gc/f7/n/004'};

nxt=length(xt);
nfn=length(fns);
h1=zeros(nfn,nxt+1);
h2=zeros(nfn,nxt+1);
dc;
figure(1);clf;ah1=jsubplot([1 nfn],[0.05 0.05],[0.02 0.02],[0.02 0.02]);
for k=1:length(xt)
  figure(k+1);
  clf;
  ah2(1:2,k)=jsubplot([2 1],[0.1 0.1],[0.05 0.05],[0.05 0.05]);
end

for j=1:nfn
  nm=fns{j};
  c=ded_coord(nm);
  Nx=c.NJx;
  Nz=2*c.NJz;
  p=ded_read_param(nm);p.dim=2;
  T=ded_convergence_T(nm);
  [a fld]=ded_read_javrg(nm,'a',[T inf],'combine');
  nmd=[];nmdd=[];nmi={'u','w','b','uu','uw','ww','p'};
  b=ded_zgrid(a,Nz,[],nmd,nmdd,nmi,p.H,1);
  b.z=b.z(:);
  b.x=c.Jx(:)';
  b.dz=b.z(2)-b.z(1);
  b.dx=b.x(2)-b.x(1);

  [psi phi w d sr rr q]=ded_helmholtz3(b,p);
  axes(ah1(j));
  imagesc(c.Jx,c.Jz,b.b);
  cnts=contourc(b.x,b.z,psi,20);
  ch1=plot_contours(cnts);
  set(ch1,'color',[1 1 1]);
  set(ah1(j),'dataaspectratio',[1 1 1],'ydir','normal','clim',[0 1]);
  gu_set_aspect_ratio(ah1(j),p.H/p.L);
  axis([0 p.L 0 p.H]);
  xlabel('x');
  ylabel('z');
  [hc hp hl]=jcolorbar(ah1(j),[20 20],'b');
  f=ded_read_hdf([ded_dedalus_data_dir '/' nm '/force/fd.hdf5']);

% $$$   axes(ah2(j));
% $$$   plot(-f.fd(:,1)/min(f.fd(:,1)),c.Az);
% $$$   hold('on');
% $$$   plot(-a.u(:,1)/min(a.u(:,1)),c.Jz);
% $$$   axis([-inf inf 0 p.H]);
  for k=1:length(xt)
    ff=findmin(abs(c.Jx-xt(k)));
    axes(ah2(1,k));hold('on');h1(j,k)=plot(a.u(:,ff),c.Jz);  
    axes(ah2(2,k));hold('on');h2(j,k)=plot(a.b(:,ff),c.Jz);  
  end
  hu(j)=p.hu;
end

lstr=cellsprintf('%3.1f',hu);
for k=1:nxt
  axes(ah2(1,k));
  xlabel('$u$','interpreter','latex')
  ylabel('$z$','interpreter','latex');
  gu_setalim;
  ax=get(ah2(1,k),'xlim');
  %  h1(nxt+1,k)=line([0 0],[0 p.H],'color',[0.7 0.7 0.7 0.5]);
  set(ah2(1,k),'ylim',[0 p.H],'xlim',ax,'box','on');
  legend(h1(:,k),lstr,'location','NorthEast');
  title(sprintf('$x=%6.3f$',xt(k)),'interpreter','latex');
  
  axes(ah2(2,k));
  xlabel('$\phi$','interpreter','latex');
  %h2(nxt+1,k)=line([0 0],[0 p.H],'color',[0.7 0.7 0.7 0.5]);
  set(ah2(2,k),'xlim',[-0.05 1.05],'yticklabels',[],'box','on');
  legend(h2(:,k),lstr,'location','NorthEast');
end




