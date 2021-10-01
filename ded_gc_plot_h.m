function [p a b c d]=ded_gc_plot_h(nm)
if nargin<2
  quick=0;
end


nm='300';
nm=['gc/' nm];
%a=ded_read(nm,{'meant'};


p=ded_read_param(nm);
b=ded_read_mom(nm);
d=ded_read_2d(nm,'yz');
g=ded_read_2d(nm,'gc');

U=p.U;
if ~isfield(p,'B')
  p.B=1;
end

if d.t(end)>g.t(end)
  dnm=setdiff(fieldnames(d),{'x','nx','sz'});
  for j=1:length(dnm)
    d.(dnm{j})=d.(dnm{j})(:,1:end-1);
  end
end   

x=d.x;
tol=1e-3;

h1 =   mean(d.b,2)/(p.W*p.B);
h2 = 2*mean(g.bz,2)./max(tol,mean(d.b,2));
h3 =   mean(d.b,2).^2./max(tol,mean(g.bs,2));
hh=plot(x,h1,x,h2,x,h3,[0 p.L],[0.5 0.5]);
legend(hh(1:3),{'hb','hz','bb'});
axis([0 inf 0 0.6]);

a=ded_read_2d(nm,'y');
imagesc(a.x,a.z,a.S(:,:,end));set(gca,'ydir','normal');set(gca,'clim',[0 50]);

E=a.Ex+a.Ey+a.Ez;
imagesc(a.x,a.z,E(:,:,end));set(gca,'ydir','normal');set(gca,'clim',[0 10]);

E=a.Ex+a.Ey+a.Ez;
imagesc(a.x,a.z,a.S(:,:,end)-E(:,:,end));set(gca,'ydir','normal');set(gca,'clim',[-50 50]);

KE=max(0,a.uu+a.vv+a.ww-a.u.^2-a.w.^2)/2;
imagesc(a.x,a.z,KE(:,:,end));set(gca,'ydir','normal');set(gca,'clim',[0 0.01]);

