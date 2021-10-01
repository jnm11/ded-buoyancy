nm='gc-ccle-024.mat'

load(['~/data/dedalus/' nm]);

U=mean(a.bv(a.bt>=a.t1 & a.bt<=a.t2));
a.H=a.param.H;
a.L=a.param.L;
nz=length(a.z)*2;
%psi=ded_helmholtz2(a,nz);
b=ded_zgrid(a,nz,[],{},{},{},a.H);
rgi=find(all(isfinite(b.u),1)&all(isfinite(b.w),1));
psi=repmat(NaN,size(b.u));
psi(:,rgi)=helmholtz_fd(b.u(:,rgi),b.w(:,rgi),b.x(2)-b.x(1),b.z(2)-b.z(1),'poisson');
u=b.u/U;
axrg=[-4 1 0 1];

clf
subplot(3,1,1);
imagesc(b.x,b.z,min(1,max(0,b.b)))
colorbar;
hold('on');
C=contourc(b.x,b.z,b.b,[0.03 0.1 0.5 0.9 0.97]);
hC=plot_contours(C);set(hC,'color',[1 1 1]);
set(gca,'ydir','normal','dataaspectratio',[1 1 1],'clim',[0 1]);
axis(axrg);

subplot(3,1,2);
imagesc(b.x,b.z,u);
colorbar;
hold('on');
C=contourc(b.x,b.z,u,[-1.01 -1 -0.99 -0.01 0 0.01 0.1 0.2 0.3 0.4]);
hC=plot_contours(C);set(hC,'color',[1 1 1]);
set(gca,'ydir','normal','dataaspectratio',[1 1 1]);
axis(axrg);

subplot(3,1,3);
contour(b.x,b.z,psi,20);
set(gca,'ydir','normal','dataaspectratio',[1 1 1]);
axis(axrg);
