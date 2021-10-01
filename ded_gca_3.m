function [a p]=ded_gca_3(a,p,nm)
% ded_gca_3 analysis of gravity current using javrg
%[a p]=ded_gca_3('gc/f6/g/26');
%[a p]=ded_gca_3('gc/f6/k/01');
%nm='gc/test/50';

if nargin==1
  nm=a;
  [a p]=ded_gc_read_javrg(nm,[0 inf]);
end
s=ded_read_stats(nm);
X=round(s.X(end));
T=a.dt;
t1=a.t1;
t2=a.t2;

dim =2 + isfield(a,'vv');

trg=[a.t1 a.t2];
PS=p.B*p.g*p.hb;
nu= 1/p.Re;
xdmax=max([p.x0 p.x1 p.x2 p.x3 p.x4 p.x5 p.x6 p.x7]);
xrg=[0 p.L];
prg=[0 p.L -0.5 0.5];
fx=find(a.x>=xdmax & a.x<=X+1);

x=a.x(fx);
xrg=[x(1) x(end)];
dimz=1;
dimx=2;
dx=a.x(2)-a.x(1);

fnm={'u','w','b','p','S','Wx','Wy','Wz','bu','bw','bb','uu','vv','ww','uw','dudz','dwdz','dpdz','dbdz','udwdx','wdudz','dudzz','duwdz','P','dudx','dwdx','dpdx','dPdx','dwdxx','duwdx','duudx','dvvdx','dwwdx','dudxx','dwdxx'};
fnd={'u','w','b','p','bu','uu','vv','ww','uw','P'};
fndd={'u','w','b','p'};
fnI={'u','w','b','uu','vv','ww','uw','p','duwdx','udwdx','wdudz','dwdxx'};

a.P     = a.p+a.uu/2+a.ww/2;

a.dwwdz = ichebdifff(a.ww,dimz,p.H,1);
a.duwdz = ichebdifff(a.uw,dimz,p.H,1);

a.dudx  = pr_diff( a.u,dx,dimx,1);
a.dwdx  = pr_diff( a.w,dx,dimx,1);
a.dpdx  = pr_diff( a.p,dx,dimx,1);
a.dPdx  = pr_diff( a.P,dx,dimx,1);
a.dwdxx = pr_diff( a.w,dx,dimx,2);
a.duwdx = pr_diff(a.uw,dx,dimx,1);
a.duudx = pr_diff(a.uu,dx,dimx,1);
a.dvvdx = pr_diff(a.vv,dx,dimx,1);
a.dwwdx = pr_diff(a.ww,dx,dimx,1);

a.dudxx = pr_diff( a.u,dx,dimx,2);
a.dwdxx = pr_diff( a.w,dx,dimx,2);

a.wdudz = a.duwdz+a.duudx/2;
a.udwdx = a.duwdx+a.dwwdz/2;

a=ded_zgrid(a,length(a.z)*2,fnm,fnd,fndd,fnI,p.H,1);

a.p=a.p-a.p(end,end);
a.P=a.P-a.P(end,end);

z=a.z;



a.wdudzIx = pr_int(a.wdudz,dx,dimx);
a.dudzzIx = pr_int(a.dudzz,dx,dimx);
a.dudzzIx = a.dudzzIx-repmat(a.dudzzIx(:,end),1,size(a.dudzzIx,2));

a.dpgdz=a.dpdz+p.g*a.b;
a.pg  = a.p  +p.g*a.bIz;
Mxa =  a.duudx  +a.duwdz+a.dpdx -nu*(a.dudxx+a.dudzz);
Mxb =  a.duudx/2+a.wdudz+a.dpdx -nu*(a.dudxx+a.dudzz);
Mxc = -a.dwwdx/2+a.wdudz+a.dPdx -nu*(a.dudxx+a.dudzz);
Mza =  a.dwwdz  +a.duwdx+a.dpgdz-nu*(a.dwdxx+a.dwdzz);
Mzb =  a.dwwdz/2+a.udwdx+a.dpgdz-nu*(a.dwdxx+a.dwdzz);
Mzc = -a.duudz/2+a.udwdx+a.dPdz -nu*(a.dwdxx+a.dwdzz)+p.g*a.b;
max(max(abs(Mxa-Mxb)))
max(max(abs(Mxa-Mxc)))
max(max(abs(Mza-Mzb)))
max(max(abs(Mza-Mzc)))

dc;
figure;clf
subplot(4,1,1);
h=plot(x,Mxa([1 end],fx));
axis('tight');
legend(h,{'Bottom','Top'},'location','best');
title(sprintf('nm: %s, $dt$: %6.2f, $M_x$',nm,T),'interpreter','latex');
ylabel('$M_x$','interpreter','latex');
subplot(4,1,2);h=plot(x,    a.dPdx([1 end],fx));axis('tight');ylabel('$P_x$','interpreter','latex');
subplot(4,1,3);h=plot(x,nu*a.dudxx([1 end],fx));axis('tight');ylabel('$\nu\, u_{xx}$','interpreter','latex');
subplot(4,1,4);h=plot(x,nu*a.dudzz([1 end],fx));axis('tight');ylabel('$\nu\, u_{zz}$','interpreter','latex');
xlabel('$x$','interpreter','latex');


Ix = -a.ww/2+a.wdudzIx+a.P-nu*(a.dudx+a.dudzzIx);Ix=Ix-repmat(Ix(:,end),1,size(Ix,2));
Bx = a.P(  1,:)-nu*(a.dudx(  1,:)+a.dudzzIx(  1,:));Bx=Bx-Bx(end);
Tx = a.P(end,:)-nu*(a.dudx(end,:)+a.dudzzIx(end,:));Tx=Tx-Tx(end);
figure;clf
subplot(4,1,1);
h=plot(x,Bx(:,fx),x,Tx(:,fx));
axis('tight');
legend(h,{'Bottom','Top'},'location','best');
title(sprintf('nm: %s, dt: %6.2f, $\\int M_x\\,dx$',nm,T),'interpreter','latex');
ylabel('$M_x$','interpreter','latex');
subplot(4,1,2);h=plot(x,   a.P(     [1 end],fx));axis('tight');ylabel('$P$','interpreter','latex');
subplot(4,1,3);h=plot(x,nu*a.dudx(   [1 end],fx));axis('tight');ylabel('$\nu\, u_x$','interpreter','latex');
subplot(4,1,4);h=plot(x,nu*a.dudzzIx([1 end],fx));axis('tight');ylabel('$\nu\,\int u_{zz}\,dx$','interpreter','latex');
xlabel('$x$','interpreter','latex');

ff1=find(a.x>=4 & a.x<=16); % Average region
ff2=length(a.x);

Mzb =  a.dwwdz/2+a.udwdx+a.dpgdz-nu*(a.dwdxx+a.dwdzz);
Izb =  a.ww/2+a.udwdxIz+a.pg-nu*(a.dwdxxIz+a.dwdz);
T1=[mean(        Mzb(:,ff1),2),mean(       Mzb(:,ff2),2)];
T2=[mean(0.5*a.dwwdz(:,ff1),2),mean(0.5*a.dwwdz(:,ff2),2)];
T3=[mean(    a.udwdx(:,ff1),2),mean(   a.udwdx(:,ff2),2)];
T4=[mean( nu*a.dwdxx(:,ff1),2),mean( nu*a.dwdxx(:,ff2),2)];
T5=[mean( nu*a.dwdzz(:,ff1),2),mean( nu*a.dwdzz(:,ff2),2)];

figure;
clf;
subplot(5,1,1);h=plot(z,T1);axis('tight');ylabel('$M_z$','interpreter','latex');
subplot(5,1,2);plot(z,T2);axis('tight');ylabel('$ww_z$','interpreter','latex');
subplot(5,1,3);plot(z,T3);axis('tight');ylabel('$uw_x$','interpreter','latex');
subplot(5,1,4);plot(z,T4);axis('tight');ylabel('$\nu\,w_{xx}$','interpreter','latex');
subplot(5,1,5);plot(z,T5);axis('tight');ylabel('$\nu\,w_{zz}$','interpreter','latex');
subplot(5,1,1);
legend(h,{'left','right'},'location','best');
title(sprintf('nm: %s, dt: %6.2f, $M_z$',nm,T),'interpreter','latex');
subplot(5,1,5);xlabel('$z$','interpreter','latex');

T1=[mean(       Izb(:,ff1),2),mean(       Izb(:,ff2),2)];T1=T1-repmat(T1(1,:),size(T1,1),1);
T2=[mean(  0.5*a.ww(:,ff1),2),mean(  0.5*a.ww(:,ff2),2)];T2=T2-repmat(T2(1,:),size(T2,1),1);
T3=[mean(  a.udwdxIz(:,ff1),2),mean(  a.udwdxIz(:,ff2),2)];T3=T3-repmat(T3(1,:),size(T3,1),1);
T4=[mean(nu*a.dwdxxIz(:,ff1),2),mean(nu*a.dwdxxIz(:,ff2),2)];T4=T4-repmat(T4(1,:),size(T4,1),1);
T5=[mean(  nu*a.dwdz(:,ff1),2),mean(  nu*a.dwdz(:,ff2),2)];T5=T5-repmat(T5(1,:),size(T5,1),1);;

figure;
clf;
subplot(5,1,1);h=plot(z,T1);axis('tight');ylabel('$\int\,M_z\,dz$','interpreter','latex');
subplot(5,1,2);plot(z,T2);axis('tight');ylabel('$\int\,ww/2\,dz$','interpreter','latex');       
subplot(5,1,3);plot(z,T3);axis('tight');ylabel('$\int\,uw_x\,dz$','interpreter','latex');       
subplot(5,1,4);plot(z,T4);axis('tight');ylabel('$\int\,\nu\,w_{xx}\,dz$','interpreter','latex');
subplot(5,1,5);plot(z,T5);axis('tight');ylabel('$\nu\,w_{z}|$','interpreter','latex');
subplot(5,1,1);
legend(h,{'left','right'},'location','best');
title(sprintf('nm: %s, dt: %6.2f, $\\int M_z\\,dx$',nm,T),'interpreter','latex');
subplot(5,1,5);xlabel('$z$','interpreter','latex');


if isfield(a,'vv')
  TE=max(0,a.uu-a.u.^2+a.vv+a.ww-a.w.^2)/2;
else
  TE=max(0,a.uu-a.u.^2+a.ww-a.w.^2)/2;
end
vy=a.dwdx-a.dudz;
for k=1:2
  figure;clf;  
  sp=get(0,'screensize');
  set(gcf,'units','pixels');
  fp=get(gcf,'position');
  set(gcf,'position',[10 sp(4)-650 780 590],'defaultaxesfontsize',8);
  
  for j=1:5;
    aa(j)=axes('units','pixels','position',[30 30+(5-j)*110 700 100]);
  end
  
  switch(k)
    case 1
      if isfield(a,'S');axes(aa(1));imagesc(x,a.z,a.S(:,fx)); ylabel('S'); bh=jcolorbar;end
      axes(aa(2));imagesc(x,a.z,TE(:,fx));  ylabel('TE');jcolorbar;
      if isfield(a,'Wx');axes(aa(3));imagesc(x,a.z,a.Wx(:,fx));ylabel('Wx');jcolorbar;end
      if isfield(a,'Wy');axes(aa(4));imagesc(x,a.z,a.Wy(:,fx));ylabel('Wy');jcolorbar;end
      if isfield(a,'Wz');axes(aa(5));imagesc(x,a.z,a.Wz(:,fx));ylabel('Wz');jcolorbar;end
    case 2
      axes(aa(1));imagesc(x,a.z,a.b(:,fx));ylabel('b');set(aa(1),'clim',[0 1]);jcolorbar;
      hold('on');C=contourc(x,a.z,a.b(:,fx),0.05:0.1:0.95);hc=plot_contours(C);set(hc,'color',[1 1 1]);
      axes(aa(2));imagesc(x,a.z,a.p(:,fx)-mean(a.p(:,end)));ylabel('p');jcolorbar;
      hold('on');C=contourc(x,a.z,a.p(:,fx),[-1 -0.5 0 0.5]);hc=plot_contours(C);set(hc,'color',[1 1 1]);
      axes(aa(3));imagesc(x,a.z,a.u(:,fx));ylabel('u');jcolorbar; 
      hold('on');C=contourc(x,a.z,a.u(:,fx),[-2 -1 0]);hc=plot_contours(C);set(hc,'color',[1 1 1]);
      axes(aa(4));imagesc(x,a.z,a.w(:,fx));ylabel('w');jcolorbar;
      hold('on');C=contourc(x,a.z,a.w(:,fx),[-0.1 0 0.1]);hc=plot_contours(C);set(hc,'color',[1 1 1]);
      axes(aa(5));imagesc(x,a.z,vy(:,fx));ylabel('vy');jcolorbar;
      hold('on');C=contourc(x,a.z,vy(:,fx),[-1 0 1]);hc=plot_contours(C);set(hc,'color',[1 1 1]);
  end
  set(aa(1:4),'ydir','normal','xticklabel',[],'yticklabel',[]);
  set(aa(5),  'ydir','normal','yticklabel',[]);
end
a.fx=fx;
a.T=T;
a.t1=t1;
a.t2=t2;

return

