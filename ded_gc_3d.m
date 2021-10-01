function [p a b c d]=ded_gc_3d(nm)
if nargin<2
  quick=0;
end


nm=['gc/' nm];

p=ded_read_param(nm); b=ded_read_mom(nm);
d=ded_read_g(nm,'yz');
a=ded_read_yavg(nm);
g=ded_read_g(nm,'gc');

U=p.U;
if ~isfield(p,'B')
  p.B=1;
end

tol=1e-3;
switch(p.version)
  case '1.30'
    p.B=1;
    if isfield(d,'u')
      d.u   =   d.u/(p.U);
    end
    d.b   =   d.b/(p.B);
    d.uu  =  d.uu/(p.U);
    d.e   =   NaN;
    d.s   =   NaN;
    a.b   = a.b/p.B;
    a.u   = a.u/p.U;
  case '1.31'
    
    n=min(length(d.t),length(g.t));
    d.t=d.t(1:n);
    g.t=g.t(1:n);
    dnm=setdiff(fieldnames(d),{'t','x','y','z','nx','ny','nz','nt','sz'});
    for j=1:length(dnm)
      d.(dnm{j})=d.(dnm{j})(:,1:n);
    end   
    dnm=setdiff(fieldnames(g),{'t','x','y','z','nx','ny','nz','nt','sz'});
    for j=1:length(dnm)
      g.(dnm{j})=g.(dnm{j})(:,1:n);
    end   

    Zhb   =     max(0,d.b)/(p.W*p.B);
    Zhbz  = max(0,2*g.bz)./max(tol,d.b);
    Zhbb  =        d.b.^2./max(tol,g.bs);
    d.ev  =  d.bw./max(tol,d.b);
    d.ev  = d.ev/max(abs(d.ev(:)));
    d.bu  =  d.bu/(p.W*p.h*p.U1/2);
    d.bu  =  d.bu/max(abs(d.bu(:)));
    d.u   =   d.u/(p.W*p.h*p.U1);
    d.uu  =  d.uu/(p.W*p.h*p.U );
    a.b   =   a.b/(p.W*p.B);
    a.u   =   a.u/(p.W*p.h*p.U);
    d.E   =   d.E/max(d.E(:));
    d.S   =   d.S/max(d.S(:));
end

S=1024/p.L;
HT=min(2,p.H);
nz=round(p.H/p.L*p.Nx);
z=filter_midpoint(linspace(0,p.H,nz+1));

T=p.T;


figure(1);
clf;
a1=axes('position',[0.01 0.79 0.9 0.20]);
a2=axes('position',[0.01 0.58 0.9 0.20]);
a3=axes('position',[0.05 0.10 0.9 0.45],'FontSize',8);
set(gcf,'units','pixels');
pf=get(gcf,'position');
set(gcf,'position',[pf(1:2),round(S*p.L),4*round(S*HT)]);

xp=([1:ceil(p.L)-1])';


vs=0.4;vt=vs;
bs=0.4;
NULLxz=repmat(NaN,[nz,p.Nx]);
NULLx =repmat(NaN,[p.Nx,1]);

if isempty(a)
  return;
end

axes(a1);
hb=imagesc(a.x,z,NULLxz,[ 0   1  ]); 
set(a1,'clim',[0 1]);
c1=colorbar('position',[0.93 0.79 0.02 0.2],'ytick',[0 1]);
axes(a2);hu=imagesc(a.x,z,NULLxz,[-1.2 0.5]); c2=colorbar('position',[0.93 0.58 0.02 0.2],'xtick',[-1 0]);
set([a1 a2],'ydir','normal','dataaspectratio',[1 1 1],'xlim',[0 p.L],'ylim',[0 HT],'ytick',[],'xtick',[]);
axes(a1);ht=text(p.L/200,HT,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
line(repmat([xp'-bs/2 xp'+bs/2],2,1),repmat([0;HT],1,2*length(xp)),'color',[0.7 0 0]);
hr=line(repmat(xp,1,p.Nz)+NaN,a.z,'color',[0 0 1],'linewidth',2);

axes(a2);
line(repmat([xp'-vt xp' xp'+vt],2,1),repmat([0;HT],1,3*length(xp)),'color',[0.7 0 0]);
hp=line(repmat(xp,1,p.Nz)+NaN,a.z,'color',[0 0 1],'linewidth',2);
nl=5;
for j=1:nl
  axes(a2);
  for k=1:4
    hc1(k,j)=line([0 0],[0 0],'color',[1 1 1]);
  end
  axes(a1);
  for k=1:9
    hc2(k,j)=line([0 0],[0 0],'color',[1 1 1]);
  end
end
axes(a3);cla;
set(a3,'xlim',[0 p.L],'ylim',[-1 1],'box','on');

hh1=line(a.x,NULLx,'color',[1 0 0],'LineStyle','-');
hh2=line(a.x,NULLx,'color',[0 0 1],'LineStyle','-');
hh3=line(a.x,NULLx,'color',[0 1 0]/2,'LineStyle','-');
hh4=line(a.x,NULLx,'color',[0 1 1]/3,'LineStyle','-');
hh5=line(a.x,NULLx,'color',[1 1 0]/3,'LineStyle','-');
% $$$ hh6=line(a.x,NULLx,'color',[0   0 0,'LineStyle','--']);
line([0 p.L],[0 0],'color',0.7*[1 1 1]);
line([0;p.L],[[1;1] [-1;-1] [1;1]/2],'color',0.7*[1 1 1],'LineStyle','--');

%legend([hh1 hh2 hh3 hh4 hh5 hh6],{'hb','hbz','hbb','p','bu','up'},'location','eastoutside');
legend([hh1 hh2 hh3 hh4 hh5],{'b','z','bb','bu','ev'},'location','eastoutside');

if isempty(g)
  mina=min([d.b(:);d.p(:)]);
  maxa=max([d.b(:);d.p(:)]);
else
  mina=min([d.b(:);d.p(:);d.bu(:);d.ev(:);g.bz(:)]);
  maxa=max([d.b(:);d.p(:);d.bu(:);d.ev(:);g.bz(:)]);
end
rga=(mina+maxa)/2+1.05*(maxa-mina)/2*[-1 1];


tt=unique([a.t;d.t]);
tt=unique([a.t]);

fao=0;
fbo=0;
for jj=1:length(tt)
  fa=findmin(abs(tt(jj)-a.t));
  fb=findmin(abs(tt(jj)-d.t));
  ts=sprintf('%s: t=%7.3f, U=%5.3f, H=%3.1f, W=%3.1f',p.name,tt(jj),U,p.H,p.W);  set(ht,'string',ts);
  
  if fa~=fao
    fao=fa;
    j=fa;
    b=interp1(a.z,a.b(:,:,j),z,'linear');
    u=interp1(a.z,a.u(:,:,j),z,'linear');
    pu=vs*interp1(a.x,permute(a.u(:,:,j),[2 1 3]),xp,'linear');
    pb=bs*(interp1(a.x,permute(a.b(:,:,j),[2 1 3]),xp,'linear')-0.5);

    [cx,cy]=longest_contours(a.x,z,u,-0.99, nl); for k=1:length(cx);set(hc1(1,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,u,-1.01, nl); for k=1:length(cx);set(hc1(2,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,u,-0.01, nl); for k=1:length(cx);set(hc1(3,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,u, 0.01, nl); for k=1:length(cx);set(hc1(4,k),'Xdata',cx{k},'Ydata',cy{k});end
    %    [cx,cy]=longest_contours(a.x,z,b, 0.999,nl); for k=1:length(cx);set(hc5(k),'Xdata',cx{k},'Ydata',cy{k});end
    %[cx,cy]=longest_contours(a.x,z,b, 0.99, nl); for k=1:length(cx);set(hc6(k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.9,  nl); for k=1:length(cx);set(hc2(1,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.8,  nl); for k=1:length(cx);set(hc2(2,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.7,  nl); for k=1:length(cx);set(hc2(3,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.6,  nl); for k=1:length(cx);set(hc2(4,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.5,  nl); for k=1:length(cx);set(hc2(5,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.4,  nl); for k=1:length(cx);set(hc2(6,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.3,  nl); for k=1:length(cx);set(hc2(7,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.2,  nl); for k=1:length(cx);set(hc2(8,k),'Xdata',cx{k},'Ydata',cy{k});end
    [cx,cy]=longest_contours(a.x,z,b, 0.1,  nl); for k=1:length(cx);set(hc2(9,k),'Xdata',cx{k},'Ydata',cy{k});end
    set(hb,'CData',b);
    set(hu,'CData',u);
    for k=1:length(hp); set(hp(k),'xdata',xp(k)+pu(k,:));end
    for k=1:length(hr); set(hr(k),'xdata',xp(k)+pb(k,:));end
  end
  if fb~=fbo
    fbo=fb;
    j=fb;
    set(hh1,'ydata', Zhb(:,j));
    set(hh2,'ydata', Zhbz(:,j));
    set(hh3,'ydata', Zhbb(:,j));
    %set(hh4,'ydata', d.p(:,j));
    set(hh4,'ydata', d.bu(:,j));
    set(hh5,'ydata', d.ev(:,j));
  end
  drawnow;
  %  pause(0.1);
end

return







figure(2);
%subplot(2,1,1);
plot(a.x,a.b(1,:,end));
axis([0 16 0 1]);
zlabel('B');
subplot(2,1,2);
h=plot(a.x,a.Q(:,end),a.x,a.u(1,:,end),a.x,a.u(end,:,end));
legend(h,{'Q','u0','u1'});
axis([0 16 -inf inf]);


return;

tol=0.02;
uu=(abs(u+0.0)<tol)+2*(abs(u+0.5)<tol)+3*(abs(u+1.0)<tol);


b=a.b(:,:,end);



figure(2);
clf
dx=(a.x(end)-a.x(1))/(length(a.x)-1);
dz=(  z(end)-  z(1))/(length(z)-1);
j=length(a.t);
u=interp1(a.z,a.u(:,:,j),z,'linear');
v=interp1(a.z,a.u(:,:,j),z,'linear');
b=interp1(a.z,a.b(:,:,j),z,'linear');
subplot(3,1,3);
f=helmholtz_fd(u,v,dx,dz,'fourier');
subplot(3,1,1);

uu=repmat(0,size(u));
tol=0.02;
uu(abs(u)<tol)=1;
uu(abs(u+0.5)<tol)=2;
uu(abs(u+1)<tol)=3;
imagesc(a.x,z,uu);

contour(a.x,z,u,[(-0.02:0.01:0.02) (-0.52:0.01:-0.48) -(0.98:0.01:1.02)]);% -0.02:0.01:0.02]);
axis([0 5 0 1]);
subplot(3,1,2);
contour(a.x,z,b,[0:0.1:0.9 0.98:0.01:1.02]);
axis([0 5 0 1]);
subplot(3,1,3);
contour(a.x,z,f,50);



rg=round(linspace(1,size(u,2),64));
uu=u(:,rg);
uu=uu/max(abs(uu(:)));
plot(repmat(a.x(rg)',nz,1)+uu,z);
axis([0 16 0 1]);

