function ded_gc(n)
if nargin<1
  fns=sort(cellstr_ls('~/data/dedalus/gc/*/*.h5'));
  fn=fns{end};
else
  fn=sprintf('~/data/dedalus/gc/gc_s%i/gc_s%i_p0.h5',n,n);
end

if ~isfile(fn)
  disp(fn);
  return;
end
disp(fn);
%h5disp(fn);
%a=h5info(fn);
try
  a.x=h5read(fn,'/scales/x/2');
  a.y=h5read(fn,'/scales/y/2');
catch
  a.x=h5read(fn,'/scales/x/1');
  a.y=h5read(fn,'/scales/y/1');
end
a.dt=h5read(fn,'/scales/timestep');
a.t=h5read(fn,'/scales/sim_time');
a.i=h5read(fn,'/scales/iteration');
a.b=h5read(fn,'/tasks/b');
a.u=h5read(fn,'/tasks/u');
a.v=h5read(fn,'/tasks/v');

a.nm=h5readatt(fn,'/','set_number');

H=round(10*a.y(end))/10;
L=round(10*a.x(end))/10;
S=1024/L;
nx=length(a.x);
nt=length(a.t);
HT=min(2,H);

szxt=[nx nt];

a.B=reshape(h5read(fn,'/tasks/B'),szxt);
a.P=reshape(h5read(fn,'/tasks/P'),szxt);
a.UU=reshape(h5read(fn,'/tasks/UU'),szxt);
a.UV=reshape(h5read(fn,'/tasks/UV'),szxt);
a.VV=reshape(h5read(fn,'/tasks/UV'),szxt);
a.Q=reshape(h5read(fn,'/tasks/Q'),szxt);

szt=[nt 1];
a.B00=reshape(h5read(fn,'/tasks/B00'),szt);
a.B01=reshape(h5read(fn,'/tasks/B01'),szt);
a.B10=reshape(h5read(fn,'/tasks/B10'),szt);
a.B11=reshape(h5read(fn,'/tasks/B11'),szt);
a.B20=reshape(h5read(fn,'/tasks/B20'),szt);
a.B02=reshape(h5read(fn,'/tasks/B02'),szt);


figure(1);
clf;
a1=axes('position',[0.01 0.69 0.9 0.3],'linewidth',p.lw);
a2=axes('position',[0.01 0.38 0.9 0.3],'linewidth',p.lw);
a3=axes('position',[0.01 0.07 0.9 0.3],'linewidth',p.lw,'FontSize',8);
set(gcf,'units','pixels');
pf=get(gcf,'position');
set(gcf,'position',[pf(1:2),round(S*L),3*round(S*HT)]);
NY=length(a.y);
ny=round(H/L*length(a.x));
y=linspace(a.y(1),a.y(end),ny);


xp=([0.5:15.5])';

maxB =max(a.B(:));
maxUV=max(abs(a.UV(:)));
maxP =max(abs(a.P(:)));
maxQ =max(abs(a.Q(:)));

vs=0.4;vt=0.5*vs;
NULLxy=repmat(NaN,[ny,nx]);
NULLx =repmat(NaN,[nx,1]);

axes(a1);hb=imagesc(a.x,y,NULLxy,[ 0   1  ]); c1=colorbar('position',[0.93 0.69 0.02 0.3],'yticklabel',[]);
axes(a2);hu=imagesc(a.x,y,NULLxy,[-1.2 0.5]); c2=colorbar('position',[0.93 0.38 0.02 0.3],'xtick',[-1 0]);
set([a1 a2],'ydir','normal','dataaspectratio',[1 1 1],'xlim',[0 L],'ylim',[0 HT],'ytick',[],'xtick',[]);
axes(a1);ht=text(L/200,HT,'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
axes(a2);
line(repmat([xp'-vt xp' xp'+vt],2,1),repmat([0;HT],1,3*length(xp)),'color',[0.7 0 0]);
hp=line(repmat(xp,1,NY)+NaN,a.y,'color',[0 0 1]);
axes(a3);
set(a3,'xlim',[0 L],'ylim',[-1 1]);
hB=line(a.x,NULLx,'color',[1 0 0]);
hU=line(a.x,NULLx,'color',[0 0 1]);
hP=line(a.x,NULLx,'color',[0.7 0.7 .7]);
hQ=line(a.x,NULLx,'color',[0   0.5 0]);
hL=legend([hB hU hP hQ],{'B','UV','P','Q'},'FontSize',8,'position',[0.93 0.01 0.067 0.32]); 

for j=1:length(a.t)
  b=interp1(a.y,a.b(:,:,j),y,'linear');
  u=interp1(a.y,a.u(:,:,j),y,'linear');
  pu=vs*interp1(a.x,permute(a.u(:,:,j),[2 1 3]),xp,'linear');
  ts=sprintf('%u: t=%6.2f, u=%6.3f, H=%3.1f',a.nm,a.t(j),-mean(mean(u(:,a.x<=1))),H);
  set(hb,'CData',b);
  set(hu,'CData',u);
  set(ht,'string',ts);
  set(hB,'ydata',a.B(:,j)/maxB);
  set(hU,'ydata',a.UV(:,j)/maxUV);
  set(hP,'ydata',a.P(:,j)/maxP);
  set(hQ,'ydata',a.Q(:,j)/maxQ);
  for j=1:length(hp); set(hp(j),'xdata',xp(j)+pu(j,:));end
  drawnow;
  pause(0.1);
end

figure(2);
subplot(2,1,1);
plot(a.x,a.b(1,:,end));
axis([0 16 0 1]);
ylabel('B');
subplot(2,1,2);
h=plot(a.x,a.Q(:,end),a.x,a.u(1,:,end),a.x,a.u(end,:,end));
legend(h,{'Q','u0','u1'});
axis([0 16 -inf inf]);


M=a.B00./max(a.B00);
X=a.B10./a.B00/L;
Y=a.B01./a.B00/H*2;
XX=sqrt(a.B20./a.B00)/L;
YY=sqrt(a.B02./a.B00)/H*2;

figure(3);
h=plot(a.t,M/max(M),a.t,X/max(X),a.t,Y/max(Y),a.t,XX/max(XX),a.t,YY/max(YY));
legend(h,{'M','X','Y','XX','YY'},'location','best');
axis([a.t(1) a.t(end) -inf 1]);

return;

tol=0.02;
uu=(abs(u+0.0)<tol)+2*(abs(u+0.5)<tol)+3*(abs(u+1.0)<tol);


b=a.b(:,:,end);



figure(2);
clf
dx=(a.x(end)-a.x(1))/(length(a.x)-1);
dy=(  y(end)-  y(1))/(length(y)-1);
j=length(a.t);
u=interp1(a.y,a.u(:,:,j),y,'linear');
v=interp1(a.y,a.u(:,:,j),y,'linear');
b=interp1(a.y,a.b(:,:,j),y,'linear');
subplot(3,1,3);
f=helmholtz_fd(u,v,dx,dy,'fourier');
subplot(3,1,1);

uu=repmat(0,size(u));
tol=0.02;
uu(abs(u)<tol)=1;
uu(abs(u+0.5)<tol)=2;
uu(abs(u+1)<tol)=3;
imagesc(a.x,y,uu);

contour(a.x,y,u,[(-0.02:0.01:0.02) (-0.52:0.01:-0.48) -(0.98:0.01:1.02)]);% -0.02:0.01:0.02]);
axis([0 5 0 1]);
subplot(3,1,2);
contour(a.x,y,b,[0:0.1:0.9 0.98:0.01:1.02]);
axis([0 5 0 1]);
subplot(3,1,3);
contour(a.x,y,f,50);



rg=round(linspace(1,size(u,2),64));
uu=u(:,rg);
uu=uu/max(abs(uu(:)));
plot(repmat(a.x(rg)',ny,1)+uu,y);
axis([0 16 0 1]);

