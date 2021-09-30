function [a b]=ded_plume_plot(nm,typ,mfn,fnt)

if nargin<2
  typ=[];
end
if nargin<3
  mfn=[];
end
if nargin<4
  fnt=[];
end
if isempty(typ)
  typ={'C','T','S','u','w'};
end

dpi=get(0,'ScreenPixelsPerInch');

if isempty(fnt)
  fnt=12*72/dpi;
end

nf=length(typ);

if 0
  figure(1);clf
  b=ded_plume_plot_forcing(nm);
end

p=ded_read_param(['plume/' nm]);
L=p.L;
H=p.H;

a=ded_plume_read_yavg(nm);
nx=length(a.x);
nz=round(2*p.Nz);
z=filter_midpoint(linspace(0,p.L,nz+1));


figure(1);
clf;
ylim(1,:)=[0 p.H];
if isfield(p,'IHeight');
  ylim(end+1,:)=[0 p.H] + p.IHeight*[1 0];
end
if isfield(p,'Buffer');
  ylim(end+1,:)=[0 p.H] + p.Buffer/2*[1 -1];
end
ylim=[max(ylim(:,1))  min(ylim(:,2))];

xlim=[0 L];
f2=find(a.x>=ylim(1) & a.x<=ylim(2));
f1=find(a.z>=xlim(1) & a.z<=xlim(2));

if L>2*H
  asz=[[1 1];[1 2];[1 3];[1 4];[1 5];[2 3];[2 4];[2 4]];
else
  asz=[[1 1];[1 2];[1 3];[2 2];[2 3];[2 3];[2 4];[2 4]];
end

asz=asz(nf,:);

ta=axes('position',[0 0 1 1],'visible','off','xlim',[-1 1],'ylim',[0 1]);

Q=p.IFlux;
IV=p.IFlux/p.IWidth;
IRe=IV*p.Re*p.IWidth;
ts=sprintf('Re=%5.0f, Nx=%u, Ny=%u, Nz=%u, L=%3.1f, H=%3.1f, W=%3.1f, V=%4.2f\n IRe=%5.1f, Q=%4.2f, IC=%4.2f, IT=%4.2f, IS=%4.2f, IW=%4.2f, IV=%4.2f',...
           p.Re,p.Nx,p.Ny,p.Nz,p.L,p.H,p.W,p.V,IRe,Q,p.IC,p.IT,p.IS,p.IWidth,IV);
htt=text(0,1,ts,'horizontalalignment','center','verticalalignment','top','color',0*[1 1 1]);

NN=1;
aa=jsubplot(asz,[0.02 0.01],[0.01 0.01],[0.02 0.05]);
delete(aa(nf+1:end));
aa(nf+1:end)=[];
%axes(aa(1));hc=imagesc(z,a.x,a.T(:,:,1)',[minT maxT]);
for j=1:nf
  X=a.(typ{j});
  minx(j)=min(min(min(X(f1,f2,:))));
  maxx(j)=max(max(max(max(X(f1,f2,:)))),minx(j)+1e-3);
  axes(aa(j));
  hc(j)=imagesc(z,a.x,NN*X(:,:,1)',[minx(j) maxx(j)]);
  c(j)=colorbar;
  ht(j)=text(xlim(1),ylim(2),'','horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
end
delete(c);

%set(aa(:,1:end-1),'xticklabel',[]);
%set(aa(2:end,:),'yticklabel',[]);
set(aa,'xtick',xlim,'ytick',ylim);
set(aa,'xtick',[],'ytick',[]);

set(aa,'xlim',xlim,'ylim',ylim)
set(aa,'box','on','ydir','normal');
set(aa,'DataAspectRatio',[1 1 1]);
set(aa,'FontSize',fnt);
set(ht,'BackgroundColor',[0 0 0 .8]);


FX=(0.2+asz(1))*L;
FY=(0.2+asz(2))*H ;
if FX>FY
  S=1024/FX;
else
  S=768/FY;
end
FX=8*round(S*FX/8);
FY=8*round(S*FY/8);
set(gcf,'units','pixels','menubar','none');
szp=get(gcf,'position');
szp(3:4)=[FX FY];
set(gcf,'position',szp);
set(gcf,'units','inches');
szi=get(gcf,'position');
szp=round(szp(3:4));
szi=szi(3:4);
if ~isempty(mfn)
  preprint(szi,fnt);colour_lines;
  %adjust_paper_size(szp,dpi); % Check that we get the right size
  fnmp4=[mfn '/plume-' p.name '.mp4'];
  fnyuv=[mfn '/plume-' p.name '.yuv'];
  fp=fopen(fnyuv,'w');
end
psz=get(gcf,'paperposition');psz=psz(3:4);
disp(sprintf('Figure size %7.3fx%7.3f pixels',szp(1),szp(2)));
disp(sprintf('Figure size %7.3fx%7.3f inches',szi(1),szi(2)));
disp(sprintf('Paper  size %7.3fx%7.3f inches',psz(1),psz(2)));
disp(sprintf('root dpi =%7.3f, measured dpi=%7.3f',dpi, szp/szi))
drawnow;
sz=[];
esf=[1 1];
for j=1:length(a.t);
  set(htt,'string',sprintf('%s, t=%7.2f',ts,a.t(j)));
  for k=1:nf
    X=interp1(a.z,a.(typ{k})(:,:,j),z,'linear');
    set(ht(k),'string',sprintf('%s=%6.4f',typ{k},mean(X(:))));
    set(hc(k),'CData',X');
  end
  if fp==0
    drawnow;
    pause(.1);
  else
    gu_write_frame(fp,szp,dpi)
  end
end
if ~isempty(mfn)
  fclose(fp);
  fps=5;
  make_mp4(fnyuv,fnmp4,szf,fps);%,'yuv444p');
end
a.aa=aa;


return;

figure(3);
nt=size(a.S,3);
mw=round(max(a.wg(:)));
mw=100;
for j=1:nt
  subplot(2,1,1);plot(a.x,mean(a.S(:,:,j),1),'b',a.x,mean(a.fs(:,:,j),1),'b--',a.x,mean(a.wg(:,:,j),1)/mw,'r');
  subplot(2,1,2);plot(a.x,mean(a.T(:,:,j),1),'b',a.x,mean(a.fT(:,:,j),1),'b--',a.x,mean(a.wg(:,:,j),1)/mw,'r');
  pause(0.1);
end
return;


for j=6:6
  ded_plume_plot(sprintf('%03i',j),{'C','S','T'},'/tmp/a',12);
end

for j=1:12;  ded_plume_plot(sprintf('%03i',j),{'C','S','T'},'/tmp/ham',8);end
