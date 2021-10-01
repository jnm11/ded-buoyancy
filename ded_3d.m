function ded_3d(nm,p,force)
%p.sz=[512 1024]/4;p.aa=1;ded_3d('pm/test','~/pm/test/a',p);
%p.sz=[320 768];ded_3d('pm/test','~/pm/test/a',p);
%p.sz=[320 512];ded_3d('pm/019','~/pm/019/a',p);
%p.sz=[320 768];p.aa=4;ded_3d('pm/011','~/pm/011/a',p);
%p.sz=[640 320];p.aa=1;ded_3d('gc/test','~/gc/test/a',p); 
%p.sz=[640 320];p.aa=1;ded_3d('gc/pos/f6/084','~/gc/mv/084/a',p);
%p.sz=2*[640 320];p.aa=4;ded_3d('gc/pos/f6/084','~/gc/mv/084/a',p);%p.sz=2*[640 320];p.aa=4;p.rgx=[-7 1];ded_3d('gc/pos/f6/064','~/films/Dedalus/gc-pos-064',p);p.cva=3.5;
%p.sz=[512 768];p.aa=2;p.sz=[512 768];p.ctar=0.4;p.cva=3.1;ded_3d('pm/qpmf0','~/films/Dedalus/pm-qpmf0',p)
%nm='pm/020';p=ded_3d_param(nm);p.ftype='s';ded_3d(nm,p);p.ftype='b';ded_3d(nm,p);
%nm='pm/ms/001';p=ded_3d_param(nm);ded_3d(nm,p,1);
DD=[ded_dedalus_data_dir '/'];  
if nargin<3
  force=0;
end

if nargin==0
  nms=findname(DD,'param.h5');
  nms=cellstrremove(nms,'/param.h5');
  nms=cellstrremoveprefix(nms,DD);
  nms=cellstrremoveprefix(nms,'/');
  nms=sort(nms);
  for j=1:length(nms)
    ded_3d(nms{j});
  end
  return;
end
if nargin==1
  if ~iscell(nm)
    nms=ded_regexp_names(nm);
    if ~isempty(nms)
      nm=cellstrremoveprefix(nms,[ded_dedalus_data_dir '/']);
    end
  end
  if ~iscell(nm)
    nms={nm};
  else
    nms=nm;
  end
  nms=sort(nms);
  for j=1:length(nms)
    nm=nms{j};
    p=ded_3d_param(nm);
    sp=ded_read_param(nm);
    if isempty(p)
      %disp(sprintf('ded_3d: no parameters for %s',nm));
      nms{j}=[];
      continue;
    end
    typs={'b','s','d','u','v','w'};
    typs={'b'};
    for k=1:length(typs)
      if isdir([DD '/' nm '/' typs{k}]) 
        p.ftype=typs{k};
        switch(typs{k})
          case 'b'
            p.isoval=0.05*abs(sp.B);
          case 's'
            p.isoval=0.05*abs(sp.S);
          case {'u','v','w','d'}
            if isfield(sp,'hb')
              p.isoval=0.05*abs(sqrt(sp.g*sp.hb));
            else
              p.isoval=0.05*sp.V;
            end
        end
        ded_3d(nm,p);
      end
    end
  end
  return;
end


pd.text='true';
pd.reduce=0;
pd.ftype='b';
pd.nmbnm=[];
pd.xtick=[];
pd.ytick=[];
pd.ztick=[];

p=combine_struct(pd,p);

[nmb t]=ded_fix_time(nm,p.ftype);

bnm=p.ftype;

if isempty(p.nmbnm)
  p.nmbnm=[nm '-' bnm];
end

fnmp4=['~/films/Dedalus/' p.nmbnm '.mp4'];
dd=['~/tmp/Dedalus/' p.nmbnm];
if ~isdir(dd)
  mkdir(dd);
end
fnyuv=[dd '.yuv'];
fnffmpeg=[dd '.ffmpeg'];
pnpp=[dd '.mat'];


if isempty(nmb)
  disp(sprintf('ded_3d: No density files found "%s"',nm));
  return;
end
%if ~isfile(nmb{1})
%  disp(sprintf('ded_3d: Not a file "%s"',nmb{1}));
%  return;
%end
if isfile(fnmp4)
  if all(file_nt(fnmp4,nmb))
    disp(sprintf('ded_3d: "%s" mp4 is up to date, %s',nm,fnmp4));
    return;
  end
end

disp(sprintf('ded_3d: Making  %s',nm));

fh=figure(1);clf
set(gcf, 'Windowstyle', 'normal')

fp=fopen(fnyuv,'w');
if fp<0
  error(sprintf('ded_3d: Failed to open "%s"',fnyuv));
end

name=nm;
nmin=nm;

a=ded_read_hdf(nmb{1});
if isempty(a)
  return;
end

b=ded_read_param(nm);
c=ded_coord(nm);
if isempty(c)
  return;
end
if isfield(a,'b');sz=size(a.b);end;
if isfield(a,'s');sz=size(a.s);end;
if isfield(a,'u');sz=size(a.u);end;
if isfield(a,'v');sz=size(a.v);end;
if isfield(a,'w');sz=size(a.w);end;

if sz(c.dim)==c.NAx
  c.Jx=c.Ax;
  c.NJx=c.NAx;
end
if c.dim==3
  if sz(2)==c.NAy
    c.Jy=c.Ay;
    c.NJy=c.NAy;
  end
end
if sz(1)==c.NAz
  c.Jz=c.Az;
  c.NJz=c.NAz;
end

YFR=false;
if c.dim==2
  y=[0;b.W];
else
  YFR = strcmp(b.T{2},'Fourier');
  YSC = strcmp(b.T{2},'SinCos');
  y=c.Jy(:);
  if YFR y(end+1)=b.W;end;
end
x=c.Jx(:);
z=c.Jz(:);


if ~isfield(p,'transparent')
  p.transparent=0;
end

%f=ded_read_g(nm,'force');
bnml=char(bnm+'A'-'a');
if isfield(b,'bnml');
  fnm=b.(bnml);
else
  fnm=1;
end
if isfield(b,'inlet')
  switch(b.inlet)
    case 'gaussian'
      fnm=fnm*2;
  end
end

if b.Ny>1
  ry=0;
else
  a.y=[0 b.W];
  ry=1;
end

if isfield(p,'copyx')
  if p.copyx
    x=[x;x+b.L];
  end
end
if isfield(p,'flipy')
  if p.flipy
    y=[y-b.W;y];
    p.rgy(1)=-p.rgy(2);
    p.y0=0;
    b.W=2*b.W;
  end
end
if isfield(p,'flipz')
  if p.flipz
    z=[z-b.H;z];
    p.rgz(1)=-p.rgz(2);
    b.H=2*b.H;
    p.z0=0;
  end
end


switch b.sType
  case {'pm','plume'}
    arg=[p.rgz p.rgy p.rgx-p.x0];
    zz=z;
    z=x-p.x0;
    y=y-p.y0;
    x=zz-p.z0;
    zz=[];
  case 'gc'
    arg=[p.rgx p.rgy p.rgz];
    zz=linspace(0,b.H,round(b.Nx*b.H/b.L));
    x=x-p.x0;
    y=y;
    z=zz;
end

rgx=arg(1:2);
rgy=arg(3:4);
rgz=arg(5:6);

Lx=diff(rgx);
Ly=diff(rgy);
Lz=diff(rgz);


[szfp szfi szpi dpi]=gu_film_setup(p.sz,p.aa,p.rd);

set(fh,'inverthardcopy','off');
set(fh,'defaulttextfontsize',8);
set(fh,'clipping','on');

if ~isfield(p,'lw')
  p.lw=p.aa;
end

axt=axes('position',[0 0 1 1],'color','none','visible','off','xlim',[0 1],'ylim',[0 1],'linewidth',p.lw);
ss=ded_read_stats(nmin);
if ~isfinite(p.tsn) 
  p.tsn=0;
  while 1
    ts=ded_interp_stats(ss,b,[],p.tsn);
    if isempty(ts)
      return;
    end
    nts=zeros(size(ts));
    for j=1:length(ts)
      nts(j)=length(ts);
    end
    j=findmax(nts);
    the=[1 1 0 0];
    th=text(1,1,ts{j},'fontsize',10,'fontweight','bold','HorizontalAlignment','right','verticalalignment','top');
    thej=get(th,'extent');
    the=[min(the(1),thej(1)) max(the(2),thej(2)) min(the(3),thej(3)) max(the(4),thej(4))];  
    if p.tsn>5
      break;
    end
    %patch(the(1)+[0 0 the([3 3])], the(2)+[0 the([4 4]) 0],[0 0 0])
    delete(th);
    if the(1)>=0 & the(3)<=1
      break;
    end
    p.tsn=p.tsn+1;
  end
end

ts=ded_interp_stats(ss,b,0,p.tsn);
th=text(1,1,ts{1},'fontsize',10,'fontweight','bold','HorizontalAlignment','right','verticalalignment','top');
the=get(th,'extent');
%patch(the(1)+[0 0 the([3 3])], the(2)+[0 the([4 4]) 0],[0 0 0])
set(th,'String','');

ax=axes('position',[0 0 1 1],'clipping','on','box','off','linewidth',p.lw);
daspect([1 1 1]);
hold('on');
set(ax,'ActivePositionProperty','position');
set(ax,'Projection','orthographic');

switch b.sType
  case {'pm','plume'}
    cpos=[-10 -10 3]*b.L;
    ctar=[0 0 (b.L-p.x0)/2*p.ctar+(1-p.ctar)*(b.L+p.x0)/2];
  case 'gc'
    ctar=[mean(rgx) mean(rgy) mean(rgz)];
    cpos=ctar+Lx/2*[1 -1 .5];
end
set(ax,'xtick',p.ztick,'xlim',rgx);
set(ax,'ytick',p.ytick,'ylim',rgy);
set(ax,'ztick',p.xtick,'zlim',rgz);
cup=[sin(p.za)*[1 1]/sqrt(2) cos(p.za)];
set(ax,'CameraPosition',cpos);
set(ax,'CameraTarget',ctar);
set(ax,'CameraUpVector',cup);

if isfinite(p.cva)
  cva=p.cva;
  set(ax,'CameraViewAngle',cva);
else 
  cva=get(ax,'CameraViewAngle');
  rgx=rgx(:);rgy=rgy(:);rgz=rgz(:);zzz=[-p.bh;rgz(2)];
  p1=patch(rgx([1 1 2 2]),rgy([1 2 2 1]),zzz([1 1 1 1]),[0 0 0]);
  p2=patch(rgx([1 1 2 2]),rgy([1 2 2 1]),zzz([2 2 2 2]),[0 0 0]);
  p3=patch(rgx([1 1 1 1]),rgy([1 2 2 1]),zzz([1 1 2 2]),[0 0 0]);
  p4=patch(rgx([2 2 2 2]),rgy([1 2 2 1]),zzz([1 1 2 2]),[0 0 0]);
  p5=patch(rgx([1 1 2 2]),rgy([1 1 1 1]),zzz([1 2 2 1]),[0 0 0]);
  p6=patch(rgx([1 1 2 2]),rgy([2 2 2 2]),zzz([1 2 2 1]),[0 0 0]);
  rgx=rgx';rgy=rgy';rgz=rgz';zzz=zzz';
  F=[];
  cvr=1.2;
  nn=0;
  msk=repmat(0==1,p.sz);
  msk([1 2 end-1 end],:)=1==1;
  msk(:,[1:round((1-the(2))*szfp(2))+2 end-1 end])=1==1;
  msk=msk';
  for k=1:50
    disp(sprintf('%i cva:%8.4f cvr:%8.8f',k,cva,cvr));
    ff=gu_write_frame([],fh,szfp,dpi,p.aa,[]);
    ff=~all(ff==255,3);
    F(k)=~any(any(ff & msk));
    if F(k)==1 % Make it larger
      cva=cva/cvr;
    else % Make it smaller
      cva=cva*cvr;
    end
    set(ax,'CameraViewAngle',cva);
    if k<2
      continue;
    end
    if F(end-1)~=F(end)
      nn=nn+1;
    end
    if nn>2
      cvr=(1+cvr)/2;
    end
    if cvr<1.001 & F(end)==1
      break;
    end
  end
  delete([p1 p2 p3 p4 p5 p6]);
  p.cva=cva;
  save(pnpp,'p');
end

upq0=atan2(cup(2),cup(1));upf0=acos(cup(3));
cpr=sqrt(sum(cpos.^2));
cpq0=atan2(cpos(2),cpos(1));cpf0=acos(cpos(3)/cpr);


if p.bottom
  arg(5)=-p.bh;
  axis(arg);
  [hp,lp]=ded_3d_bottom(p,rgx,rgy);
end

if p.top
  p.lw=1.01;
  lc=[0.8 0.8 1.0]/2;
  [X Y]=ndgrid(-b.H/2:b.H/2,-b.W/2:b.W/2);
  Z=repmat(-0.5,size(X));
  line([X  X ],[Y   Y],0*[Z   Z],'color',lc,'linewidth',p.lw);
  line([X' X'],[Y' Y'],0*[Z' Z'],'color',lc,'linewidth',p.lw);
  [Z X Y]=ndgrid([0 -0.5],-b.H/2:b.H/2,[-b.W/2 b.W/2]);nsz=[2 length(Z(:))/2];
  line(reshape(X,nsz), reshape(Y,nsz), reshape(Z,nsz),'color',lc,'linewidth',p.lw);
  [Z X Y]=ndgrid([0 -0.5],-b.H/2*[-1 1],-b.W/2: b.W/2);
  nsz=[size(X,1) size(X,2)*size(X,3)];
  line(reshape(X,nsz), reshape(Y,nsz), reshape(Z,nsz),'color',lc,'linewidth',p.lw);
end

set(ax,'visible','on');
camlight(80,40);
camlight('left');
lightangle(-45,30);
set(ax,'ambientlightcolor',[1 1 1]);

nmb={nmb{1:min(length(nmb),p.maxfn)}};
j=0;
maxx=0;
maxy=0;
maxz=0;
hb=[];hc=[];

if p.box
  X=repmat(rgx([1 1 2 2]),2,1);
  Y=repmat(rgy([1 2 2 1]),2,1);
  Z=repmat([-p.bh;rgz(2)],1,4);
  hh=line(X,Y,Z,'color',[0 0 0],'linewidth',p.lw);
  
  X=rgx([1 1 2 2 1])';
  Y=rgy([1 2 2 1 1])';
  Z=repmat(rgz(2),5,1);
  hh=line(X,Y,Z,'color',[0 0 0],'linewidth',p.lw);
end

%set(gcf,'position',[1 -400 512        1024])

%set(ax,'BoxStyle','full','Box','on');
%set(ax,'PlotBoxAspectRatio',get(ax,'PlotBoxAspectRatio'));

fnpng=[];

if p.cont
  fh2=figure(2);
  clf
  figure(fh);
end
if p.transparent
  set(fh,'color','none');
end

x=x(:);
y=y(:);
z=z(:);

for kkk=1:p.reduce
  x=imhalf(x,1);
  y=imhalf(y,1);
  z=imhalf(z,1);
end

fno=[]; 
ts=ded_interp_stats(ss,b,t,p.tsn);
xx=x;
for j=1:length(t)
  fn=nmb{j};
  if p.png
    fnpng=sprintf('%s/%05u.png',dd,j);
  end
  if p.cont 
    if file_nt(fnpng,fn)
      disp(sprintf('loading: %s',fnpng))
      ff=imread(fnpng);
      figure(fh2);
      imshow(ff,'InitialMagnification','fit','Border','tight');
      drawnow;
      figure(fh);
      if fp>0
        yuv_write(fp,ff,'4:2:0','xy');%4:2:0'
      end
      continue;
    else
      delete(fh2);
      p.cont=0;
    end
  end
  

  if ~strcmp(fn,fno)
    i=1;
    try
      a=ded_read_hdf(fn);
      fno=fn;
    catch
      disp(sprintf('ded_3d: ded_read_hdf(''%s'') failed',fn));
      break;
    end
  else
    i=i+1;
  end
  f=a.(bnm)(:,:,:,i);
  for kkk=1:p.reduce
    f=imhalf(f,2);
    f=imhalf(f,3);
  end
  if ~isempty(zz)
    if size(f,1)==length(c.Jz)
      f=interp1(c.Jz,f,zz,'linear','extrap');
    elseif size(f,1)==length(c.Az)
      f=interp1(c.Az,f,zz,'linear','extrap');
    end
  end
  for kkk=1:p.reduce
    f=imhalf(f,1);
  end
  if isfield(p,'copyx')  
    if p.copyx
      f=cat(ndims(f),f,f);
    end
  end
  f=permute(f,p.pcx);
  f=max(p.rg(1),min(p.rg(2),f));
  if p.rg(2)==0 & p.rg(1)<0
    f=-f;
  end
  if ry
    f=repmat(permute(f,[2 1 3]),[2 1 1]);
  end
  if isfield(p,'flipy')  
    if p.flipy
      f=cat(1,f(end:-1:1,:,:),f);
    end
  end
  if isfield(p,'flipz')  
    if p.flipz
      f=cat(2,f(:,end:-1:1,:),f);
    end
  end
  if YFR
    if c.dim==2
      f=cat(1,f,f(1,:));
    else
      f=cat(1,f,f(1,:,:)); 
    end
  end
  if prod(size(f))~=prod([length(y) length(x) length(z)])
    disp(sprintf('ded_3d: %s Non-matching sizes [%u %u %u] [%u %u %u]',nm,size(f,1),size(f,2),size(f,3),length(y),length(x),length(z))); 
    return;
  end
  f=reshape(f,[length(y) length(x) length(z)]);
    
  qqq=pi/180*p.dq*t(j);
  cpos=cpr*[sin(cpf0)*[cos(cpq0+qqq) sin(cpq0+qqq)] cos(cpf0)];
  cup=[sin(p.za)*[cos(cpq0+qqq) sin(cpq0+qqq)] cos(p.za)];
  set(ax,'CameraPosition',cpos,'CameraUpVector',cup);
  if ~isempty(hb);if ishandle(hb);delete(hb);end;end;
  if ~isempty(hc);if ishandle(hc);delete(hc);end;end;
  ff=p.isoval; %*max(f(:));     %f(:,:,z<0)=0;
  if isfield(p,'X')
    x=xx-p.X(a.t);
  end
  if p.trackfront
    XX=interp1(ss.t,ss.X,double(a.t),'linear','extrap');
    XX=mod(XX-b.L/2,b.L)+b.L/2;
    x=xx-XX;
    if p.bottom
      delete(hp);
      delete(lp);
      [hp,lp]=ded_3d_bottom(p,rgx+[0 1]-mod(XX,1),rgy);
    end
  end

  for i=1:length(ff)
    try
      hb(i) = patch(isosurface(x,y,z,f,ff(i)),'FaceColor',p.col{i},'EdgeColor','none','clipping','on');
      isonormals(x(:),y(:),z(:),f,hb(i));
    end
  end
  if ~isempty(p.isocap)
    hc=patch(isocaps(x,y,z,f,ff,p.isocap),'FaceColor','interp','EdgeColor','none');
    set(hc,'ambientstrength',0.9);
  end
  if ~isempty(p.text)
    set(th,'string',ts{j});
  end
  if fp==0
    drawnow;
    pause(.1);
  else
    ff=gu_write_frame(fp,fh,szfp,dpi,p.aa,fnpng);
  end
end

if fp>0
  fclose(fp);
  if length(t)>1
    fps=max(5,min(50,round(5/median(diff(t)))));
  else
    fps=1;
  end
  make_mp4(fnyuv,fnmp4,szfp,fps,[],fnffmpeg);%,'yuv444p');
end

return;


function [hp,lp]=ded_3d_bottom(p,rgx,rgy)

[X,Y,Z]=ndgrid(rgx,rgy,[-p.bh 0]);
V=[X(:) Y(:) Z(:)];
F=[[1 2 4 3];[1 2 6 5];[1 3 7 5];[2 4 8 6];[5 6 8 7];[3 4 8 7];];
C=[linspace(0.8,0.9,8)' linspace(0.8,0.9,8)' linspace(0.9,1,8)'];
hp=patch('Faces',F,'Vertices',V,'FaceVertexCData',C);

lc=[0.8 0.8 1.0]/2;
gx=linspace(rgx(1),rgx(2),1+round(diff(rgx)/p.gdx));
gy=linspace(rgy(1),rgy(2),1+round(diff(rgy)/p.gdy));
[X Y]=ndgrid(gx,gy);
Z=repmat(-p.bh,size(X,1),size(X,2));

lp=line([X X],[Y Y], [0*Z Z],'color',lc,'linewidth',p.lw);
lp=[lp;line([X' X'],[Y' Y'], [0*Z' Z'],'color',lc,'linewidth',p.lw)];
[Z X Y]=ndgrid([0 -p.bh],gx,gy);nsz=[2 length(Z(:))/2];
lp=[lp;line(reshape(X,nsz), reshape(Y,nsz), reshape(Z,nsz),'color',lc,'linewidth',p.lw)];
[Z X Y]=ndgrid([0 -p.bh],gx,gy);
lp=[lp;line(reshape(X,nsz), reshape(Y,nsz), reshape(Z,nsz),'color',lc,'linewidth',p.lw)];

hp.EdgeColor = 'none';
hp.FaceColor = 'interp';
%hp.EdgeLighting = 'gouraud';
hp.FaceLighting = 'gouraud';
hp.AmbientStrength = 0.3;
hp.DiffuseStrength = 0.8;
hp.SpecularStrength = 0.9;
hp.SpecularExponent = 25;
hp.BackFaceLighting = 'unlit';
return

% $$$       axis('tight');
% $$$       aa=axis;
% $$$       maxx=max(maxx,max(abs(aa(1:2))));
% $$$       maxy=max(maxy,max(abs(aa(3:4))));
% $$$       maxz=max(maxz,aa(6));
% $$$       %      axis([maxx*[-1 1] maxy*[-1 1]  0 maxz]);

clf
b.H=8;
b.W=8;


hold('on');
hold('on');
view([40 20]);

C=[linspace(0.4,0.8,4)' linspace(0.4,0.8,4)' linspace(0.8,1,4)'];
h=0.05;
X=b.H/2*[-1 1 1 -1]';
Y=b.W/2*[1 1 -1 -1]';
Z=h*[1 1 1 1]';

hp=patch('xdata',[X X],'ydata',[Y Y],'zdata',[0*Z -Z],'FaceVertexCData',col);
%hp.EdgeColor = [0 0 0];
%hp.EdgeLighting = 'gouraud';
hp.FaceColor = 'interp';
hp.FaceLighting = 'gouraud';
hp.AmbientStrength = 0.3;
hp.DiffuseStrength = 0.8;
hp.SpecularStrength = 0.9;
hp.SpecularExponent = 25;
hp.BackFaceLighting = 'unlit';



lighting('gouraud');

daspect([1 1 1]);

set(ax,'ambientlightcolor',[1 1 1]);

nmb={nmb{1:min(length(nmb),p.maxfn)}};

clf;
h=patch('xdata',X,'ydata',Y,'zdata',Z,'FaceVertexCData',col);
view(0,75)
%shading interp
lightangle(-45,30)
h.EdgeColor = 'interp';
h.EdgeLighting = 'gouraud';
h.FaceColor    = 'interp';
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.3;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';


clf




%mpiexec -n 4 ded_gc.py --preset ccle --Nx 256 --Ny 1 -W None -T 5 --Re 50 --dtjavrg 1 --dtjb 1 --dtju 1 --dtjv 1 --dtjw 1 --dtjy 0.1 --dtforce 1 --dtstats 1 --dtxyz 0 --dtmomb 0 --dtavrg 0 --dtgc 0 --dtleft 1 --dtright 1 --dtslice 0 --slicex 0 --slicey 0 --slicez 0 --dtb 0 --dtu 0 --dtv 0 --dtw 0 --dtx 0 --dty 0 --dtz 0 --dtxy 0 --dtyz 0 --dtxz 0 --lbc slip --ubc slip gc/test/01
%mpiexec -n 4 ded_gc.py --preset ccle --Nx 256 --Ny 1 -W None -T 5
%--Re 50 --dtjavrg 1 --dtjb 0 --dtju 1 --dtjv 1 --dtjw 1 --dtjy 0.1
%--dtforce 1 --dtstats 1 --dtxyz 0 --dtmomb 0 --dtavrg 0 --dtgc 0
%--dtleft 1 --dtright 1 --dtslice 0 --slicex 0 --slicey 0 --slicez
%0 --dtb 1 --dtu 0 --dtv 0 --dtw 0 --dtx 0 --dty 0 --dtz 0 --dtxy 0
%--dtyz 0 --dtxz 0 --lbc slip --ubc slip gc/test/01




p=ded_3d_param('gc/test/09');
p.copyx=1;
p.rgx=[-8 8];
p.trackfront=true;
ded_3d('gc/test/09',p);
