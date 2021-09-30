function ded_movie_slice(q)
%ded_gc_slice('gc/emle/013');ded_gc_slice('gc/emle/014');ded_gc_slice('gc/emle/015');ded_gc_slice('gc/emle/016');ded_gc_slice('gc/emle/017');ded_gc_slice('gc/emle/018');ded_gc_slice('gc/emle/019');ded_gc_slice('gc/ccle/024');ded_gc_slice('gc/ccle/046');

if nargin==0
  ded_movie_slice_all;
  return;
end

ttol=1000;
for i=1:length(q.typ)
  if strcmp(q.typ{i},'rho')
    fnB=ded_get_fn(q.nm,'b');
    fnS=ded_get_fn(q.nm,'s');
    tB=round(ttol*ded_get_times(fnB));
    tS=round(ttol*ded_get_times(fnS));
    t{i}=intersect(tB,tS);
  else
    fns{i}=ded_get_fn(q.nm,q.typ{i});
    t{i}=round(ttol*ded_get_times(fns{i}));
  end
  if i==1
    ut=t{i};
  else
    ut=intersect(ut,t{i});
  end
end
p=ded_read_param(q.nm);
c=ded_coord(q.nm);

if ~isfield(p,'rgx');  p.rgx = [0 p.L]; end
if ~isfield(p,'rgy');  p.rgy = [0 p.L]; end
if ~isfield(p,'rgz');  p.rgz = [0 p.L]; end


dd=['~/tmp/Dedalus/' q.nm '-slice'];
if ~isdir(dd)
  mkdir(dd);
end
fq.nmp4=['~/films/Dedalus/' q.nm '-slice.mp4'];
fnyuv=[dd '.yuv'];
fnffmpeg=[dd '.ffmpeg'];
fp=fopen(fnyuv,'w');


wx=c.dJx;
wy=c.dJy;
wz=c.dJz;

if strcmp(p.Tz,'Chebychev'); wz=ichebIw(p.NJz); end


ntype=length(q.stype);
for i=1:ntype
% $$$   r1{i}=find(c.Jz>=q.rgz(1) & c.Jz<=q.rgz(2));
% $$$   r2{i}=find(c.Jy>=q.rgy(1) & c.Jy<=q.rgy(2));
% $$$   r3{i}=find(c.Jx>=q.rgx(1) & c.Jx<=q.rgx(2));
  r1{i}=1:length(c.Jz);
  r2{i}=1:length(c.Jy);
  r3{i}=1:length(c.Jx);
  switch(q.stype{i})
    case {'x','Ix'}
      r3{i}=findmin(abs(c.Jx-q.sloc(i)));
      Iw=wx;
      c1{i}=c.Jy(r2{i});
      c2{i}=c.Jz(r1{i});
      xlim{i}=q.rgy;
      ylim{i}=q.rgz;
    case {'y','Iy'}
      r2{i}=findmin(abs(c.Jy-q.sloc(i)));
      Iw{i}=wy;
      c1{i}=c.Jx(r3{i});
      c2{i}=c.Jz(r1{i});
      xlim{i}=q.rgx;
      ylim{i}=q.rgz;
    case {'z','Iz'}
      r1{i}=findmin(abs(c.Jz-q.sloc(i)));
      sc{i}=c.Jz-q.sloc(i);
      Iw{i}=wz;
      c1{i}=c.Jx(r3{i});
      c2{i}=c.Jy(r2{i});
      xlim{i}=q.rgx;
      ylim{i}=q.rgy;
  end
  if q.flip(i);
    c3=c2{i};
    c2{i}=c1{i};
    c1{i}=c3;
    tlim=xlim{i};
    xlim{i}=ylim{i};
    ylim{i}=tlim;
  end
end

%s=ded_read_stats(q.nm);
clf;
[szfp szfi szpi dpi]=gu_film_setup(q.sz,q.aa,q.rd);
fh=gcf;
ah=jsubplot(q.jsp,[0.05 0.05],[0.02  0.02],[0.02  0.02])';
for i=1:length(c1)
  axes(ah(i));
  ih(i)=pcolor(c1{i},c2{i},zeros(length(c2{i}),length(c1{i})));
  shading('interp');
  set(ah(i),'dataaspectratio',[1 1 1],'xlim',xlim{i},'ylim',ylim{i},'ydir','normal','xtick',[],'ytick',[]);
  th(i)=text(mean(xlim{i}),ylim{i}(2),'','verticalalignment','bottom','horizontalalignment','center');
  su=q.zerof(i);
  n1=max(0,round(-1000*(q.minf(i)+su)));
  n2=max(0,round( 1000*(su+su)));
  n3=max(0,round( 1000*(q.maxf(i)-su)));
  n1=hot(n1);
  n3=parula(n3);
  s3=repmat(max(0,linspace(-1, 1,n2))',1,3);
  s1=repmat(max(0,linspace( 1,-1,n2))',1,3);
  if ~isempty(n1)
    n2=s3.*n3(1,:)+s1.*n1(end,:)+1-s1-s3;
  else
    n2=s3.*n3(1,:)+1-s3;
  end
  set(ah(i),'colormap',[n1;n2;n3],'clim',[q.minf(i),q.maxf(i)]);
  line(xlim{i}([1 2 2 1 1]),ylim{i}([1 1 2 2 1]),'color',[0 0 0],'linewidth',2);
end
set(ah,'box','off');

if ~isfield(q,'start');q.start=1;      end;
if ~isfield(q,'end');  q.end=length(ut);end;
for j=1:length(ut)%q.start:q.end
  fnpng=sprintf('%s/%06u.png',dd,j);
  if isfile(fnpng)
    continue
  end
  for i=1:length(q.stype)
    jj=find(ut(j)==t{i});
    if isempty(jj)
      continue;
    end
    if strcmp(q.typ{i},'rho')
      jb=find(ut(j)==tB);
      js=find(ut(j)==tS);
      ab=ded_read_hdf(fnB{jb});
      as=ded_read_hdf(fnS{js});
      if isempty(ab) | isempty(as);break;end
      f=ab.b*p.B+as.s*p.S;
      ct=(ab.t+as.t)/2;
      clear('ab');
      clear('as');
    else
      a=ded_read_hdf(fns{i}{jj});
      if isempty(a);break;end
      f=a.(q.typ{i});
      ct=a.t;
      clear('a');
    end
    axes(ah(i));
    f=squeeze(f(r1{i},r2{i},r3{i}));
    set(th(i),'string',sprintf('%s %s t:%6.2f',q.typ{i},q.stype{i},ct));
    if q.flip(i);f=f';end
    if strcmp(q.P1,'odd' );f=[-flip(f,1);f];end
    if strcmp(q.P1,'even');f=[ flip(f,1);f];end
    if strcmp(q.P2,'odd' );f=[-flip(f,2),f];end
    if strcmp(q.P2,'even');f=[ flip(f,2),f];end
    if isfinite(q.scalef(i));f=q.scalef(i)*f;end;
    if isfinite(q.minf(  i));f=max(q.minf(i),f);end;
    if isfinite(q.maxf(  i));f=min(q.maxf(i),f);end;
    %set(th(i),'string',sprintf('%s: t=%6.2f, X=%6.2f',q.nm,a.t,X));
    set(ih(i),'cdata',f,'xdata',c1{i},'ydata',c2{i});
    clear('f');
  end
  if fp>0
    gu_write_frame(fp,fh,szfp,dpi,q.aa,fnpng);
  end
end

fps=15;
make_mp4(fnyuv,fq.nmp4,szfp,fps,[],fnffmpeg);%,'yuv444p');



return;
function ded_movie_slice_all
q.sz=[1280 1024];
q.aa=2;
q.rd='zbuffer';
q.stype={'y','z','y','z','y','z','y','z'};
q.sloc=[0 0 0 0 0 0 0 0];
q.P1='';
q.P2='';
q.flip=[1 1 1 1 1 1 1 1];
q.jsp=[4 2];
q.typ={'u','u','s','s','b','b','rho','rho'};
q.rgx=4+[0 15];
q.rgy=5*[-1 1];
q.rgz=5*[-1 1];
nms={'pm/f7/e/27','pm/f7/e/28','pm/f7/e/26','pm/f7/e/25','pm/f7/e/24'};
%q.start=197;
%q.end=199;
for j=1:length(nms)
  q.nm=nms{j};
  p=ded_read_param(q.nm);
  rhor=p.B-p.S;
  q.minf=   [-0.1 -0.1  0    0    0    0    -0.50 -0.50];
  q.maxf=   [ 1.0  1.0  1    1    1    1     1.00  1.00];
  q.zerof=  [ 0.05 0.05 0.1  0.1  0.1  0.1   0.20  0.20];
  q.scalef= [ 1.0  1.0  5    5    5    5    30   30];
  ded_movie_slice(q)
end
return;








  
q=rmfield(q,'minf');
q=rmfield(q,'maxf');
q=rmfield(q,'scalef');
q.typ='u';ded_movie_slice(q)



  k=findmin(abs(a.z-p.H/2));
  if size(a.b,1)==p.Nz
    bI=squeeze(sum(w.*a.b,1));
    b1=squeeze(a.b(k,:,:));
  elseif size(a.b,2)==p.Nz
    bI=squeeze(sum(w.*a.b,2));
    b1=squeeze(a.b(:,k,:,:));
  elseif size(a.b,3)==p.Nz
    bI=squeeze(sum(w.*a.b,2));
    b1=squeeze(a.b(:,:,k));
  end
  fy=find(size(a.b)==p.Ny);
  b2=squeeze(mean(a.b,fy));
  fz=find(size(a.b)==p.Nz);
  b2=squeeze(mean(a.b,fy));
  b2=ichebc2f(ichebf2c(b2,fz),fz,2*z/p.H-1);
  
  if size(b1,1)==p.Nx
    b1=b1';
    bI=sum(bI,2);
  else
    bI=sum(bI,1);
  end
  bI=bI-0.05*max(bI(:));
  f=max(find(bI>0));
  X=(a.x(f)*bI(f+1)-a.x(f+1)*bI(f))/(bI(f+1)-bI(f));
  %X=interp1(s.t,s.X,max(s.t(1),a.t));
  ar1=2*p.W/(p.H+2*p.W)*0.9;
  ar2=  p.H/(p.H+2*p.W)*0.9;
  if j==1;
    y=[-flip(a.y);a.y];
    ah1=axes('position',[0.1 0.05+ar2 0.89 ar1]);
    ih1=imagesc(a.x-X,y,[b1;flip(b1,1)]);
    set(ah1,'dataaspectratio',[1 1 1]);
    set(ah1,'ylim',[-p.W p.W]);
    set(ah1,'xlim',[-4.8 0.2]);
    set(ah1,'ydir','normal');
    th=text(-1.5,p.W,'','verticalaligq.nment','top','horizontalaligq.nment','center');
    ah2=axes('position',[0.1 0.05 0.89 ar2]);
    ih2=imagesc(a.x-X,z,b2);
    set(ah2,'dataaspectratio',[1 1 1]);
    set(ah2,'ylim',[0 p.H]);
    set(ah2,'xlim',[-4.8 0.2]);
    set(ah2,'ydir','normal');
  end
