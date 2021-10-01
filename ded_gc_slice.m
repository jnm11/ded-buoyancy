function ded_gc_slice(nm)
%ded_gc_slice('gc/emle/013');ded_gc_slice('gc/emle/014');ded_gc_slice('gc/emle/015');ded_gc_slice('gc/emle/016');ded_gc_slice('gc/emle/017');ded_gc_slice('gc/emle/018');ded_gc_slice('gc/emle/019');ded_gc_slice('gc/ccle/024');ded_gc_slice('gc/ccle/046');

p=ded_read_param(nm);

fns=ded_get_fn(nm,'b');

pd.rd='zbuffer';
pd.aa=1;
pd.sz=[1024 512];


dd=['~/tmp/dedalus/' nm '-slice'];
if ~isdir(dd)
  mkdir(dd);
end
fnmp4=['~/films/dedalus/' nm '-slice.mp4'];
fnyuv=[dd '.yuv'];
fnffmpeg=[dd '.ffmpeg'];
fp=fopen(fnyuv,'w');

clf;
[szfp szfi szpi dpi]=gu_film_setup(pd.sz,pd.aa,pd.rd);
fh=gcf;

w=ichebIw(p.Nz);

s=ded_read_stats(nm);
z=linspace(0,p.H,2*p.Nz);
for j=1:length(fns)
  fnpng=sprintf('%s/%05u.png',dd,j);
  a=ded_read_hdf(fns{j});
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
    th=text(-1.5,p.W,'','verticalalignment','top','horizontalalignment','center');
    ah2=axes('position',[0.1 0.05 0.89 ar2]);
    ih2=imagesc(a.x-X,z,b2);
    set(ah2,'dataaspectratio',[1 1 1]);
    set(ah2,'ylim',[0 p.H]);
    set(ah2,'xlim',[-4.8 0.2]);
    set(ah2,'ydir','normal');
  end
  set(th,'string',sprintf('%s: t=%6.2f, X=%6.2f',nm,a.t,X));
  set(ih1,'cdata',[b1;flip(b1,1)]);
  set(ih1,'xdata',a.x-X);
  set(ih2,'cdata',b2);
  set(ih2,'xdata',a.x-X);
  if fp>0
    gu_write_frame(fp,fh,szfp,dpi,pd.aa,fnpng);
  end
end

fps=15;
make_mp4(fnyuv,fnmp4,szfp,fps,[],fnffmpeg);%,'yuv444p');
