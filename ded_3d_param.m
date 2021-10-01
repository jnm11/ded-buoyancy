function p=ded_3d_param(nm,p)
%Return a parameter structure for a simulation for 3d plotting

if nargin<2
  p=struct();
end

dd=['~/tmp/dedalus/' nm];
pnpp=[dd '.mat'];
if ~isdir(dd)
  mkdir(dd);
end

b=ded_read_param(nm);
if isempty(b);p=[];return;end

dim = (b.Nx>1)+(b.Ny>1)+(b.Nz>1);
  
if isempty(b)
  disp(sprintf('ded_3d_param: No parameter file found for %s',nm'));
  keyboard
  p=[];
  return
end
nmb=ded_get_fn(nm,'b');
if isempty(nmb)
  nmb=ded_get_fn(nm,'s');
end
if isempty(nmb)
  disp(sprintf('ded_3d: No density files for "%s"',nm));
  p=[];
  return;
end
if ~isfile(nmb{1})
  disp(sprintf('ded_3d: Not a file "%s"'),fns{1});
  p=[];
  return;
end
a=ded_read_hdf(nmb{1});
f=ded_read_g(nm,'force');
c=ded_coord(nm);
if isempty(c)
  p=[];;
  return;
end

if dim==3
  ry=0;
else
  c.Jy=[0 b.W];
  ry=1;
end

if isfield(b,'Ty'); pd.flipy=strcmp(b.Ty,'SinCos'); end;
if isfield(b,'Tz'); pd.flipz=strcmp(b.Tz,'SinCos'); end;
pd.aa=1;
pd.png=1;
pd.rgx=[0 b.L];
pd.rgy=[0 b.W] - (c.Jy(1)<0)*b.W/2;
pd.rgz=[0 b.H] - (c.Jz(1)<0)*b.H/2;
pd.cva=NaN;

switch(b.sType)
  case 'gc'
    pd.sz=[1024 576]; % 16/9
    pd.rgx=[0 b.L];
    if b.H<=2 & b.W<=2 & b.L<=48
      pd.cva=32;
    end
  case {'pm','plume'}
    b.sType='pm';
    sz=[sqrt(b.H.^2+b.W.^2) b.L];
    pd.sz = 32+16*round(sz/max(sz)*(512-32)/16); % 16/9
    if b.L==32 & b.H==8 & b.W==8
      pd.cva=2;
    end
    
end
[nm2 nm1]=fileparts(nm);
[nm3 nm3]=fileparts(nm2);      

switch(nm)
  case 'gc/f6/f/10'
    pd.reduce=2;
    pd.flipy=1;
end
switch(nm1)
end
switch(nm2)
end   
fnm2=max(find(nm2=='/'));

if ~isempty(fnm2)
  nm2b=nm2(1:fnm2-1);
  switch(nm2b)
    case 'gc/f7/i'
      pd.cva=34;
    case 'gc/L'
      pd.cva=38.5;
  end   
end

switch(nm3)
  case {'qgc2','qgc2b','qgc3','qgc3b','qgc4','qgc4b','qgc5','qgc5b','qgc6','qgc6b'}
    pd.cva=29.9;
  case {'qgc7a','qgc7b','qgc7b','qgc7c','qgc7d','qgc7e','qgc8a','qgc8b','qgc8b','qgc8c','qgc8d','qgc10a'}
    pd.cva=29.9;
  case {'gc2d7a','gc2d7b','gc2d7b','gc2d7c','gc2d7d','gc2d7e','gc2d7f','gc2d8a'}
    pd.cva=30;
    pd.zmax=0;
    pd.sz=[1024 480];
  case 'f63'
    pd.cva=30;
    pd.zmax=0;
    pd.sz=[1024 480];
  case{'ccle'}
    pd.cva=32.1156;
    pd.sz=[1024 576]; % 16/9
  case 'f3'
    pd.cva=30;
    pd.sz=[1024 480];
    pd.zmax=0;
  case 'qpmf4'
    pd.cva=2.6;
    pd.sz=[256 512];
    pd.zmax=0;
end   

switch([nm3 '/' nm1])
  case {'pm/005','pm/006','pm/007'}
    pd.rgx=[0 30];
    pd.cva=4.8;
    pd.zmax=0;
    pd.sz=[480 512];
  case {'pm/008'}
    pd.rgx=[0 40];
    pd.cva=2.0;
    pd.zmax=0;
    pd.sz=[256 576];
  case {'pm/020','pm/021','pm/022','pm/023'}
    pd.cva=5;
end

pd.gdx=0.5;
pd.gdy=0.5;
p=combine_struct(pd,p);


pd.box=1;
pd.bottom=1;
pd.png=0;
pd.isoval=0.01;
pd.isocap='all';  

pd.rg=[-inf inf];

fmaxx=-inf;
fmaxy=-inf;
fmaxz=-inf;
fminx=+inf;
fminy=+inf;
fminz=+inf;
X=NaN;
Y=NaN;
Z=NaN;

if ~isempty(f)
  nms=fieldnames(f);
  for j=1:length(nms)
    tol=1e-3;
    nm=nms{j};
    if nm(1)=='w' & all(nm~='_')
      F=f.(nm);
      F=F(:,:,:,1);
      F=F>max(F(:))*tol;
      if isfield(f,'y')
        X=f.x(squeeze(any(any(F,2),1)));fmaxx=max(fmaxx,max(X));fminx=min(fminx,min(X));
        Y=f.y(squeeze(any(any(F,3),1)));fmaxy=max(fmaxy,max(Y));fminy=min(fminy,min(Y));
        Z=f.z(squeeze(any(any(F,3),2)));fmaxz=max(fmaxz,max(Z));fminz=min(fminz,min(Z));
      else
        X=f.x(squeeze(any(F,1)));fmaxx=max(fmaxx,max(X));fminx=min(fminx,min(X));
        Z=f.z(squeeze(any(F,2)));fmaxz=max(fmaxz,max(Z));fminz=min(fminz,min(Z));
        fminy=0;fmaxy=b.W;
      end
    end
  end
end

pd.fminx=fminx;
pd.fmaxx=fmaxx;
pd.fminy=fminy;
pd.fmaxy=fmaxy;
pd.fminz=fminz;
pd.fmaxz=fmaxz;

switch b.sType
  case {'pm','plume'}
    pd.sz=[256 512];   % Output size
    pd.dq=10; % rotate velocity
    pd.top=0;
    pd.bottom=0;
    pd.x0=0;
    pd.y0 = (c.y(1)>=0)*b.W/2;
    pd.z0 = (c.z(1)>=0)*b.H/2;
    pd.top=1;
    pd.bh=0.5;
    pd.pcx=[2 1 3];
  case 'gc'
    pd.sz=[1024 512];   % Output size
    pd.dq=10; % rotate velocity
    pd.x0 = 0;
    pd.y0 = 0;
    pd.z0 = 0;
    pd.top=0;
    pd.bh=b.H/20;
    pd.dq=0;
    pd.pcx=[2 3 1];
 end

pd.incremental=0;
pd.col={'red','blue'};
pd.rgx=[2 0.2];
pd.maxfn=inf;
pd.nrot=1;
pd.zmax=0.7;
pd.trackfront=0;
pd.aa=2;
pd.rd='zbuffer';
pd.rd='opengl';
pd.rgx=[-inf inf];
pd.tsn=NaN;
pd.za=0;
pd.ctar=1;
pd.cont=1;


p=combine_struct(pd,p); % Overwrite defaults with parameters

if isfield(b,'Forcing')
  if b.Forcing==8 & strcmp(b.sType,'gc')% & 0 
    p.copyx=1;
    p.rgx=[-1 1]*b.L/2;
    p.trackfront=true;
  end
end

if p.cont
  p.png=1;
  if isfile(pnpp)
    pp=load(pnpp);
    if ~isfield(pp,'rg')
      pp.rg=[0 1];
    end
    rmfld={'time'};
    p=rmfieldifexist(p,rmfld);
    pp.p=rmfieldifexist(p,rmfld);
    if (struct_cmp(p,pp.p)~=0)
      disp(sprintf('ded_3d_param: %s continuation structures do not match. deleting pngs',nm));
      unix(sprintf('/bin/rm -rf ~/tmp/dedalus/%s*',nm'));
    else
      %disp(sprintf('ded_3d_param: %s continuation structures match',nm));
      p=pp.p;
    end
  end
end
save(pnpp,'p');


if ~isfield(p,'cva')
  keyboard
end