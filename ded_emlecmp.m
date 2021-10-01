function [a nms nmb times]=ded_emlecmp(tt,oneplot,typ,nms,nmb,times)

if nargin==0
  %[a nms nmb times]=ded_emlecmp(0.5,0,{'b','bb','vv'},{'ccle/025'});
  %save('~/025-time.mat','nms','nmb','times');
  nms=[];nmb=[];times=[];
  nms={'emle/019','emle/017','ccle/024','ccle/046'};nmb=[];times=[];
  for tt=8.5:0.5:25.5;
    [a nms nmb times]=ded_emlecmp(tt,0,{'b','bb','vv'},nms,nmb,times);
  end
  return;
end
%ded_emlecmp( 0.5,0,{'b','bb','vv'},{'ccle/024'});
%ded_emlecmp(12.5,0,{'b','bb','vv'},{'emle/017'});
if nargin<1
  tt=[];
end
if nargin<2
  oneplot=0;
end
if nargin<3
  typ=[];
end
if nargin<4
  nms=[];
end
if nargin<5
  nmb=[];
end
if nargin<6
  times=[];
end

if isempty(typ)
  typ={'b','bb','uu','vv','ww','uw'};
  typ={'b','bb','vv'};
end
if ~iscell(typ)
  typ={typ};
end

if isempty(nms)
  nms={'ccle/046','emle/001','emle/002','emle/003','emle/004','emle/005','emle/006','emle/007'};
  nms={'emle/010','emle/009','emle/008','emle/004','emle/005','emle/006','emle/007','emle/011','emle/012','emle/013','emle/014','emle/015','emle/016'};
  nms={'ccle/046','emle/001','emle/002','emle/003','emle/004',...
       'emle/005','emle/006','emle/007','emle/008','emle/009',...
       'emle/010','emle/011','emle/012','emle/013','emle/014',...
       'emle/015','emle/016','emle/017'};
  nms={'emle/001','emle/002','emle/003','emle/004','emle/005',...
       'emle/006','emle/007','emle/008','emle/009','emle/010',...
       'emle/011','emle/012','emle/013','emle/014','emle/015',...
       'emle/016','emle/017',...
       'ccle/016','ccle/017','ccle/019','ccle/023','ccle/024','ccle/025',...
       'ccle/026','ccle/027','ccle/028','ccle/029','ccle/030','ccle/031',...
       'ccle/046','ccle/082'};
  fn1=cellstr_ls('~/gc/emle/*');
  fn2=cellstr_ls('~/gc/ccle/*');
  fns={fn1{:},fn2{:}};
  for j=1:length(fns)
    [dd nm1]=fileparts(fns{j});
    [dd nm2]=fileparts(dd);
    dd=['~/Dropbox/Jim-Claudia-GC/mat/' nm2 '/' nm1];
    if ~isdir(dd); mkdir(dd); end
    nms{j}=[nm2 '/' nm1];
  end
end

n=length(nms);

quiet=1;
if isempty(times) | isempty(nmb)
  jj=0;
  p=[];
  nnms={};
  nmb={};
  times={};
  for j=1:n
    [ttimes nnmb]=ded_get_times_typ(['gc/' nms{j}],'y');
    nmbb=ded_get_fn(['gc/' nms{j}],'y');
    if isempty(nmbb)
      disp(sprintf('ded_emlecmp: no y files %s',nms{j}));
      continue;
    end
    if isempty(ttimes)
      disp(sprintf('ded_emlecmp: no times in y files %s',nms{j}));
      continue;
    end
    jj=jj+1;
    nnms{jj}=nms{j};
    times{jj}=ttimes; 
    nmb{jj}=nnmb;
  end
  nms=nnms;
end

p=ded_read_param(cellstrprefix('gc/',nms));

ttol=0.05;

tol=1e-7;

n=length(times);
DHOME=ded_dedalus_data_dir;
e=zeros(n,1);
eh=zeros(n,1);
hfn={};
for j=1:n
  [e(j), m]=min(abs(tt-times{j}));
  hfn{j}=[nmb{j}{m}];
  if isfile(hfn{j})
    eh(j)=1;
  elseif isfile([DHOME '/' hfn{j}])
    hfn{j}=[DHOME '/' hfn{j}];
    eh(j)=1;
  else
    eh(j)=0;
  end
end

% $$$ f=find(e<=ttol & eh);
% $$$ if isempty(f)
% $$$   a=[];
% $$$   disp(['ded_emlecmp: e<ttol & eh' sprintf(' %8.6f',ttol)])
% $$$   keyboard
% $$$   return;
% $$$ end
f=1:length(hfn);


hfn={hfn{f}};
nnms={nms{f}};
W=[p(f).W];
Wx=[p(f).Wx];
dx=[p(f).L]./[p(f).Nx];
Wn=Wx./dx;
n=length(hfn);
a=[];
for j=1:n
  [dd fn]=fileparts(hfn{j});
  fnmat=[dd '/' fn '.mat'];
  if ~isfile(fnmat) 
    aa=ded_read_hdf(hfn{j});
    if isempty(aa)
      disp(sprintf('ded_emlecmp: %s is empty',hfn{j}));
      continue;
    end
    if isempty(aa.b)
      disp(sprintf('ded_emlecmp: %s is empty',hfn{j}));
      continue;
    end
    if isfield(aa,'sim_time')
      fff=findmin(abs(aa.sim_time-tt));
      t1(j)=aa.sim_time(fff);
      fnms=fieldnames(aa);
      for k=1:length(fnms)
        if ndims(aa.(fnms{k}))==4
          aa.(fnms{k})=squeeze(aa.(fnms{k})(:,:,:,fff));
        end
      end
    else
      t1(j)=aa.t1;
    end
    aa.t1=t1(j);
    aa=ded_zgrid(aa,400,{},[],[],[],1,1);
    aa.nms=nnms{j};
    aa.t1=t1(j);
    save(fnmat,'aa');
    disp(sprintf('Saving %s',fnmat'));
  else
    disp(sprintf('Loading %s',fnmat'));
    load(fnmat);
  end
  aa.nms=nnms{j};
  a=struct_array_append(a,aa,[],quiet);
end

t=[a.t1];
% $$$ for j=1:n
% $$$   
% $$$ end

for j=1:n
  [res b]=unix(sprintf('grep -- --dimple ~/gc/%s/*out',a(j).nms));
  if res==0
    f=findstr('--dimple',b);
    dimple(j)=str2num(strtok(b(f+8:end)));
  else
    dimple(j)=0;
  end
  W(j)=median(a(j).b(a(j).b>mean(a(j).b(:)))); % instead of W(j);
  xb(j)=0.1+solve_first(a(j).x,max(a(j).b)/W(j)-0.05);
end

set(0,'defaulttextfontsize',8)
X=4;
fd='~/Dropbox/Jim-Claudia-GC/ofigs';


ff=zeros(n,1);
for k=1:length(typ)
  for j=1:n
    ff(j)=~isfield(a(j),typ{k});
  end
  if all(ff) ; continue; end
  
  if oneplot
    figure(k);clf;
    
    nsh=floor(sqrt(n));
    nsh=[nsh ceil(n/nsh)];
    preprint(nsh.*[X+0.1 1.1]*0.75+[0.1 0.1],6);
    ah=jsubplot(nsh,[0.01 0.05],[0.01 0.05],[0.01 0.05]);
  else
    for j=1:n
      figure(j);
      clf;
      preprint([X 1]+11/72*[0 2],7);
      ah(j)=axes;
    end
  end
  maxf=-inf;
  minf=0;
  for j=1:n
    if ff(j); continue; end
    switch(typ{k})
      case 'b'
        f=a(j).b/W(j);
        maxf=1;
      case 'bb'
        f=max(0,a(j).bb/W(j)-(a(j).b/W(j)).^2);
      case 'uu'
        f=max(0,a(j).uu/W(j)-(a(j).u/W(j)).^2);
      case 'vv'
        f=max(0,a(j).vv/W(j));
      case 'ww'
        f=max(0,a(j).ww/W(j)-(a(j).w/W(j)).^2);
      case 'uw'
        f=a(j).uw/W(j)-a(j).u.*a(j).w/W(j)^2;
        minf=min(minf,min(f(:)));
    end
    maxf=max(maxf,max(f(:)));
    axes(ah(j));
    imagesc(a(j).x,a(j).z,f);
    set(ah(j),'xlim',[-X 0]+xb(j),'ylim',[0 1]);    
    set(ah(j),'xtick',(ceil(2*(xb(j)-X)):floor(2*xb(j)))/2);
    title(sprintf('%s %s, t=%5.2f, Wx=%7.4f %4.1f, dimp=%5.3f, Sc=%3.1f',typ{k},a(j).nms,t(j),Wx(j),Wn(j),dimple(j),p(j).Scb))
  end
  set(ah,'dataaspectratio',[1 1 1])
  set(ah,'clim',[minf max(minf+1e-5,maxf)],'ydir','normal');
  set(ah,'yticklabels',[],'ytick',[]);
  if oneplot
    set(ah(:,1:end-1),'xticklabels',[],'xtick',[-X 0 X]);
    fn=sprintf('%s/emle-%s-%04u.eps',fd,typ{k},tt*100);
    disp(['writing figure ' fn]);
    print(k,'-depsc2',fn);
  else
    pp=zeros(n,4);
    for j=1:n
      make_tight_axes;
      pp(j,:)=get(gca,'position');
    end
    pp=[max(pp(:,1:2),[],1) min(pp(:,3:4),[],1)]; 
    set(ah,'position',pp);
    for j=1:n
      if ff(j); continue ; end
      fn=sprintf('%s/%s/%s-%04u.eps',fd,a(j).nms,typ{k},100*tt);
      dd=fileparts(fn);
      if ~isdir(dd);mkdir(dd);end;
      disp(['writing figure ' fn]);
      print(j,'-depsc2',fn);
    end
  end  
end


return;

clear x b 
for j=1:n
  bb=ded_read_hdf(sprintf('~/gc/%s/b/b-00000.hdf5',a(j).nms));
  b{j}=squeeze(mean(mean(bb.b,1),2));
  x{j}=bb.x;
end

figure(6);clf;
hold('on');
for j=1:n
  plot(x{j},b{j},'s-');
end
axis([-0.05 0.05 0 1]);


ht=findobj('FontName','Helvetica');
set(ht,'fontsize',6,'fontweight','normal');

return;

if 0
  for j=1:n
    hfn{j}=sprintf('~/gc/%s/y/y-%05u.hdf5',nms{j},m);
    if isfile(hfn{j})
      disp(['Exists: ' hfn{j}]);
    else
      if ~strcmp(ded_dedalus_data_dir('HOSTNAME'),'hamilton')
        disp(['Copying: ' hfn{j}]);
        unix(sprintf('rsync -vap --progress tpos:gc/%s/y/y-%05u.hdf5 ~/gc/%s/y 2> /dev/null',nms{j},m,nms{j}));
      end
    end
    if ~isfile(hfn{j})
      disp(['Failed to copy: ' hfn{j}]);
    end
    pfn=sprintf('~/gc/%s/param.h5',nms{j});
    if ~strcmp(ded_dedalus_data_dir('HOSTNAME'),'hamilton')
      if ~isfile(pfn)
        unix(sprintf('rsync -vap tpos:gc/%s/param.h5 ~/gc/%s',nms{j},nms{j}));
      end
    end
  end
end
