function ded_gc_ccle_cplots(fns,T,pfn)

if nargin<1
  fns=[];
end
if nargin<2
  T=[];
end
if nargin<3
  pfn='';
end

% $$$ fns=sort({'mat-ccle-017/017-time.mat','mat-ccle-022/022-time.mat','mat-ccle-046/046-time.mat',...
% $$$      'mat-ccle-020/020-time.mat','mat-ccle-024/024-time.mat','mat-ccle-021/021-time.mat',...
% $$$      'mat-ccle-025/025-time.mat','mat-emle-001/001-time.mat','mat-emle-002/002-time.mat','mat-emle-003/003-time.mat'});
% $$$ fns={'mat-ccle-046/046-time.mat','mat-emle-001/001-time.mat','mat-emle-002/002-time.mat','mat-emle-003/003-time.mat',...
% $$$      'mat-emle-004/004-time.mat','mat-emle-005/005-time.mat','mat-emle-006/006-time.mat','mat-emle-007/008-time.mat',...
% $$$      'mat-emle-009/009-time.mat','mat-emle-010/010-time.mat','mat-emle-011/011-time.mat','mat-emle-012/012-time.mat',...
% $$$      'mat-emle-013/013-time.mat','mat-emle-014/014-time.mat','mat-emle-015/015-time.mat','mat-emle-016/016-time.mat'};
% $$$ fns={'mat-emle-004/004-time.mat','mat-emle-005/005-time.mat','mat-emle-006/006-time.mat','mat-emle-007/007-time.mat',...
% $$$      'mat-emle-008/008-time.mat','mat-emle-010/010-time.mat','mat-emle-009/009-time.mat',...
% $$$      'mat-emle-011/011-time.mat','mat-emle-012/012-time.mat','mat-emle-013/013-time.mat','mat-emle-014/014-time.mat',...
% $$$      'mat-emle-015/015-time.mat','mat-emle-016/016-time.mat','mat-emle-017/017-time.mat'};
%ded_gc_ccle_cplots([],6,'all-6-')
%ded_gc_ccle_cplots({'emle/017','ccle/024','ccle/046'},33.5,'017-33-')
if isempty(fns)
  fn1=cellstr_ls('~/gc/emle/*');
  fn2=cellstr_ls('~/gc/ccle/*');
  fns={fn1{:},fn2{:}};

  fns={'ccle/024','emle/004','emle/005','emle/006','emle/007',...
       'emle/008','emle/010','emle/009',...
       'emle/011','emle/012','emle/013','emle/014',...
       'emle/015','emle/016','emle/017','emle/018','emle/019'};
end
if isempty(T)
  T=6;
end
trg=[0 T];


dfig=sprintf('~/Dropbox/Jim-Claudia-GC/ofigs/%s',pfn);

sz=[6 3];
fnt=10;

k=0;
n=length(fns);
a=[];
quiet=1;
for j=1:n
  [dd nm1]=fileparts(fns{j});
  [dd nm2]=fileparts(dd);
  nms{j}=['gc/' nm2 '/' nm1];
  ffn=['~/Dropbox/Jim-Claudia-GC/mat/' nm2 '/' nm1 '/time.mat'];
  if ~isfile(ffn)
    keyboard
    continue;
  end
  aa=load(ffn);
  aa=ded_emle_info(nms{j},aa.a);
  a=struct_array_append(a,aa,[],quiet);
  k=k+1;
  fn{k}=fns{j};
end
nm={};
n=length(a);
for j=1:n
  [dd nm{j}]=fileparts(fn{j});
  nm{j}=nm{j}(1:3);
end

[ff f]=sortrows([[a.Sc];[a.dimple];[a.Wx]]');

a=a(f);
nm={nm{f}};
fns={fns{f}};
nms={nms{f}};

n=length(a);
for j=1:n
  f=find(a(j).t>=trg(1) & a(j).t<=trg(2));
  fnm=setdiff(fieldnames(a),'nms');
  for k=1:length(fnm)
    if length(a(j).(fnm{k}))>1
      a(j).(fnm{k})=a(j).(fnm{k})(f);
    end
  end
end

Wx=[a.Wx];
dimple=[a.dimple];
Wx(dimple>0)=Wx(dimple>0);

for j=1:n
  nm{j}=sprintf('%s %5.3f %6.2f',nm{j},Wx(j),dimple(j));
end


for j=1:n
  E{j}=abs(a(j).E1);
  if length(E{j})~=length(a(j).t)
    E{j}=repmat(NaN,size(a(j).t));
  end
end
col=[0    0    1
     1    0    0
     0.25 0.25 0.25
     0    0.75 0.75
     0.75 0    0.75
     0.75 0.75 0
     0    0.75    0];
col=[col;col;col];

Sc=[a.Sc];
lc=jet(n);
f=find(dimple==0 & Sc==1);[fff ff]=sort(Wx(f));    lc(f(ff),1:3)=col(1:length(f),:);
f=find(dimple>0);         [fff ff]=sort(dimple(f));lc(f(ff),1:3)=col(1:length(f),:);

for j=1:n
  ls{j}='-';
  lw(j)=a(j).Sc;
  if a(j).noise==0.005 & a(j).Sc==1 & a(j).dimple==0
    lw(j)=1;
    ls{j}='--';
   end
end

f=figure;clf;preprint(sz,fnt);colour_lines;hold('on');
for j=1:n
  h(j)=plot(a(j).t,E{j},'linewidth',lw(j),'color',lc(j,:),'linestyle',ls{j});
end

axis([trg 0 inf]);
legend(h,nm,'location','best','fontsize',6);
xlabel('t');
ylabel('Entropy');
set(gca,'box','on');
print(f,'-depsc2',[dfig 'Entropy.eps']);

f=figure;clf;preprint(sz,fnt);colour_lines;hold('on');
for j=1:n
  m=max(1,round(0.1*length(a(j).t)/(max(a(j).t)-min(a(j).t))));
  mt=(a(j).t(1:end-m)+a(j).t(1+m:end))/2;
  dE=(E{j}(1:end-m)-E{j}(1+m:end))./(a(j).t(1:end-m)-a(j).t(1+m:end));
  h(j)=plot(mt,dE,'linewidth',lw(j),'color',lc(j,:),'linestyle',ls{j});
end
axis([trg 0 inf]);
legend(h,nm,'location','best','fontsize',6);
xlabel('t');
ylabel('Entrainment');
set(gca,'box','on');
print(f,'-depsc2',[dfig 'Entrainment.eps']);


f=figure;clf;preprint(sz,fnt);colour_lines;hold('on');
for j=1:n
  h(j)=plot(a(j).t,a(j).vv,'linewidth',lw(j),'color',lc(j,:),'linestyle',ls{j});
end
axis([trg 0 inf]);
legend(h,nm,'location','best','fontsize',6);
xlabel('t');
ylabel('vv');
set(gca,'box','on');
print(f,'-depsc2',[dfig 'vv.eps']);


f=figure;clf;preprint(sz,fnt);colour_lines;hold('on');
for j=1:n
  h(j)=plot(a(j).t,a(j).bz(1)-a(j).bz,'linewidth',lw(j),'color',lc(j,:),'linestyle',ls{j});
end
axis([trg 0 inf]);
legend(h,nm,'location','best','fontsize',6);
xlabel('t');
ylabel('-pe');
set(gca,'box','on');
print(f,'-depsc2',[dfig 'bz.eps']);

f=figure;clf;preprint(sz,fnt);colour_lines;hold('on');
for j=1:n
  s=ded_read_stats(nms{j});
  e=max(max(0,s.maxb-1),-s.minb);
  h(j)=plot(s.t,e,'linewidth',lw(j),'color',lc(j,:),'linestyle',ls{j});
end
axis([trg 0 inf]);
legend(h,nm,'location','best','fontsize',6);
xlabel('t');
ylabel('b overshoot');
set(gca,'box','on');
print(f,'-depsc2',[dfig 'bover.eps']);


f=figure;clf;preprint(sz,fnt);colour_lines;hold('on');
for j=1:n
  s=ded_read_stats(nms{j});
  p=ded_read_param(nms{j});
  e=max(s.maxv,-s.minv);
  h(j)=plot(s.t,e,'linewidth',lw(j),'color',lc(j,:),'linestyle',ls{j});
end
axis([trg 0 inf]);
legend(h,nm,'location','best','fontsize',6);
xlabel('t');
ylabel('v max');
set(gca,'box','on');
print(f,'-depsc2',[dfig 'vmax.eps']);


