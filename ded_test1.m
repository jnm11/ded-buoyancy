clear('all');
fns={'gc-ccle-024-01.mat','gc-ccle-024-04.mat','gc-ccle-024-07.mat','gc-ccle-024-10.mat'...
     'gc-ccle-024-02.mat','gc-ccle-024-05.mat','gc-ccle-024-08.mat',...
     'gc-ccle-024-03.mat','gc-ccle-024-06.mat','gc-ccle-024-09.mat','gc-ccle-024.mat'};

fns={'gc-ccle-024-07.mat',...
     'gc-ccle-024-05.mat','gc-ccle-024-08.mat',...
     'gc-ccle-024-06.mat','gc-ccle-024-09.mat','gc-ccle-024.mat'};
d1='~/Dropbox/Jim-Claudia-GC/mat/old/';

fns={'024-18-steady.mat','024-22-steady.mat','024-26-steady.mat'};
d2='~/Dropbox/Jim-Claudia-GC/mat/mat-ccle-024';


fns={'gc-ccle-024-04.mat','gc-ccle-024-05.mat','gc-ccle-024-06.mat','gc-ccle-024-07.mat','gc-ccle-024-08.mat',...
     'gc-ccle-024-09.mat','gc-ccle-024.mat','024-18-steady.mat','024-22-steady.mat','024-26-steady.mat'};

clear('a');
n=length(fns)
for j=1:n
  fn1=[d1 '/' fns{j}];
  if isfile(fn1)
    aa=load(fn1);
  end
  fn2=[d2 '/' fns{j}];
  if isfile(fn2)
    aa=load(fn2);
  end
  
  a(j).z=aa.a.z;
  a(j).x=aa.a.x;
  a(j).b=aa.a.b;
  a(j).t1=aa.a.t1;
  a(j).t2=aa.a.t2;
end

x=-6:0;
nx=length(x);
clf;
ah=jsubplot([nx 1],[0.02 0.05],[0.02 0.02],[0.02 0.05]);

%b=load('~/Dropbox/Jim-Claudia-GC/mat/024-steady.mat');

for j=1:nx
  axes(ah(j));
  cla;
  hold('on');
  for k=1:n
    f=findmin(abs(a(k).x-x(j)));
    h(j,k)=plot(a(k).b(:,f),a(k).z);
  end
  title(sprintf('%4.1f',x(j)));
end
set(ah,'ylim',[0 1],'xlim',[0 1],'yticklabels',[]);
%set(h(:,end),'linewidth',3,'color',[0 0 0]);
set(h(:,end-2:end),'linewidth',2);
legend(h(end,:),{'04','05','06','07','08','09','a','18','22','26'});
%legend(h(end,:),{'18','22','26'});





