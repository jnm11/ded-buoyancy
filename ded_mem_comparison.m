%#!/bin/bash
%mkdir -p ~/ded-mem/poseidon
%a=$(ls ~/pm/test-stampede/test)
%for b in $a
%  do
%  cd ~/pm/test-stampede/test/$b
%  c=$(ls -1t *.out | head -n 1)
%  /bin/cp -f $c  ~/ded-mem/S5-$b.out
%  done
%  
%for a in $(seq -f %02g 1 36);do  rsync -vap ~/pm/poseidon/$a/out   asahi:ded-mem/poseidon-$a.out ; done
%for a in $(seq -f %02g 1 54);do  rsync -vap ~/pm/hamilton/$a/out   asahi:ded-mem/hamilton-$a.out ; done
function a=ded_mem_comparison
%for a in $(seq -f %02g 1 4);do for b in $(seq -f %02g 1 24);do  rsync -vap archer:pm/A2/$a-$b/pm-A2-$a-$b.o* ~/ded-mem/A2-$a-$b.out;done;done
%for a in $(seq -f %02g 1 4);do for b in $(seq -f %02g 1 96);do  rsync -vap ~/pm/S2/$a-$b/*.out ~/ded-mem/S2-$a-$b.out;done;done
%for a in $(seq -f %02g 1 4);do for b in $(seq -f %02g 1 48);do  rsync -vap ~/pm/S3/$a-$b/*.out ~/ded-mem/S3-$a-$b.out;done;done;rsync -uvap ~/ded-mem asahi:
%for a in $(seq -f %02g 1 4);do for b in $(seq -f %02g 1 48);do  rsync -vap ~/pm/S4/$a-$b/*.out ~/ded-mem/S4-$a-$b.out;done;done;rsync -uvap ~/ded-mem asahi:
%for a in $(seq -f %02g 1 4);do for b in $(seq -f %02g 1 48);do  rsync -vap ~/pm/T2/$a-$b/*.out ~/ded-mem/T2-$a-$b.out;done;done
%for a in $(seq -f %02g 1 8);do for b in $(seq -f %02g 1 16);do  rsync -vap ~/pm/H2/$a-$b/*.out ~/ded-mem/H2-$a-$b.out;done;done
%for a in $(seq -f %02g 1 4);do for b in $(seq -f %02g 1 36);do  rsync -vap ~/pm/P2/$a-$b/*.out ~/ded-mem/P2-$a-$b.out;done;done
% rsync -uvap hamilton:ded-mem ~/
% rsync -vap asahi:ded-mem ~/
cp={'H2','P2','S3','S4'};
cp={'S5'};

for k=1:length(cp)
  fns=cellstr_ls(sprintf('~/ded-mem/%s-*.out',cp{k}));
  nfns=length(fns);
  j=0;
  for i=1:length(fns)
    f=fopen(fns{i},'r');
    if f==-1
      disp(sprintf('could not open %s',fns{i}));
      continue
    end
    [dd fn]=fileparts(fns{i});
    j=j+1;
    a(k).nm{j}=fn;
    a(k).gc(j)=NaN;    
    a(k).nodes(j)=NaN;    
    a(k).tasks(j)=NaN;    
    a(k).tpn(j)=NaN;    
    a(k).ttime(j)=NaN;    
    a(k).t(j)=NaN;    
    a(k).Nx(j)=NaN;    
    a(k).sw(j)=NaN;    
    while ~feof(f)
      l=fgetl(f);
      [A C]=sscanf(l,'JOB_NUM_NODES: %f');      if C>0 ; a(k).nodes(j)=A ; end
      [A C]=sscanf(l,'JOB_NUM_TASKS: %f');      if C>0 ; a(k).tasks(j)=A ; end
      [A C]=sscanf(l,'JOB_TASKS_PER_NODE: %f'); if C>0 ; a(k).tpn(j)=A ;   end
      a(k).Nx(j)    = getx(l,'Nx:',a(k).Nx(j));      
      a(k).gc(j)    = getx(l,'Grid cells',a(k).gc(j));      
      a(k).nodes(j) = getx(l,'NUM_NODES',a(k).nodes(j));      
      a(k).tasks(j) = getx(l,'NUM_TASKS',a(k).tasks(j));      
      a(k).tpn(j)   = getx(l,'TASKS_PER_NODE',a(k).tpn(j)); 
      a(k).ttime(j) = gett(l,'main loop',a(k).ttime(j));      
      a(k).t(j)     = getx(l,' t; ',a(k).t(j));   
      a(k).sw(j)    = getx(l,'sim/wall',a(k).sw(j));
      g=strfind(l,'size:');
      if ~isempty(g)
        a(k).nt(j)=str2num(l(g+6:end));
      end
      g=strfind(l,'GB ,');
      if ~isempty(g)
        x=str2num(l(g-9:g-1));
        a=cstats(a,k,j,'G',x);
      end
      g=strfind(l,', e:');
      if ~isempty(g)
        x=str2num(l(g+4:g+11));
        a=cstats(a,k,j,'e',x);
      end  
    end
    fclose(f);
    if ~isfinite(a(k).tasks)
      keyboard;
    end
  end
  a(k).G=a(k).G./a(k).Gn;
  a(k).Gs=sqrt(a(k).Gs./a(k).Gn-a(k).G.^2);
  a(k).e=a(k).e./a(k).en;
  a(k).es=sqrt(a(k).es./a(k).en-a(k).e.^2);
  disp(sprintf('%s read %u',cp{k},j));
end
return;
for j=1:length(a)
  [t f]=sort(a(j).tasks);
  N=length(t);
  nm=fieldnames(a(j));
  for k=1:length(nm)
    b=nm{k};
    if iscell(a(j).(b))
      a(j).(b)={a(j).(b){end+1:N}};
      a(j).(b)={a(j).(b){f}};
    else
      a(j).(b)(end+1:N)=NaN;
      a(j).(b)=a(j).(b)(f);
    end
  end
end






na=length(a);
col=jet(na);
col=[     1     0     0
          0     1     0
          0     0     1
          0     0     0
          0     1     1];

p.lc='-';
[ugc xxx ukgc]=unique(a.gc);
[und xxx uknd]=unique(a.nodes);
uk=ukgc+(max(ukgc)+1)*uknd;
ugc=ugc(isfinite(ugc));
for j=1:length(ugc)
  figure(j);
  f=find(a.gc==ugc(j));
  clf;groupplot(a.tpn(f),a.G(f),a.nodes(f),p);
  title(sprintf('Millions of grid cells %f',ugc(j)));
end

keyboard


clf;hold('on');h=zeros(1,na);
for j=1:na
  u=unique(a(j).nodes);
  for k=1:length(u)
    f=find(a(j).nodes==u(k));
    h(j)=plot(a(j).tasks(f),a(j).G(f)./a(j).tpn(f),'color',col(j,:));
  end
end
legend(h,cp);

figure(2);
clf;hold('on');h=zeros(1,na);
for j=1:na
  u=unique(a(j).nodes);
  for k=1:length(u)
    f=find(a(j).nodes==u(k));
    h(j)=plot(a(j).tasks(f),a(j).G(f),'color',col(j,:));
  end
end
axis('tight');aa=axis;aa([1 3])=0;
axis(aa);
plot(aa(1:2),64*[[1;1] [2;2] [3;3]]);
legend(h,cp);




figure(3);
clf;hold('on');h=zeros(1,na);
for j=1:na
  u=unique(a(j).nodes);
  for k=1:length(u)
    f=find(a(j).nodes==u(k));
    h(j)=plot(a(j).tasks(f),a(j).e(f)./a(j).nodes(f),'color',col(j,:));
  end
end
figure(4);
clf;hold('on');h=zeros(1,na);
for j=1:na
  u=unique(a(j).nodes);
  for k=1:length(u)
    f=find(a(j).nodes==u(k));
    h(j)=plot(a(j).tasks(f),24*60^2./a(j).ttime(f),'color',col(j,:));
  end
end
legend(h,cp);



return;

function y=getx(l,m,y)
f=strfind(l,m);
if ~isempty(f)
  x=str2num(strtok(l(f+length(m):end)));
else
  x=NaN;
end
if isempty(x)
  x=NaN;
end
if x>0
  y=x;
end

function y=gett(l,m,y)
L=l;
f=strfind(l,m);
x=NaN;
if ~isempty(f)
  l=l(f+length(m):end);
  if ~isempty(l)
    if l(1)==':';l(1)=[]; end
    x=[];
    while ~isempty(l)
      [tok l]=strtok(l,' :');
      x(end+1)=str2num(tok);
    end
    x=flip(x);
    x(end+1:4)=0;
    x=sum( x.*[1 60 60^2 24*60^2]);
  end
end
if isempty(x)
  x=NaN;
end
if x>0
  y=x;
end

function a=cstats(a,k,j,nm,x)
if ~isfield(a(k),nm)
  a(k).([nm ''])(j)    = x;
  a(k).([nm 's'])(j)   = x^2;
  a(k).([nm 'n'])(j)   = 1;
  a(k).([nm 'min'])(j) = x;
  a(k).([nm 'max'])(j) = x;
elseif length(a(k).(nm))<j 
  a(k).([nm ''])(j)    = x;
  a(k).([nm 's'])(j)   = x^2;
  a(k).([nm 'n'])(j)   = 1;
  a(k).([nm 'min'])(j) = x;
  a(k).([nm 'max'])(j) = x;
else
  a(k).([nm ''])(j)    = a(k).([nm ''])(j)+x;
  a(k).([nm 's'])(j)   = a(k).([nm 's'])(j)+x^2;
  a(k).([nm 'n'])(j)   = a(k).([nm 'n'])(j)+1;
  a(k).([nm 'min'])(j) = min(a(k).([nm 'min'])(j),x);
  a(k).([nm 'max'])(j) = max(a(k).([nm 'max'])(j),x);
end  
return;
