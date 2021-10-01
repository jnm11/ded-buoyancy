function [nmb t]=ded_fix_time(nm,typ)
%[nmb t]=ded_fix_time('gc/ccle/046','b');
nmb={};
t=[];

fnt='~/misc/ded-3d-rm';
unix(sprintf('touch %s',fnt));
unix(sprintf('chmod a+x  %s',fnt));

p=ded_read_param(nm);
if isempty(p)
  disp(sprintf('ded_fix_stats: %s: no param file',nm));
  return;
end
% $$$ if strcmp(p.status,'Running')
% $$$   disp(sprintf('ded_fix_general: %s: Simulation is running',nm));
% $$$   return;
% $$$ end  

nmb=ded_get_fn(nm,typ);

if isempty(nmb)
  disp(sprintf('ded_fix_general: %s: no files of type %s',nm,typ));
  return;
end

nmmat=sprintf('%s/results/%s/%s-time.mat',ded_dedalus_data_dir,nm,typ);
if isfile(nmmat) & file_nt(nmmat,nmb{end})
  load(nmmat);
  return;
end

tn =1;
tt={};
jj=[];
for k=1:length(nmb)
  fn=nmb{k};
  disp(sprintf('ded_fix_time: %s',fn));
  try
    try
      tt=h5read(fn,'/scales/sim_time');
    catch
      tt=h5read(fn,'/t');
    end
   if length(tt)~=1 
     disp(sprintf('ded_fix_time: %s contains %u times',fn,length(tt{k})));
     tn=0;
   else
     t(end+1)=tt;
     jj(end+1)=k;
   end
  end
end
nmb={nmb{jj}};
if length(t)==1
  return;
end
if tn==0
  disp(sprintf('ded_fix_time: %s: files contain more than one time',nm));
  return;
end
dtt=['dt' typ];
dt=[];
if isfield(p,dtt)
  if p.(dtt)>0
    dt=p.(dtt);
  end
end
dtt=['dtj' typ];
if isfield(p,dtt)
  if p.(dtt)>0
    dt=p.(dtt);
  end
end

if isempty(dt)
  dt=round(1e3*diff(sort(t)))/1e3;
  dt=mode(dt(dt>0));
end

rg=round(min(t)/dt):round(max(t)/dt);
f=repmat(1==1,size(nmb));
for j=1:length(rg)
  g=findmin(abs(t-rg(j)*dt));
  f(g)=0;
end
f=find(f);
n=length(f);
g=1:length(nmb);

if n>0
  DD=[ded_dedalus_data_dir '/'];
  nmr=cellstrprefix('~/',cellstrremoveprefix({nmb{f}},DD));
  for j=1:n
    cmd=sprintf('/bin/rm -f %s',nmr{j});
    disp(sprintf('%s # t=%8.3f ',cmd,t(f(j))));
    unix(sprintf('echo "%s" >> %s',cmd,fnt));
    %unix(cmd);
  end
  g(f)=[];
  nmb={nmb{g}};
  t=t(g);
end
[t f]=sort(t);
nmb={nmb{f}};

mkdirifnotexist(nmmat);
save(nmmat,'nmb','t');
