function a=ded_gc_avg_le(nm,typ,out,Xf,Uf,minx,trg)

fnmat = sprintf('%s/results/%s/%s.mat',ded_dedalus_data_dir,nm,out);

for k=1:length(typ)
  fns{k}=ded_get_fn(nm,typ{k});
end

if isfile(fnmat)
  return;
end
if file_nt(fnmat,fns{1}{end}) & false
  disp(sprintf('ded_gc_avg_le:        "%s" up to date ',fnmat));
  if nargout>1
    load(fnmat);
  end
  return;
end

disp(sprintf('ded_gc_avg_le:        "%s" making ',fnmat));
 
p=ded_read_param(nm);
c=ded_coord(nm);
parity=ded_read_parity(nm);


for k=1:length(fns)
  tt{k}=ded_get_times(fns{k});
  if k==1;
    t=tt{k};
  else
    t=intersect(t,tt{k});
  end
end

t=t(t>=trg(1)&t<=trg(2));
if length(t)<2
  return;
end
n=length(t);
if isempty(c)
   aa=ded_read_hdf(fns{1}{1});
   c=ded_coord('gc/emle/019');
   c.Jx=aa.x;
   c.Jz=aa.z;
   c.dd=[mean(diff(aa.x)) 0  mean(diff(aa.z))]; 
   %c.c='xyz';
   c.NJz=length(c.Jz);
   c.NJx=length(c.Jx);
   %s.dim=3;
   clear('aa');
end

s=c.Jx-p.L+4; % new coordinates

fnm={'b','u','w',...
     'dbdx',  'dbdz', 'dudx', 'dudz', 'dwdx', 'dwdz',...
     'udbdx','wdbdz','ududx','wdudz','udwdx','wdwdz',...
     'bb','bu','bw','uu','uw','ww',...
     'bD','S','R','Q'}; 

for j=1:length(fnm)
  a.(fnm{j})=zeros(c.NJz,c.NJx);
end
m=zeros(1,c.NJx);
A.c='xz';
dt=1;

for k=1:length(fns)
  jj(k)=1;
  aa{k}=ded_read_hdf(fns{k}{1});
  if isfield(aa{k},'sim_time')
    aa{k}.t=aa{k}.sim_time;
  end
end


for j=1:n
  disp(t(j));
  A=struct('c','xz');
  for k=1:length(fns)
    if length(tt{k})>length(fns{k})
      ii=find(aa{k}.t==t(j));
      while isempty(ii)
        jj(k)=jj(k)+1;
        aa{k}=ded_read_hdf(fns{k}{jj(k)});
        if isfield(aa{k},'sim_time');aa{k}.t=aa{k}.sim_time;end
        ii=find(aa{k}.t==t(j));
      end
      if isfield(aa{k},'dt');dt=aa{k}.dt(ii);end
      nn=intersect(fnm,fieldnames(aa{k}));
      for i=1:length(nn)
        if isfield(aa{k},nn{i});A.(nn{i})=squeeze(aa{k}.(nn{i})(:,:,:,ii))/dt/p.W;end
      end
    else
      aa{k}=ded_read_hdf(fns{k}{tt{k}==t(j)});
      if isfield(aa{k},'dt');dt=aa{k}.dt;end
      nn=intersect(fnm,fieldnames(aa{k}));
      for i=1:length(nn)
        if isfield(aa{k},nn{i});
          A.(nn{i})=aa{k}.(nn{i})/dt/p.W;end
      end
    end
  end
  
  A=ded_add_diff(A,p,parity,c,{'b','u','w'});
  X=Xf(t(j));
  U=Uf(t(j));
  x=c.Jx-X;
  f=find(s>=minx-X&s<=x(end)&s>=x(1));
  nn=intersect(fnm,fieldnames(A));
  for i=1:length(nn)
    B.(nn{i}) = interp1(x,   A.(nn{i})',s(f),'linear')';
  end
  m(f)         = m(f)+1;
  if ~isfield(A,'udbdx');  B.udbdx = B.u.*B.dbdx; end;
  if ~isfield(A,'wdbdz');  B.wdbdz = B.w.*B.dbdx; end;
  if ~isfield(A,'ududx');  B.ududx = B.u.*B.dudx; end;
  if ~isfield(A,'wdudz');  B.wdudz = B.w.*B.dudx; end;
  if ~isfield(A,'udwdx');  B.udwdx = B.u.*B.dwdx; end;
  if ~isfield(A,'wdwdz');  B.wdwdz = B.w.*B.dwdx; end;
  if ~isfield(A,'bu');     B.bu    = B.b.*B.u; end;
  if ~isfield(A,'bw');     B.bw    = B.b.*B.w; end;
  if ~isfield(A,'uu');     B.uu    = B.u.*B.u; end;
  if ~isfield(A,'uw');     B.uw    = B.u.*B.w; end;
  if ~isfield(A,'ww');     B.ww    = B.w.*B.w; end;
    
  nn=intersect(fnm,fieldnames(B));
  for i=1:length(nn)
    a.(nn{i})(:,f) = a.(nn{i})(:,f)+B.(nn{i});
  end

  
  a.bD(:,f)    = a.bD(:,f)  + B.dbdx.^2+B.dbdz.^2;
  a.S(:,f)     = a.S(:,f)   + B.dudx.^2+B.dwdz.^2 + (B.dudz+B.dwdx).^2/2;
  a.R(:,f)     = a.R(:,f)   + (B.dudz-B.dwdx).^2/2;
end


f=min(find(m>0)):max(find(m>0));
m=m(f);
for j=1:length(fnm)
  a.(fnm{j})=a.(fnm{j})(:,f)./m;
end

a.Q=a.S-a.R; % Okubo-Weiss
a.m=m;
a.t1=t(1);
a.t2=t(end);
a.n=n;
a.nm=nm;
a.z=c.Jz;
a.x=s(f)';

mkdirifnotexist(fnmat);
save(fnmat','a');

return;

cd('~/');
fns=cellstr_ls('gc/le/a/*/*/*/status',[],'dir');
minx=24;
trg=[8 24];
for j=1:length(fns)
  nm=fns{j};
  a=ded_gc_profiles_slice(fns{j});
  ded_gc_avg_le(nm,a.Xf,a.Uf,minx,trg);
end

