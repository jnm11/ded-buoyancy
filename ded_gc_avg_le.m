function a=ded_gc_avg_le(nm,typ,out,Xf,Uf,minx,trg)

fnmat = sprintf('%s/results/%s/%s.mat',ded_dedalus_data_dir,nm,out);

fns=ded_get_fn(nm,typ);

if file_nt(fnmat,fns{end})
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


tt=ded_get_times(fns);
f= find(tt>=trg(1) & tt<trg(2));
fns={fns{f}};
t=tt(f);
if length(t)<2
  return;
end
n=length(t);

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


avg=fns{j}(1)=='a';

for j=1:length(fns)
  A=ded_read_hdf(fns{j});
  if isfield(A,'t1');A.t=(A.t1+A.t2)/2;end
  disp([t(j) A.t]);
  if avg;dt=A.dt;end
  nn=intersect(fnm,fieldnames(A));
  if dt*p.W~=1 
    for i=1:length(nn)
      A.(nn{i})  = A.(nn{i})/(dt*p.W);
    end
  end

  A=ded_add_diff(A,p,parity,c,{'b','u','w'});

  % Transform u velocities
  X=Xf(t(j));
  U=Uf(t(j));

  if isfield(A,'bu');  A.bu = A.bu-A.b*U; end;
  if isfield(A,'uu');  A.uu = A.uu+U^2-2*A.u*U; end;
  if isfield(A,'u');   A.u  = A.u-U; end;
  nn = fieldnames(A);
  for i=1:length(nn)
    if length(nn{i})<2; continue; end
    n1=nn{i}(1);
    n2=nn{i}(2);
    n3=nn{i}(2:end);
    if n1~='u' | n2=='u';continue;end
    if isfield(A,n3)
      A.(nn{i})  = A.(nn{i})-A.(n3)*U;
    else
      A=rmfield(A,nn{i});
    end
  end
  
 
  
  x=c.Jx-X;
  f=find(s>=minx-X&s<=x(end)&s>=x(1));
  nn=intersect(fnm,fieldnames(A));
  for i=1:length(nn)
    B.(nn{i}) = interp1(x,   A.(nn{i})',s(f),'linear')';
  end
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
    a.(nn{i})(:,f) = a.(nn{i})(:,f)+dt*B.(nn{i});
  end
  m(f)         = m(f) + dt*p.W;

  a.bD(:,f)    = a.bD(:,f)  + dt*(B.dbdx.^2+B.dbdz.^2);
  a.S(:,f)     = a.S(:,f)   + dt*(B.dudx.^2+B.dwdz.^2 + (B.dudz+B.dwdx).^2/2);
  a.R(:,f)     = a.R(:,f)   + dt*(B.dudz-B.dwdx).^2/2;
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

