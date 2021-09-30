function a=ded_slice(nm,fnm,X,trg,typ,out)
%r=ded__slice(nm,X,trg,display) Calculate profiles at a fixed position as a they vary in time

if iscell(nm)
  for j=1:length(nm)
    ded_slice(nm{j},X,trg,typ,out,display);
  end
  return;
end

fnmat=sprintf('%s/results/%s/%s.mat',ded_dedalus_data_dir,nm,out);
p=ded_read_param(nm);
dt=1;
avrg=false;
switch(typ)
  case 'bu'
    fns{1}=ded_get_fn(nm,'u');
    fns{2}=ded_get_fn(nm,'b');
  case 'ay'
    fns{1}=ded_get_fn(nm,'ay');
    avrg=true;
  case 'y'
    fns{1}=ded_get_fn(nm,'y');
   otherwise
     error(sprintf('ded_gc_profiles_slice: unknown type %s',typ));
end

if file_nt(fnmat,fns{1}{end})
  disp(sprintf('ded_slice: %s up to date',nm));
  if nargout>0
    r=load(fnmat);
    r=r.a;
  end
  return;
end
if nargin==1 & isfile(fnmat)
  a=load(fnmat);
  a=a.a;
  return;
end

disp(sprintf('ded_slice: %s making',nm));

p=ded_read_param(nm);
if isempty(p);
  disp(sprintf('ded_slice: %s no data'));
  return;
end
c=ded_coord(nm);
if isempty(c)
  b=ded_read_hdf(fns{1}{1});
  c.Jx=b.x;
  c.Jz=b.z;
end
W=p.W;

for k=1:length(fns)
  n(k)=length(fns{k});
end
n=max(n);

a.x=c.Jx;
a.z=c.Jz;

nz=length(a.z);
nx=length(a.x);
a.t=zeros(1,n);
for i=1:length(fnm)
  a.(fnm{i})=zeros(nz,n);
end

for k=1:length(fns)
  for j=1:n
    b=ded_read_hdf(fns{k}{j});
    if avrg; dt=b.dt;end
    for i=1:length(fnm)
      if isfield(b,fnm{i})
        a.(fnm{i})(:,j)=interp1(a.x,b.(fnm{i})',X,'linear')/(dt*W);
      end
    end
    if isfield(b,'t')
      a.t(j)=b.t;
    else
      a.t(j)=(b.t1+b.t2)/2;
    end
  end
end

mkdirifnotexist(fnmat);
save(fnmat,'a');

return;
