function a=ded_read_3d(dd,nm,typ)

if nargin<3
  typ=[];
end

fns=ded_get_fn(dd,nm);

if isempty(fns)
  a=[];
  return;
end

if ~isempty(typ)
  if ~iscell(typ)
    typ={typ};
  end
end

%h5disp(fn);
b=h5info(fns{1});

for k=1:length(b.Groups)
  nmx=[b.Groups(k).Groups(1).Name '/' b.Groups(1).Groups(1).Datasets.Name];
  nmy=[b.Groups(k).Groups(2).Name '/' b.Groups(1).Groups(2).Datasets.Name];
  nmz=[b.Groups(k).Groups(3).Name '/' b.Groups(1).Groups(3).Datasets.Name];
  if strcmp(nmx(1:9),'/scales/x')
    break;
  end
end

if isempty(typ)
  tnm=b.Groups(2).Name;
  typ={b.Groups(2).Datasets.Name};
else
  tnm='/tasks';
end


a=struct;
for k=1:length(fns)
  fn=fns{k};
  if k==1
    a.x=h5read(fn,nmx);
    a.y=h5read(fn,nmy);
    a.z=h5read(fn,nmz);
    a.nx=length(a.x);
    a.ny=length(a.y);
    a.nz=length(a.z);
  end
  t{k}=h5read(fn,'/scales/sim_time');
  nt=length(t{k});
  sz=[a.nz a.ny a.nx nt];
  for j=1:length(typ)
    x=h5read(fn,[tnm '/' typ{j}]);
    c(k).(typ{j}) = reshape(x,sz);
  end
end
for j=1:length(typ)
  a.(typ{j})=cat(4,c.(typ{j}));
end
a.t=cat(1,t{:});
a.nt=length(a.t);

return

a=ded_read_3d('plume/003','noise');
plot(a.t,sqrt(squeeze(sum(sum(sum(a.psix.^2,1),2),3))));
subplot(2,1,1);plot(a.x,squeeze(max(max(max(abs(a.psix),[],1),[],2),[],4)));
subplot(2,1,2);plot(a.t,squeeze(max(max(max(abs(a.psix),[],1),[],2),[],3)));


