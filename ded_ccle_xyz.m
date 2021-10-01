function a=ded_ccle_xyz(nm)

dd=ded_dedalus_data_dir;
if nargin<1
  nms=cellstr_ls([dd '/gc/ccle/*']);
  nms=cellstrremoveprefix(nms,[dd '/']);
  for j=1:length(nms)
    a{j}=ded_ccle_xyz(nms{j});
  end
  return;
end


nmo=nm;nmo(nmo=='/')='-';
nmo=[dd '/mat/' nmo '.mat'];

fns=ded_get_fn(nm,'xyz');
a=[];
if ~isempty(fns)
  if all(file_nt(nmo,fns))
    disp(sprintf('ded_ccle_xyz: loading %s',nmo));
    load(nmo);
    return;
  end
  for j=1:length(fns)
    aa=ded_read_hdf(fns{j});
    if isempty(aa)
      break;
    end
    aa=rmfield(aa,intersect(fieldnames(aa),{'x','y','z','t2','n','i2'}));
    a=struct_array_append(a,aa);
  end
  fnm=fieldnames(a);
  for k=1:length(fnm)
    b.(fnm{k})=cat(2,[a.(fnm{k})]);
  end
  b.t=b.t1;
  b.i=b.i1;
  b=rmfield(b,{'i1','t1'});
  a=b;
  disp(sprintf('ded_ccle_xyz: saving %s',nmo));
  save(nmo,'a');
  return
end

fns=ded_get_fn(nm,'y');
if isempty(fns)
  a=[];
  return;
end

if all(file_nt(nmo,fns))
  disp(sprintf('ded_ccle_xyz: loading %s',nmo));
  load(nmo);
  return;
end
p=ded_read_param(nm);

W=[];
for j=1:length(fns)
  disp(sprintf('ded_ccle_xyz: loading %s',fns{j}));
  aa=ded_read_hdf(fns{j});
  if isempty(aa)
    continue;
  end
  if length(aa.z)~=length(W)
    W=p.L/length(aa.x)*ichebintw(length(aa.z),p.H);
  end
  fnm=intersect(fieldnames(aa),{'Ex','Ey','Ez','S','b','p','u','ub','uu','uv','uw','v','vb','vv','vw','w','ww'});
  for k=1:length(fnm)
    a(j).(fnm{k})=W*squeeze(sum(sum(aa.(fnm{k}),3),2));
  end
  a(j).bz = (W.*aa.z')*squeeze(sum(sum(aa.b,3),2));;
  a(j).uz = (W.*aa.z')*squeeze(sum(sum(aa.u,3),2));;
  if isfield(aa,'sim_time')
    a(j).t  = aa.sim_time';
    a(j).dt = aa.timestep';
    a(j).i  = aa.iteration';
  else
    a(j).t  = aa.t1;
    a(j).dt = aa.dt;
    a(j).i  = aa.i1;
  end    
  if isfield(aa,'write_number')
    a(j).wn = aa.write_number';
  elseif isfield(aa,'write_num')
    a(j).wn = aa.write_num';
  else
    a(j).wn = aa.wn;
  end
end
fnm=fieldnames(a);
for k=1:length(fnm)
  b.(fnm{k})=cat(2,[a.(fnm{k})]);
end

a=b;
disp(sprintf('ded_ccle_xyz: saving %s',nmo));
save(nmo,'a');

