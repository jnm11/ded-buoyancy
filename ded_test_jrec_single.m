function [E a c avrg]=ded_test_jrec_single(nm,typs)
%ded_test_jrec_single('gc/jrec/00');nm='gc/jrec/00';typs={'s','b','u','v','w','p'};

if nargin==0
  fns=cellstr_ls([ded_dedalus_data_dir '/gc/jrec/*');
  fns=cellstrremoveprefix(fns,ded_dedalus_data_dir('HOME'));
  for j=1:length(fns)
    E(j)=ded_test_jrec_single(fns{j});
  end
  return;
end

dd=ded_dedalus_data_dir;
if nargin==0
  nms=cellstr_ls([dd '/gc/jrec/*']);
  nms=cellstrremoveprefix(nms,[dd '/']);
  for j=1:length(nms)
    ded_test_jrec_single(nms{j});
  end
  return;
end
if nargin<2
  typs=[];
end
p=ded_read_param(nm);

if isempty(typs)
  typs={'s','b','u','v','w','p'};
end

ityp={'x','y','z','xy','xz','yz','xyz'};
k=0;
for kk=1:length(ityp)
  fns =cellstr_ls(sprintf('%s/%s/%s/%s-*.hdf5',ded_dedalus_data_dir,nm,ityp{kk},ityp{kk}));
  if isempty(fns)
    continue;
  end
  k=k+1;
  f(k)=kk;
  for j=1:length(fns)
    aa=ded_read_hdf(fns{j});
    I.(ityp{kk})(j)=aa;
  end
end
ityp={ityp{f}};

k=0;
f=[];
for kk=1:length(typs)
  fns =cellstr_ls(sprintf('%s/%s/%s/%s-*.hdf5',ded_dedalus_data_dir,nm,typs{kk},typs{kk}));
  if isempty(fns)
    continue;
  end
  k=k+1;
  f(k)=kk;
  for j=1:length(fns)
    aa=ded_read_hdf(fns{j});
    a(j).(typs{kk})=aa.(typs{kk});
    a(j).t=aa.t;
    a(j).dt=aa.dt;
  end
end
typs={typs{f}};


dx=p.L/p.Nx;
dy=p.W/p.Ny;
dz=p.H/2*ichebIw(p.Nz);
dim =2+(p.Ny>1);

for k=1:length(ityp)
  for j=1:length(I.(ityp{k}))
    ttp=ityp{k};
    x=I.(ttp)(j);
    if dim==2
      switch(ttp)
        case 'y'
          ttp='';
        case 'xy'
          ttp='x';
        case 'yz'
          ttp='z';
      end
    end
    y=a(j);
    if y.t~=x.t1 | y.t~=x.t2 
      disp(sprintf('ded_test_jrec_single: times do not match %f %f %f',y.t,x.t1,x.t2));
    end
    fld=setdiff(fieldnames(y),{'t','dt'});
    ff=fld;
    for i=1:length(fld)
      for ii=1:i
        ff{end+1}=sort([fld{i} fld{ii}]);
        y.(ff{end}) = y.(fld{i}).*y.(fld{ii});
      end
    end
    fld=intersect(ff,fieldnames(x));
    for i=1:length(fld)
      switch(ttp)
        case 'x'
          y.(fld{i})=sum(dx.*y.(fld{i}),dim);
        case 'y'
          y.(fld{i})=sum(dy.*y.(fld{i}),2);
        case 'z'
          y.(fld{i})=sum(dz.*y.(fld{i}),1);
        case 'xy'
          y.(fld{i})=sum(dy.*sum(dx.*y.(fld{i}),3),2);
        case 'xz'
          y.(fld{i})=sum(dz.*sum(dx.*y.(fld{i}),dim),1);
        case 'yz'
          y.(fld{i})=sum(dz.*sum(dy.*y.(fld{i}),2),1);
        case 'xyz'
        y.(fld{i})=sum(dz.*sum(dy.*sum(dx.*y.(fld{i}),3),2),1);
      end
    end
    for i=1:length(fld)
      e(i)=max(abs(x.(fld{i})(:)-y.(fld{i})(:)));
    end
    %disp(sprintf('%s %i %f',ityp{k},j,max(e)));
    if max(e)>1e-6
      keyboard;
    end
  end   
end



ttt={};
E=-inf;
ff=[];

for k=1:length(typs)
  typ=typs{k};
  fns1 =cellstr_ls(sprintf('%s/%s/%s/%s-*.hdf5',ded_dedalus_data_dirnm,typ,typ));
  aa=ded_avrg_single(fns1);
  if k==1
    avrg=aa;
  else
    avrg.(typ)=aa.(typ);
    avrg.typs{end+1}=aa.typs{1};
  end
  for j=1:k
    typ=typs{j};
    fns2 =cellstr_ls(sprintf('%s/%s/%s/%s-*.hdf5',ded_dedalus_data_dir,nm,typ,typ));
    aa=ded_avrg_double(fns1,fns2);
    for i=1:length(aa.typs)
      nn=aa.typs{i};
      avrg.(nn)=aa.(nn);
      avrg.typs{end+1}=nn;
    end
  end
end

for k=1:length(typs)
  typ=typs{k};
  fns1 =cellstr_ls(sprintf('%s/%s/%s/%s-*.hdf5',ded_dedalus_data_dir,nm,typ,typ));
  fns2 =cellstr_ls(sprintf('%s/%s/%s/%s_s*.hdf5',ded_dedalus_data_dir,nm,typ,typ));
   n=min(length(fns1),length(fns2));
  if n>0
    ff(end+1)=k;
  end
  for j=1:n
    b=ded_read_hdf(fns2{j});
    [d f1]=fileparts(fns1{j});
    [d f2]=fileparts(fns2{j});
    x=a(j).(typ);
    y=b.(typ);
    e1=max(abs(x(:)-y(:)));
    e2=sqrt(mean(abs((x(:)-y(:)).^2)));
    x=sort(x(:));
    y=sort(y(:));
    e3=max(abs(x(:)-y(:)));
    e4=sqrt(mean(abs((x(:)-y(:)).^2)));
    dt=b.sim_time-a(j).t;
    E=max([E e1 e2 e3 e4]);
    if any([e1 e2 e3 e4]>1e-10)
      disp(sprintf('ded_test_jrec_single: %s:%s %s, dt=%6.4f, e=[%9.7f %9.7f %9.7f %9.7f]',nm,f1,f2,dt,e1,e2,e3,e4));
    end
  end
end
if max(E)<1e-12
  disp([sprintf('ded_test_jrec_single: %s passed %6.2e',nm,E) sprintf(' %s',typs{ff}) ]);
end
nms=cellstr_ls(sprintf('%s/%s/a/*.hdf5',dd,nm));
c=ded_combine_javrg(nms);
if ~isempty(c)
  for j=1:length(c.typs)
    nm=c.typs{j};
    e=max(abs(c.(nm)(:)-avrg.(nm)(:)));
    if strcmp(c.typs{j},'Eb') | strcmp(c.typs{j},'Es')
      e=0;
    end
    if all(c.(nm)(:)==0)
      disp(sprintf('%2s all zero',nm));   
    elseif e>1e-10 
      disp(sprintf('%2s %8.6f',nm,e));
      keyboard
    end
    E=max([E e]);
  end
end

return;




return;


nm='gc/jrec/01';
ded_test_jrec_single(nm);
ded_test_jrec_single('gc/jrec/02');

a=ded_read_hdf('%s/gc/jrec/01/y/y_s1.hdf5',ded_dedalus_data_dir);
b=ded_read_hdf('%sgc/jrec/01/y/y-00001.hdf5'ded_dedalus_data_dir,);

nms=setdiff(intersect(fieldnames(a),fieldnames(b)),{'x','y','z'});
for j=1:length(nms)
  nm=nms{j};
  e=max(max(abs(a.(nm)(:,:,2)-b.(nm)(:,:))));
  disp(sprintf('%2s %8.6f',nm,e));
end

disp([min(a.uw(:)) max(a.uw(:)) min(b.uw(:)) max(b.uw(:))]);

fns=cellstr_ls('%s/pm/024/ayz/*hdf5',ded_dedalus_data_dir);
for j=1:length(fns)
  a=ded_read_hdf(fns{j});
  if max(a.b(:))==0
    disp(sprintf('/bin/rm -f %s',fns{j}));
  end
end
