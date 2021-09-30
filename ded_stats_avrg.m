function [m v p]=ded_stats_avrg(nm,ctp)


if ~iscell(nm)
  if any(nm=='*') 
    nm=ded_get_fn(nm,[],[],'stats');
    nm=cell_rm_suffix('/stats.hdf5',nm);
    nm=cell_rm_suffix('/',nm);
    nm=cell_rm_prefix(ded_dedalus_data_dir,nm);
    nm=cell_rm_prefix('/',nm);
  end
end

quiet=false;
m=[];v=[];p=[];
if iscell(nm)
  for j=1:length(nm)
    [mm vv pp]=ded_stats_avrg(nm{j},ctp);
    m=struct_array_append(m,mm,nm{j},quiet);
    v=struct_array_append(v,vv,nm{j},quiet);
    p=struct_array_append(p,pp,nm{j},quiet);
  end
  return;
end

s=ded_read_stats(nm);
p=ded_read_param(nm);
fnm=fieldnames(s);
T=ded_convergence_T(nm,ctp);
n=min(find(s.t>=T));

for j=1:length(fnm)
  b=fnm{j};
  m.(b)=mean(s.(b)(n:end));
  v.(b)=std( s.(b)(n:end));
end
m.t1=s.t(n);
m.t2=s.t(end);
v.t1=s.t(n);
v.t2=s.t(end);
m.T=T;
v.T=T;
m.dt=m.t2-m.t1;
v.dt=v.t2-v.t1;
m.nm=nm;
s.nm=nm;



