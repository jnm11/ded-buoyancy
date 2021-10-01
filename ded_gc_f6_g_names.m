function nms=ded_gc_f6_g_names
cd(ded_dedalus_data_dir);
nm1=cellstr_ls('gc/f6/g/*/*/*/[23]*/param.h5',[],'dir');
nm2=cellstr_ls('gc/f6/i/*/*/*/*b/param.h5',[],'dir');
if  isempty(nm1)
  cd('results');
  nm1=cellstr_ls('gc/f6/g/*/*/*/[23]*/profile.mat',[],'dir');
  nm2=cellstr_ls('gc/f6/i/*/*/*/*b/profile.mat',[],'dir');
end
nms=cat(1,nm1,nm2);
