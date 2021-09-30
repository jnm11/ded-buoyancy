function c=ded_hdf(nm,a,d)
if nargin<3
  d=[];
end
if isempty(d)
  fn=[ded_dedalus_data_dir '/' nm '/' a '.hdf5'];
else
  fn=[ded_dedalus_data_dir '/' nm '/' d '/' a '.hdf5'];
end 
if isfile(fn)
  c=ded_read_hdf(fn);
else
  c=[];
  disp(sprintf('ded_hdf: %s does not exist',fn));
end
