function  ded_make_param(nm)

%ded_make_param(cellstr_ls('gc/f6/[hg]/*/*/*/*/status',[],'dir'));
if nargin<3; trg=[];end;

if iscell(nm)
  for j=1:length(nm)
    ded_make_param(nm{j});
  end
  return;
end

fnp = sprintf('%s/results/%s/param.mat',ded_dedalus_data_dir,nm);
fnc = sprintf('%s/results/%s/coord.mat',ded_dedalus_data_dir,nm);
fnr = sprintf('%s/results/%s/parity.mat',ded_dedalus_data_dir,nm);
p=ded_read_param(nm);
if isempty(p)
  return;
end
mkdirifnotexist(fnp);
save(fnp,'p');

c=ded_coord(nm);
save(fnc,'c');

p=ded_read_parity(nm);
save(fnr,'p');
