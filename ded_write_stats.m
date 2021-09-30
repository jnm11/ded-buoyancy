function a=ded_write_stats(fn,a)
%a=ded_read_hdf('force/force_s1.hdf5');

fnm=fieldnames(a);
n=length(fnm);
if isfile(fn)
  unix(sprintf('/bin/mv -f %s %s.old',fn,fn));
end
for j=1:n
  c=fnm{j};
  h5create(fn,['/' c],Inf,'ChunkSize',2^20);
  h5write(fn,['/' c],a.(c),1,length(a.(c)));
end
%h5disp(fn)
%c=ded_read_hdf5(fn
