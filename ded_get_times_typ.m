function [times nms]=ded_get_times_typ(nm,typ)

DHOME=ded_dedalus_data_dir;
nmmat=[DHOME '/' nm '/' typ '/times.mat'];
nms=ded_get_fn(nm,typ);
if all(file_nt(nmmat,nms))
  disp(sprintf('ded_get_times_typ: loading %s',nmmat));
  load(nmmat);
  return
end
disp(sprintf('ded_get_times_typ: remaking %s',nmmat));
[times nms]=ded_get_times(nms);
save(nmmat,'times','nms');

return;
[times nms]=ded_get_times_typ('gc/emle/017','y');

