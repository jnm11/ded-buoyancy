function  ded_find_small_dt(nm,typ,mindt)
% Search for files only containing a small time interval
% 
% ded_find_small_dt('gc/f6/g/*/*/*/*','ay',1)
% 
% 

nms=ded_get_fn(nm,typ);
[t nms dt]=ded_get_times(nms);
nms=cellstrremoveprefix(nms,[ded_dedalus_data_dir '/']);
f=find(dt<mindt);
fp=fopen('~/misc/ded-rm','w');
fprintf(fp,'cd ~/data/dedalus\n',ded_dedalus_data_dir);
fprintf(fp,'/bin/rm -f %s\n',nms{f});
fclose(fp);
unix('chmod a+x ~/misc/ded-rm');
