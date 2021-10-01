
DD=[ded_dedalus_data_dir '/'];
fns=cellstr_ls([DD 'gc/ccle/*'],'ls -d');
fns=cellstrremoveprefix(fns,DD);
do='~/Dropbox/Jim-Claudia-GC/mat-files/';

for j=9:length(fns)
  nm=fns{j};
  s=ded_read_stats(nm);
  p=ded_read_param(nm);
  y=ded_read_g(nm,'y');
  l=ded_read_g(nm,'l');
  r=ded_read_g(nm,'r');
  gc=ded_read_g(nm,'gc');
  yz=ded_read_g(nm,'yz');
  xyz=ded_read_g(nm,'xyz');
  momb=ded_read_g(nm,'momb');
  a=ded_read_g(nm,'avrg');
  dd=[do nm '/'];
  unix(sprintf('mkdir -p %s',dd));
  save([dd 'stats.mat'],'s');
  save([dd 'param.mat'],'p');
  save([dd 'y.mat'],'y');
  save([dd 'l.mat'],'l');
  save([dd 'r.mat'],'r');
  save([dd 'yz.mat'],'yz');
  save([dd 'gc.mat'],'gc');
  save([dd 'xyz.mat'],'xyz');
  save([dd 'momb.mat'],'momb');
  save([dd 'avrg.mat'],'a');
end


%rsync -uavp ~/Dropbox/Jim-Claudia-GC/mat-files asahi:Dropbox/Jim-Claudia-GC
  
  
  
