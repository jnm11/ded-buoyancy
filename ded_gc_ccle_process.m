function ded_gc_ccle_process
%rsync -vap hamilton:mat-ccle-* /home/vzfv57/Dropbox/Jim-Claudia-GC/
dd=ded_dedalus_data_dir;
fn1=cellstr_ls([dd '/gc/emle/*']);
fn2=cellstr_ls([dd '/gc/ccle/*']);
fns={fn1{:},fn2{:}};
fns={[dd '/gc/emle/017']};
fns={[dd '/gc/ccle/024']};
fns={[dd '/gc/emle/004']};
fns={'ccle/024','emle/004','emle/005','emle/006','emle/007',...
     'emle/008','emle/010','emle/009',...
     'emle/011','emle/012','emle/013','emle/014',...
     'emle/015','emle/016','emle/017','emle/018','emle/019'};
fns={'emle/019'};
fns={'emle/017','ccle/024','ccle/024'};
tfinal=33.5;
tfinal=25.5;
for j=1:length(fns)
  [nm2 nm1]=fileparts(fns{j});
  dd=['~/Dropbox/Jim-Claudia-GC/mat/' nm2 '/' nm1];
  if ~isdir(dd); mkdir(dd); end
  nm=['gc/' nm2 '/' nm1];
  p=ded_read_param(nm);
  switch(p.L)
    case 21
      ttraj=5;
      tstart=16:4:32;
    case 30
      ttraj=5;
      tstart=18:4:34;
  end
  a=ded_gc_traj(nm,dd,[ttraj tfinal],1);
  for t=tstart
    snm=sprintf('%4.1f-%4.1f',t,tfinal);
    ded_gc_steady(nm,dd,[t tfinal],a.fx,a.fu, 'y',[-20 5],snm);
  end
  ded_gc_yint(nm,dd);
end
