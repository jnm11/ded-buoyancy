cd('~/data/dedalus/results');
%fn1=cellstr_ls('gc/emle/01[79]/*-steady-cheb.mat');
%fn2=cellstr_ls('gc/ccle/024/*-steady-cheb.mat');
%fn3=cellstr_ls('gc/ccle/046/*-steady-cheb.mat');
%fns=cat(1,fn1,fn2,fn3);   
fns={'gc/emle/017/16-steady-cheb.mat'       
     'gc/emle/017/20-steady-cheb.mat'       
     'gc/emle/017/24-steady-cheb.mat'       
     'gc/emle/017/28-steady-cheb.mat'       
     'gc/emle/017/32-steady-cheb.mat'       
     'gc/emle/019/16-steady-cheb.mat'       
     'gc/emle/019/20-steady-cheb.mat'       
     'gc/emle/019/24-steady-cheb.mat'       
     'gc/emle/019/28-steady-cheb.mat'       
     'gc/emle/019/32-steady-cheb.mat'       
     'gc/ccle/024/18-steady-cheb.mat'       
     'gc/ccle/024/22-steady-cheb.mat'       
     'gc/ccle/024/26-steady-cheb.mat'       
     'gc/ccle/024/30-steady-cheb.mat'       
     'gc/ccle/046/16-steady-cheb.mat'       
     'gc/ccle/046/20-steady-cheb.mat'       
     'gc/ccle/046/24-steady-cheb.mat'       
     'gc/ccle/046/28-steady-cheb.mat'};

for j=1:length(fns)
  a=load(fns{j});
  a=a.a;
  a.nm=fileparts(fns{j});
  a.fnmat=[fns{j}(1:end-16) '-profile.mat'];
  ded_make_param(a.nm);
  %ded_gc_fit_erf_profiles(a,[],'y');
end
fns={'gc/emle/017/16-profile.mat'       
     'gc/emle/017/20-profile.mat'       
     'gc/emle/017/24-profile.mat'       
     'gc/emle/017/28-profile.mat'       
     'gc/emle/017/32-profile.mat'       
     'gc/emle/019/16-profile.mat'       
     'gc/emle/019/20-profile.mat'       
     'gc/emle/019/24-profile.mat'       
     'gc/emle/019/28-profile.mat'       
     'gc/emle/019/32-profile.mat'       
     'gc/ccle/024/18-profile.mat'       
     'gc/ccle/024/22-profile.mat'       
     'gc/ccle/024/26-profile.mat'       
     'gc/ccle/024/30-profile.mat'       
     'gc/ccle/046/16-profile.mat'       
     'gc/ccle/046/20-profile.mat'       
     'gc/ccle/046/24-profile.mat'       
     'gc/ccle/046/28-profile.mat'};

dd=[ded_dedalus_data_dir '/results/gc/ccle'];
xmin=-inf;
a=ded_gc_f7_g_classify(fns,-25,false);    save([dd '/dns-le-3d-2.mat'],'a');
a=ded_gc_f7_g_RJ(fns,xmin,false,'front'); save([dd '/dns-le-3d-1-f.mat'],'a');
a=ded_gc_f7_g_RJ(fns,xmin,false,'tail' ); save([dd '/dns-le-3d-1-t.mat'],'a');

if false

  a=ded_slice('gc/emle/019',{'u','b'}, 5,[2 inf],'y','slice');
  %a=ded_slice('gc/ccle/024',{'u','b'},10,[2 27], 'y','slice');
  
  a=ded_gc_traj('gc/emle/019',[8 32.5],1);plot(a.t,a.X,a.t,a.fx(a.t)); % 5                                                                     
  ded_gc_avg_le('gc/emle/019',{'y'},'y',a.fx,a.fu,5,[8 inf]);
  
  a=ded_gc_traj('gc/ccle/024',[5  46],1);plot(a.t,a.X,a.t,a.fx(a.t)); % 10
  ded_gc_avg_le('gc/ccle/024',{'y'},'y',a.fx,a.fu,10,[8 47]);
  
  a=ded_gc_traj('gc/ccle/046',[8 32.5],1);plot(a.t,a.X,a.t,a.fx(a.t)); % 10
  ded_gc_avg_le('gc/ccle/046',{'y'},'y',a.fx,a.fu,5,[8 32.5]);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 3d lock release
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  rm={'avrg','bV','clipB','clipS','clipu','dnoise','dtavrg','dtb','dtcheck','dtforce','dtfpid','dtgc','dtjWx','dtjWy','dtjWz','dtjavrg','dtjb','dtjcheck','dtjp','dtjr','dtjs','dtju','dtjv','dtjw','dtjxyz','dtjy','dtjyz','dtleft','dtmomb','dtnoise','dtp','dtpm','dtright','dtslice','dtu','dtv','dtw','dtx','dtxy','dtxyz','dtxz','dty','dtyz','dtz','fpid','fpidb','fpids','fpidu','fpidv','fpidw','gmax','gmin','inlet','m1','m2','noise','noiseL','noiseT','noised','num','pmsl','radius','series','xn1','xn2','xn3','xn4','U','Umax','Umin','V','VMEM','Wx','Wy','Wz','alpha','PIDD','PIDDD','PIDG','PIDI','PIDIT','PIDP','PIDS1','PIDS2','PIDST','PIDT','PIDX','time'};
  ded_cmp_sims({'gc/ccle/046','gc/ccle/024','gc/emle/017','gc/emle/019'},rm);
  %        name      L   Nx   Ny   Nz       Peb        Re       Ty                       X
  % gc/ccle/024   30.0 2400  200  120    8000.0    8000.0  Fourier   30.00    2.50    1.00 
  % gc/ccle/046   21.0 4200   72  300   10000.0   10000.0   SinCos   21.00    0.36    1.00 
  % gc/emle/017   21.0 2800   48  216   10000.0   10000.0   SinCos   21.00    0.36    1.00 
  % gc/emle/019   21.0 4200   72  300   20000.0   10000.0   SinCos   21.00    0.36    1.00 
end
