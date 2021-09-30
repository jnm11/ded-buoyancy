% Process 2d lock exchange simulations
%2d lock exchange dns
% This part on hamilton
dd=[ded_dedalus_data_dir '/results'];
cd('~/');
fns=cellstr_ls('gc/le/[ab]/*/*/*/status',[],'dir');
%ded_make_param(fns);
minx=24;
xrg=[29 46];
%ded_gc_profiles_slice(fns,minx,xrg,'y',true);
for j=1:length(fns)
  nm=fns{j};
  a=ded_gc_traj(nm,xrg);
  ded_gc_avg_le(nm,'y','y',a.Xf,a.Uf,minx,[a.t1 a.t2]);
end
ded_gc_fit_erf_profiles(fns,[],'ry');


return;

dd=[ded_dedalus_data_dir '/results'];
a=load('~/Dropbox/GFD-2019-Gravity-Currents/data/new/dns-le-2d-2.mat');a=a.a;
b=load('~/Dropbox/GFD-2019-Gravity-Currents/data/dns-le-2d-2.mat');b=b.a;
k=3;
ga.lt='-';ga.nm=cellstrremoveprefix(a(k).nm,'results/gc/f6/g/');
gb.lt='-';gb.nm=cellstrremoveprefix(b(k).nm,'results/gc/f6/g/');
dc;
figure;clf;[h lh uc]=groupplot(a(k).nRe,a(k).R,      a(k).typ2,ga);legend(lh,uc,'location','best');ylabel('J');
[h lh uc]=groupplot(b(k).nRe,b(k).R,      b(k).typ2,gb);legend(lh,uc,'location','best');ylabel('J');

