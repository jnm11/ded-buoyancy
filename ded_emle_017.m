nm='gc/emle/017';
p=ded_3d_param(nm);;
p.isoval = 0.02;


p.aa=4;

p.sz = 1920*[1 9/16]/2; % 16x10
p.sz = 1024*[1 9/16]; % 16x10

%p.aa = 1;
%p.reduce=1;

t=load('~/Dropbox/Jim-Claudia-GC/mat/emle/017/traj.mat');
p.nmbnm='gc/emle/017-1';

p.X  = @(t) interp1(t.a.t,t.a.X,t,'linear');
p.X  = t.a.fx;
p.rgx = [-5 1];
p.cva = 35;
ded_3d(nm,p);

p.nmbnm='gc/emle/017-2';
p.rgx = [-40 8];
p.cva = 33;
ded_3d(nm,p);
