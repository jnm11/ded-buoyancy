
nm='gc/ccle/023b';

%rsync hamilton:gc/f6/f/09 -vaP ~/gc/f6/f --exclude "b*"   --exclude "check*"   --exclude "final*"  
p=ded_3d_param(nm);;

p.png = 1;
p.rgy = [0 2.5];
p.rgz = [0 1];

p.gdx = 0.5;
p.gdy = 0.5;
p.box = 1;
p.bottom = 1;
p.isoval = 0.02;
p.isocap = 'all';
p.dq = 0;
p.x0 = 0;
p.y0 = 0;
p.z0 = 0;
p.top = 0;
p.bh = 0.1;
p.pcx = [2 3 1];
p.col = {'red'  'blue'};
p.maxfn = Inf;
p.trackfront = 0;
p.rd = 'opengl';
p.tsn = NaN;
p.za = 0;
p.ctar = 1;
p.cont = 1;
p.ftype = 'b';
p.text ='';
p.lw=2;
p.reduce=1;

p.sz = 1920*[1 9/16]; 

p.aa = 4;
p.nmbnm='gc/ccle/023-1';

p.rgx = [8 18];
p.cva = 38;
ded_3d(nm,p,1);

p.aa = 4;
p.nmbnm='gc/ccle/023-2';
p.rgx = [0 30];
p.cva = 33;
ded_3d(nm,p,1);

%rsync -vap hamilton:tmp/dedalus/gc/ccle/023-1/ ~/gc-ccle-023-1
%rsync -vap hamilton:tmp/dedalus/gc/ccle/023-2/ ~/gc-ccle-023-2
