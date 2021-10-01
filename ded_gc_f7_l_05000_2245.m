
%rsync -uavp asahi:gc/f7/l/05000/2245/b/b-00*0.hdf5 ~/gc/f7/l/05000/2245/b
nm='gc/f7/l/05000/2245';

p.flipy= 1
p.flipz= 0
p.aa= 1
p.png= 1


p.rgy= [0 0.4800]
p.rgz= [0 2]
p.gdx= 0.5000
p.gdy= 0.5000
p.box= 1
p.bottom= 1
p.isoval= 0.0200
p.isocap= 'all'
p.rg= [-Inf Inf]
p.fminx= Inf
p.fmaxx= -Inf
p.fminy= Inf
p.fmaxy= -Inf
p.fminz= Inf
p.fmaxz= -Inf
p.dq= 0
p.y0= 0
p.z0= 0
p.top= 0
p.bh= 0.1000
p.pcx= [2 3 1]
p.incremental= 0
p.col= {'red'  'blue'}
p.maxfn= Inf
p.nrot= 1
p.zmax= 0.7000
p.trackfront= 0
p.rd= 'opengl'
p.tsn= NaN
p.za= 0
p.ctar= 1
p.cont= 1

p.text ='';
p.lw=2;
           
p.sz = 1920*[1 9/16]; % 16x10

p.aa = 4;
p.nmbnm='gc/f7/l/05000/2245-b-1';

%p.x0= 0
%p.rgx= [0 10.5000]
%p.cva= 32

p.x0 = 8.32;
p.rgx = [-8 0.5];
p.transparent=0;
p.cva=37.1;
ded_3d(nm,p,1);
