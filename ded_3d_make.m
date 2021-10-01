p=struct('cva',30.0721);
ded_3d('gc/qgcf3','~/films/dedalus/qgcf3',p);


nm='pgc/f6/069';
nmo='~/films/dedalus/pgc-f6-069';
p.aa=1;
p.sz=[1024 512];
p.rgx=[0 48];
ded_3d(nm,nmo,p);

nm='ppm/005';
nmo='~/films/dedalus/ppm-005';
p.aa=1;
p.sz=[576 768];
p.rgx=[0 30];
ded_3d(nm,nmo,p);

 
