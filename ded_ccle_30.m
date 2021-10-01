nm='gc/ccle/030';
p=ded_read_param(nm);
b=ded_read_hdf(['~/' nm '/b/b-00000.hdf5']);
y=ded_read_hdf(['~/' nm '/y/y-00000.hdf5']);

plot(y.x,y.b,'-s');
