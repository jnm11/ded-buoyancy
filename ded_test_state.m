

dd=ded_dedalus_data_dir;

fn1=[dd '/gc/test/004/final/state-00008.hdf5'];
c=ded_coord('gc/test/004');
a=ded_read_hdf(fn1);
b1=ded_read_hdf([dd '/gc/test/004/state1.hdf5']);
b2=ded_read_hdf([dd '/gc/test/004/state2.hdf5']);
b3=ded_read_hdf([dd '/gc/test/004/state3.hdf5']);

d=struct_cmp(a,b1);
d=struct_cmp(a,b2);
d=struct_cmp(a,b3);

plot(c.x,a.bz(1,:),c.x,b3.bz(1,:))


plot(c.x,b3.bz(1,:))


nm='gc/test/004';
