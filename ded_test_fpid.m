mpiexec -n  10 ded_gc.py --dtavrg 0 --preset gcf6 -c 1 --AA 1 --dtforce 1 --Nx 1024 --Ny 1 --Re 100 --fpidb True --fpidu True --fpidw True --fpid 0.1 --force 0 gc/fpid/01
mpiexec -n  10 ded_gc.py -r -p --force 100  gc/fpid/01
 
rm -rf $DEDALUS_DATA/gc/fpid/01 ; ded_gc.py --preset gcf6 --Re 100 --Nx 1024 --Ny 1 --fpid 10 --fpidu True --fpidw True --fpidb True --force 100 gc/fpid/01

mpiexec -n 4 ded_gc.py -c 1 --dtforce 1 -r -p --Wx 0.5 --Wz 0.5  -T 2 --x0 1 --x1 2 --x2 4 --x3 5 --x4 1 --x5 2 --x6 3 --x7  5 gc/fpid/01

nm='gc/fpid/01';

p=ded_read_param(nm);
f=ded_gc_plot_forcing(nm);
nmf=[ded_dedalus_data_dir '/' nm '/final/final_s1.hdf5'];
%
%for j=28:28
j=10;
nmf=sprintf('%s/%s/checkpoint/checkpoint_s%u.hdf5',ded_dedalus_data_dir,nm,j);
  s=ded_read_state(nmf,p);
  f=ded_read_g(nm,'force');
  
  s.x=(0.5:p.Nx-0.5)*p.L/p.Nx;
  subplot(3,2,1);
  mesh(s.x,s.z,s.fpidu);title('fpidu');axis('tight');
  subplot(3,2,3);
  mesh(s.x,s.z,s.fpidw);title('fpidw');axis('tight');
  subplot(3,2,5);
  mesh(s.x,s.z,s.fpidb);title('fpidb');axis('tight');
  subplot(3,2,2);
  mesh(s.x,s.z,f.wu.*(s.u-f.fu));title('u-fu');axis('tight');
  subplot(3,2,4);
  mesh(s.x,s.z,f.wu.*(s.w-f.fw));title('w-fw');axis('tight');
  subplot(3,2,6);
  mesh(s.x,s.z,f.wb.*(s.b-f.fb));title('b-fb');axis('tight');
  %end






f=ded_zgrid(f,nz,nm,nmd,nmdd,nmi,H)



mesh(f.wu.*(f.fu-s.u));

a=ded_read_g(nm,'fpid');mesh(a.fpidu(:,:,end));


