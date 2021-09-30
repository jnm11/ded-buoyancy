
opt={'--Ty SinCos --Tz SinCos','--Ty Fourier --Tz Fourier'};
for k=1:2
  nm = sprintf('pm/sawtooth/%02u',k);  
  fn=[ded_dedalus_data_dir '/' nm '/force/diffdom.hdf5'];
  if ~isfile(fn)
    unix(sprintf('rm -rf ~/%s/*',nm))
    unix(sprintf('mpiexec -n 4 ded_gc.py --preset pmf8 --Nx 32 -L 1 -H 1 -W 1 -T 0.01 %s %s',opt{k},nm));
  end
  p  = ded_read_param(nm);
  c  = ded_coord(nm);
  a  = ded_read_hdf(fn);
  figure(2*k-1);clf;
  subplot(2,2,1);mesh(c.z,c.y,a.ddfy(:,:,1));axis('tight');
  subplot(2,2,2);mesh(c.z,c.y,a.ddfz(:,:,1));axis('tight');
  subplot(2,2,3);mesh(c.z,c.y,a.ddf(:,:,1));axis('tight');
  subplot(2,2,4);mesh(c.x,c.y,squeeze(a.ddfx(1,:,:)));axis('tight');
  figure(2*k);clf;
  %  fy=reshape(permute(a.ddfy,[2 1 3]),[c.Ny c.Nx*c.Nz]);
  %fz=reshape(        a.ddfz         ,[c.Nz c.Nx*c.Ny]);
  fy=a.ddfy(1,:,1)';
  fz=a.ddfz(:,1,1);
  subplot(3,1,1);plot(c.y,fy,'-s',c.y,c.y);xlabel('y');axis('tight');
  subplot(3,1,2);plot(c.z,fz,'-s',c.z,c.z);xlabel('z');axis('tight');
  subplot(3,1,3);plot(c.z,fz-fy,'-s',c.z,0*c.z);xlabel('z');axis('tight');
end
