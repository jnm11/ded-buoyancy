nm='gc/emle/014'
nm='gc/test/01'


p=ded_read_param(nm);
a=ded_read_hdf([ded_dedalus_data_dir '/' nm '/dimple.hdf5']);
contour(a.y,a.z,a.xb,50);
xlabel('y');
ylabel('z');
set(gca,'dataaspectratio',[ 1 1 1]);



% $$$ kx=ded_read_hdf([ded_dedalus_data_dir '/' nm '/kx.hdf5']);
% $$$ ky=ded_read_hdf([ded_dedalus_data_dir '/' nm '/ky.hdf5']);
% $$$ kz=ded_read_hdf([ded_dedalus_data_dir '/' nm '/kz.hdf5']);
% $$$ 
% $$$ nx=(0:p.Nx-1)';
% $$$ ny=(0:p.Ny-1)';
% $$$ nz=(0:p.Nz-1)';
% $$$ 
% $$$ subplot(3,1,1);
% $$$ plot(nx,kx.kx-nx*pi/p.L);
% $$$ subplot(3,1,2);
% $$$ plot(ny,ky.ky-ny*pi/p.W);
% $$$ subplot(3,1,3);
% $$$ plot(ny,kz.kz-nz);






figure(1);
mesh(a.y,a.z,a.xb);
xlabel('y');
ylabel('z');


figure(2);
pz=mean(abs(fft(a.xb,[],1)).^2,2);
py=mean(abs(fft(a.xb,[],2)).^2,1);
h=loglog((1:p.Ny)/p.Ny,py,(1:p.Nz)/p.Nz,pz);
legend(h,{'y','z'});
axis([0 0.5 1e-20 1e1]);


disp(std(a.xb(:)))

a=ded_read_hdf([ded_dedalus_data_dir '/' nm '/b/b-00000.hdf5']);
p=ded_read_param(nm);

[xb2 e]=solve_first(a.x,a.b-0.5,3);
