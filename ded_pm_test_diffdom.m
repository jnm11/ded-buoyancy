function ded_pm_test_diffdom(nm)
nm='pm/dd/01';
p=ded_read_param(nm);
c=ded_coord(nm);
dAA=c.dAy*c.dAz;
dA=c.dJy*c.dJz;
fns=cellstr_ls([ded_dedalus_data_dir '/' nm '/force/*.hdf5']);
for j=1:length(fns)
  [dd fn]=fileparts(fns{j});
  a.(fn)=ded_read_hdf(fns{j});
end

%subplot(2,2,1);plot(c.Ax,a.ddfx,'s-'); axis('tight');title(sprintf('I=%5.3f %5.1e',sum(a.ddfx)*c.dAx,min(a.ddfx(:))/max(a.ddfx(:))));
%subplot(2,2,2);mesh(c.Ay,c.Az,a.ddf);  axis('tight');title(sprintf('I=%5.3f %5.1e',sum(a.ddf(:))*dA,min(a.ddf(:) )/max(a.ddf(:) )));
%subplot(2,2,3);mesh(c.Ay,c.Az,a.ddfy); axis('tight');title(sprintf('I=%5.3f %5.1e',sum(a.ddfy(:))*dA,min(a.ddfy(:) )/max(a.ddfy(:) )));
%subplot(2,2,4);mesh(c.Ay,c.Az,a.ddfz); axis('tight');title(sprintf('I=%5.3f %5.1e',sum(a.ddfz(:))*dA,min(a.ddfz(:) )/max(a.ddfz(:) )));

figure(1)
subplot(3,2,1);mesh(c.Ay,c.Az,a.fyz.nry);
subplot(3,2,2);mesh(c.Ay,c.Az,a.fyz.nrz);
subplot(3,2,3);mesh(c.Ay,c.Az,a.fyz.nry.*a.fyz.rr);
subplot(3,2,4);mesh(c.Ay,c.Az,a.fyz.nrz.*a.fyz.rr);
subplot(3,2,5);mesh(c.Ay,c.Az,a.fyz.nrz.^2+a.fyz.nry.^2);
subplot(3,2,6);mesh(c.Ay,c.Az,a.fyz.wdrn);


fns=ded_get_fn(nm,'a');
fnyz=ded_get_fn(nm,'ayz');
a=ded_read_hdf(fns{end});
yz=ded_read_hdf(fnyz{end});
div=(a.dudx+a.dvdy+a.dwdz)/a.dt;
disp(sum(div(:)));


%calc_diff=1
%[a b]=ded_diff_int(a,p,dx,calc_diff,calc_int,calc_NS,fd);
 
dudn=a.dudx.*b.ddfx+a.dudy.*b.ddfy+a.dudz.*b.ddfz;
dvdn=a.dvdx.*b.ddfx+a.dvdy.*b.ddfy+a.dvdz.*b.ddfz;
dwdn=a.dwdx.*b.ddfx+a.dwdy.*b.ddfy+a.dwdz.*b.ddfz;

w=1/max(abs(b.ddfx(:)));
Nx=findmin(abs(c.x-p.L/2));
Ny=findmin(abs(c.y));
Nz=findmin(abs(c.z));
figure(1);
subplot(3,1,1);h=plot(c.x,w*(squeeze(b.ddfx(Nz,Ny,:))),'s-',c.x,squeeze(b.ddf(Nz,Ny,:)),'s-');
axis('tight');legend(h,{'fx','f'},'location','best');xlabel('x');
subplot(3,1,2);h=plot(c.y,w*(squeeze(b.ddfy(Nz,:,Nx))),'s-',c.y,squeeze(b.ddf(Nz,:,Nx)),'s-');
axis('tight');legend(h,{'fy','f'},'location','best');xlabel('y');
subplot(3,1,3);h=plot(c.z,w*(squeeze(b.ddfz(:,Ny,Nx))),'s-',c.z,squeeze(b.ddf(:,Ny,Nx)),'s-');
axis('tight');legend(h,{'fz','f'},'location','best');xlabel('z');

figure(2);
subplot(2,2,1);mesh(c.y,c.z,dvdn(:,:, Nx));axis('tight');title('dvdn');
subplot(2,2,2);mesh(c.y,c.z,dwdn(:,:, Nx));axis('tight');title('dwdn');
subplot(2,2,3);mesh(c.y,c.z,dvdn(:,:,end));axis('tight');
subplot(2,2,4);mesh(c.y,c.z,dwdn(:,:,end));axis('tight');

figure(3);
subplot(2,2,2);mesh(c.x,c.z,squeeze(b.ddfx(:,Ny, :)));axis('tight');
subplot(2,2,3);mesh(c.y,c.z,squeeze(b.ddfy(:, :, 1)));axis('tight');
subplot(2,2,4);mesh(c.y,c.z,squeeze(b.ddfz(:, :, 1)));axis('tight');


figure(4);clf;
subplot(2,2,1);mesh(c.y,c.z,squeeze(b.ddf( :,:,1)));axis('tight');
subplot(2,2,2);mesh(c.x,c.y,squeeze(b.ddf(Nz,:,:)));axis('tight');
subplot(2,2,3);mesh(c.x,c.z,squeeze(b.ddf(:,Ny,:)));axis('tight');


figure(5);clf;
subplot(3,1,1);plot(squeeze(cat(3,b.ddfx(Nz,Ny,:),-flip(b.ddfx(Nz,Ny,:)),b.ddfx(Nz,Ny,:))),'s-')
subplot(3,1,2);plot(squeeze(cat(2,b.ddfy( Nz,:,1),-flip(b.ddfy(Nz, :,1)),b.ddfy(Nz, :,1))),'s-')
subplot(3,1,3);plot(squeeze(cat(1,b.ddfz( :,Ny,1),-flip(b.ddfz( :,Ny,1)),b.ddfz( :,Ny,1))),'s-')

figure(6);clf;
subplot(3,1,1);plot(c.x,reshape(a.u([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-');ylabel('u');
subplot(3,1,2);plot(c.x,reshape(a.v([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-');ylabel('v');
subplot(3,1,3);plot(c.x,reshape(a.w([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-');ylabel('w');

figure(7);clf;
subplot(4,1,1);plot(c.x,reshape(a.dudx([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-');ylabel('u');
subplot(4,1,2);plot(c.x,reshape(a.dvdx([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-');ylabel('v');
subplot(4,1,3);plot(c.x,reshape(a.dwdx([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-');ylabel('w');
subplot(4,1,4);plot(c.x,reshape(   div([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-');ylabel('div');

figure(8);clf;
subplot(2,2,1);mesh(c.y,c.z,div(:,:,1*Nx/2));axis('tight');title('div 1*x/4');
subplot(2,2,2);mesh(c.y,c.z,div(:,:,2*Nx/2));axis('tight');title('div 2*x/4');
subplot(2,2,3);mesh(c.y,c.z,div(:,:,3*Nx/2));axis('tight');title('div 3*x/4');
subplot(2,2,4);mesh(c.y,c.z,div(:,:,end));axis('tight');title('div end');

figure(9);clf;
[zz yy]=ndgrid(c.z,c.y);
r=sqrt(yy.^2+zz.^2);
[w f xg]=jgrid(r(:)',reshape(div,[c.Ny*c.Nz,c.Nx])',c.dy,'quadratic');
divr=(f./w);
mesh(xg,c.x,divr);
xlabel('r');
ylabel('x');

figure(10);clf;
div1=squeeze(sum(sum(   b.ddf .*div*c.dy*c.dz,1),2));
div2=squeeze(sum(sum((1-b.ddf).*div*c.dy*c.dz,1),2));
subplot(2,1,1);plot(c.x,div1);axis('tight');title('inside'); axis([0 p.L -1e-1 1e-1]);
subplot(2,1,2);plot(c.x,div2);axis('tight');title('outside');axis([0 p.L -1e-1 1e-1]);


%subplot(2,1,2);plot(c.x,reshape(a.dvdx([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-')
%subplot(4,1,3);plot(c.x,reshape(a.dwdx([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-')
%subplot(4,1,4);plot(c.x,reshape(   div([1 Nz end],[1 Ny end],:)/a.dt,[9 c.Nx]),'s-')


%[z y]=ndgrid(c.z,c.y);
%r=max(eps,sqrt(y.^2+z.^2));
%Ix=-squeeze(sum(b.ddfx,3))*c.dx;
%Ir=-squeeze(sum(sum(b.ddfy.*y./r+b.ddfz.*z./r,1),2))*c.dy*c.dz/(2*pi);
%ddfr=sqrt(b.ddfy.^2+b.ddfz.^2);
%Is=squeeze(sum(sum(ddfr,1),2))*c.dy*c.dz/(2*pi);


keyboard


yz.div=yz.dudx+yz.dvdy+yz.dwdz;

plot(c.x,yz.div(:));axis([0 p.L -0.1 0.5]);

