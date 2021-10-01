% look at vertical momentum balance

num='128';
xrg=[6 24];
a=ded_gc_load_tint(num);
b=ded_gc_read_param(num);
dx=b.L/b.Nx;
dz=(b.H/2);
dudx = fft_diff(a.au,1,2)/dx;
dwdx = fft_diff(a.aw,1,2)/dx;
dudz = ichebdiff(a.au,1,1)/dz;
dwdz = ichebdiff(a.aw,1,1)/dz;

w=dudz-dudx;
imagesc(a.x,a.z,w);set(gca,'ydir','normal');




plot(dudx,-dwdz,'.');

fr=ded_read_hdf([ded_dedalus_data_dir '/gc/' num '/force/force_s1.hdf5']);

div=dudx+dwdz;
cdiv=ichebf2c(squeeze(fr.fd),1);
fd=ichebc2f(cdiv(1:64,:));
%fd=squeeze(mean(reshape(fd,[64,2,768]),2));
fd=fd(:,1:2:end);
imagesc(a.x,a.z,div-fd);set(gca,'ydir','normal');
max(max(abs(div-fd)))
e=div(:)-fd(:);

e=log10(abs(div-fd));

imagesc(a.x,a.z,log10(abs(div-fd)));set(gca,'ydir','normal','clim',[-3 -1]);


fdiv=

imagesc(a.x,a.z,div+fd);set(gca,'ydir','normal');


imagesc(fr.x,fr.z,squeeze(fr.fd));set(gca,'ydir','normal')


B=ichebint(a.ab,1)*dz;
p=a.ap+a.auu.^2/2;
p=p-p(end,end);


figure(1);clf;
h=plot(a.x,p(end,:),a.x,p(1,:));
legend(h,{'top','bottom'});




f=384;
plot(a.z,a.au(:,f),a.z,sqrt(a.auu(:,f)-a.au(:,f).^2));

plot(a.z,a.auw(:,f));



Lw=(dz*ichebint(fft_diff(a.au,2,2),1)/dx^2+ichebdiff(a.aw,1,1)/dz)/b.Re;
Ruw=dz*ichebint(fft_diff(a.auw,1,2)/dx) ;


for f=300:600
figure(2);clf;
subplot(2,1,1);
plot(a.z,p(1,f)-p(:,f),a.z,B(:,f)-B(1,f));%,a.z,a.aww(:,f),a.z,Ruw(:,f),a.z,Lw(:,f));
axis('tight');
title(sprintf('x=%.3f',a.x(f)));
subplot(2,1,2);
plot(a.z,p(1,f)-B(:,f)-p(:,f)+B(1,f));
axis('tight');
drawnow
pause;
end




clf;imagesc(a.x,a.z,a.uu);set(gca,'ydir','normal');

