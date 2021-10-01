


a=ded_read_g('gc/411','y');
h=ded_read_param('gc/411');

clf;
imagesc(a.x,a.z,mean(a.w,3));set(gca,'ydir','normal');
colorbar;

b=ded_read_g('gc/411','force');
b.wb=squeeze(mean(b.wb,2));
b.wu=squeeze(mean(b.wu,2));
b.fu=squeeze(mean(b.fu,2));
b.fb=squeeze(mean(b.fb,2));

hold('on');
[Cx Cy]=longest_contours(b.x,b.z,b.wb,max(b.wb(:))*.99,2);
for j=1:length(Cx);plot(Cx{j},Cy{j});end






[Cx Cy]=longest_contours(b.x,b.z,b.fb,max(b.fb(:))*.99,2);
for j=1:length(Cx);plot(Cx{j},Cy{j});end

hold('on');
[Cx Cy]=longest_contours(b.x,b.z,b.wu,max(b.wu(:))*.99,2);
for j=1:length(Cx);plot(Cx{j},Cy{j});end
[Cx Cy]=longest_contours(b.x,b.z,b.wu,max(b.wu(:))*.01,2);
for j=1:length(Cx);plot(Cx{j},Cy{j});end

plot(a.x,squeeze(a.b(1,:,:))/h.W,[0 24],[1 1],'--');

figure(2);
clf;
c=ded_read_g('gc/411','avrg');
zz=linspace(0,1,c.nz*2);
b=interp1(c.z,c.ab,zz,'linear','extrap')/c.t;
u=interp1(c.z,c.au,zz,'linear','extrap')/c.t;
w=interp1(c.z,c.aw,zz,'linear','extrap')/c.t;
w([1 end],:)=0;
psi=helmholtz_fft(u',w',c.x(2)-c.x(1),zz(2)-zz(1))';

imagesc(c.x,zz,b);set(gca,'ydir','normal');

hold('on');
% $$$ dx=20;dz=5;
% $$$ quiver(c.x(1:dx:end),c.z(1:dz:end),c.au(1:dz:end,1:dx:end),c.aw(1:dz:end,1:dx:end));
C=contourc(c.x,zz,psi,50);
hh=plot_contours(C);
set(hh,'color',[0 0 0]);



[Cx Cy]=longest_contours(c.x,c.z,c.au,0,1);
plot(Cx,Cy,'b');

[Cx Cy]=longest_contours(c.x,c.z,c.aw,0,4);
for j=1:length(Cx);plot(Cx{j},Cy{j},'k--');end


cc=0.05:.1:0.95;
for k=1:length(cc)
[Cx Cy]=longest_contours(c.x,zz,b,cc(k),2);
for j=1:length(Cx);plot(Cx{j},Cy{j},'w');end
end
get(gcf,'renderer','opengl');
set(gca,'position',[0 0 1 1]);
set(gcf,'paperposition',8*[0 0 12 1]);
set(gcf,'units','pixels','position',[10 500 2410 200]);
print('-depsc2','-r1200','/home/vzfv57/Dropbox/Jim-Claudia-GC/a.eps');


