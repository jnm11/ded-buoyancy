dc;
nm='pm/030';
nm='pm/test/01';
cd([ded_dedalus_data_dir '/' nm]);
p=ded_read_param(nm);
Nxx=round(p.Nx*p.AA);

p.Nxx=round(p.Nx*p.AA);
p.Nyy=round(p.Ny*p.AA);
p.Nzz=round(p.Nz*p.AA);
V=p.L*p.H*p.W;
N=p.Nx*p.Ny*p.Nz;
NN=p.Nxx*p.Nyy*p.Nzz;
dV=V/N;
dVV=V/NN;
dx=p.L/p.Nx;dxx=p.L/p.Nxx;
dy=p.W/p.Ny;dyy=p.W/p.Nyy;
dz=p.H/p.Nz;dzz=p.H/p.Nzz;
dA=dy*dz;
dAA=dyy*dzz;

wux=ded_read_hdf('wux.hdf5');
Iux=ded_read_hdf('Iux.hdf5');
xx=wux.x;
figure(1);
subplot(2,1,1);
W=[wux.wux1 wux.wux2 wux.wux3 wux.wux4 wux.wux5];
I=[Iux.Iux1 Iux.Iux2 Iux.Iux3 Iux.Iux4 Iux.Iux5];
h=plot(xx,W);  
axis([0 p.L -0.05 1.05]);
title('Weight functions');
legend(h,{'1 mask','2 wu','3 gaussian','4 switch','5 injection'});
subplot(2,1,2);
h=plot(xx,wux.wux1,filter_midpoint(xx),diff(wux.wux4));  
NW=sum(~isfinite(W));
NI=sum(~isfinite(W));
if NW>0
  disp(sprintf('%u non finite wux',NW));
end
if NI>0
  disp(sprintf('%u non finite Iux',NI));
end




figure(2);
h=plot(xx,I);  
axis([0 p.L 0 inf]);
title('Density functions');
h=plot(xx,I);
II=sum(I)*dVV*p.Nyy*p.Nzz;
disp(II-1);


wr=ded_read_hdf('nwdr.hdf5');


% $$$ figure(3);
% $$$ mesh(wr.y,wr.z,wr.nwdr);
% $$$ A=sum(wr.nwdr(:))*dyy*dzz;
% $$$ title(sprintf('A=%f',A));
% $$$ [yy zz]=ndgrid(wr.y,wr.z);
% $$$ q=atan2(yy,zz);
% $$$ dq=2*pi/p.Ny;
% $$$ [w f qg]=jgrid1(q(:)'/dq,wr.nwdr(:)',-p.Ny/2,p.Ny/2,'cubic');
% $$$ w(3:4)=w(3:4)+w(end-1:end);w(end-3:end-2)=w(end-3:end-2)+w(1:2);w([1:2 end-1:end])=[];
% $$$ f(3:4)=f(3:4)+f(end-1:end);f(end-3:end-2)=f(end-3:end-2)+f(1:2);f([1:2 end-1:end])=[];
% $$$ qg([1:2 end-1:end])=[];
% $$$ qg=qg*dq;
% $$$ figure(4);
% $$$ subplot(2,1,1);plot(qg,w);
% $$$ subplot(2,1,2);plot(qg,f);
% $$$ 



figure(5);
fnd=ded_get_fn(nm,'d');
fnu=ded_get_fn(nm,'u');
d=ded_read_hdf(fnd{end});

nw0 = interpft(interpft(wr.nwdr,p.Nz ,1),p.Ny, 2);
uu  = interpft(interpft(wr.nwdr,p.Nzz,1),p.Nyy,2);

AU=pi*p.radius^2*p.U;
xx1=wux.wux5; 
for j=1:length(fnu)
  u=ded_read_hdf(fnu{j});
  d=ded_read_hdf(fnd{j});
  [dd nmu]=fileparts(fnu{j});
  [dd nmd]=fileparts(fnd{j});
  qv1=squeeze(sum(sum(u.u.*repmat(nw0,[1 1 size(u.u,3)]))))*dA/AU;
  subplot(3,1,1);
  Iu=squeeze(sum(sum(u.u)))*dA/AU;
  Iuu=squeeze(sum(sum(u.u.^2)))*dA/AU/p.U;
  Id=squeeze(sum(sum(d.d)))*dAA/AU/p.U;
  maxd=max(Id);
  maxqv=max(qv1);
  maxu=max(max(Iu),max(Iuu));
  h=plot(u.x,Iu,u.x,Iuu,u.x,ones(size(u.x)));
  line([p.x0([1;1]),p.x1([1;1]),p.x2([1;1]),p.x3([1;1])],[0;maxu],'color',0.7*[1 1 1]);
  ylabel('int u');
  legend(h,{'u','uu','1'},'location','se','fontsize',6);
  title(sprintf('%s',nmu));axis('tight');
  subplot(3,1,2);
  h=plot(u.x,qv1,xx,xx1*maxqv,u.x,u.x*0);axis('tight');
  line([p.x0([1;1]),p.x1([1;1]),p.x2([1;1]),p.x3([1;1])],[0;maxqv],'color',0.7*[1 1 1]);
  ylabel('qv');
  legend(h,{'qv1','xx1','0'},'location','se','fontsize',6);
  subplot(3,1,3);
  h=plot(d.x,Id);axis([0 p.L 0 maxd]);
  line([p.x0([1;1]),p.x1([1;1]),p.x2([1;1]),p.x3([1;1])],[0;maxd],'color',0.7*[1 1 1]);
  ylabel('int d');
  title(sprintf('%s %6.4f',nmd,sum(d.d(:))));
  drawnow;
end

figure(9);
x=u.x;
y=u.y;
z=u.z;
ff=u.u;
mind=min(ff(:));maxd=max(ff(:));
for j=1:length(x)
  mesh(y,z,ff(:,:,j));
  title(sprintf('%5.2f',x(j)));
  axis([-inf inf -inf inf mind maxd]);
  drawnow;
end




figure(6);
subplot(2,1,1);plot(d.x,squeeze(sum(sum(d.d)))*dAA/AU);ylabel('d');
subplot(2,1,2);plot(u.x,squeeze(sum(sum(u.u)))*dAA/AU,u.x,ones(size(u.x)));ylabel('u');

if 0
  figure(7);
  fd=ded_read_hdf('fd.hdf5');
  fd0=ded_read_hdf('fd0.hdf5');
  f=ded_read_hdf('force/force_s1.hdf5');
  dA=(f.y(2)-f.y(1))*(f.z(2)-f.z(1));
  Ifu=dA*squeeze(sum(sum(f.fu)))/AU;
  Ifd1=dA*squeeze(sum(sum(f.fd)))/AU;
  Ifd2=dAA*squeeze(sum(sum(fd.fd)))/AU;
  Ifd0=dAA*squeeze(sum(sum(fd0.fd)))/AU;
  h=plot(f.x,squeeze(max(max(f.wu))),f.x,Ifu,f.x,Ifd1,fd.x,Ifd2,fd.x,Ifd0,xx,x1)
  legend(h,{'wu','Ifu','Ifd1','Ifd2','x1'});


  figure(8);
  plot( fd.x,dAA*squeeze(sum(sum( fd.fd))),fd0.x,squeeze(min(min( fd.fd))),fd0.x,squeeze(max(max( fd.fd))));
end



figure(10);
wr=ded_read_hdf('nwdr.hdf5');
j=-1;
k=0;
while 1
  j=j+1;
  fna=sprintf('qv1/qv1-%05i.hdf5',j);
  fnb=sprintf('qvx/qvx-%05i.hdf5',j);
  %fnc=sprintf('qv2/qv2-%05i.hdf5',j);
  %fnu=sprintf('uu/u-%05i.hdf5',j);
  if ~isfile(fna) | ~isfile(fna)
    if k==1 
      break;
    end 
    continue;
  end
  k=1;
  a=ded_read_hdf(fna);
  b=ded_read_hdf(fnb);
  %c=ded_read_hdf(fnc);
  %u=ded_read_hdf(fnu);
  %dAA=(u.y(2)-u.y(1))*(u.z(2)-u.z(1));
  %qv2=squeeze(dAA*sum(sum(u.u.*wr.nwdr)));p
  %plot(u.x,qv1+qv2);
  %qv1=a.qv1;
  subplot(2,1,1);plot(xx,a.qv1(:));ylabel('qv1');title(fna);
  subplot(2,1,2);plot(xx,b.qvx(:));ylabel('qvx');title(fnb);
  drawnow;
end

if 0
  wr=ded_read_hdf('nwdr.hdf5');
  j=0;
  k=0;
  j=0
  fnu=sprintf('fu/fu-%05i.hdf5',j);
  fnd=sprintf('fd/fd-%05i.hdf5',j);
  u=ded_read_hdf(fnu);
  d=ded_read_hdf(fnd);
  
  Iu=squeeze(sum(sum(u.u)))*dyy*dzz;
  Id=squeeze(sum(sum(d.d)))*dyy*dzz;
  subplot(2,1,1);plot(xx,Iu);
  subplot(2,1,2);plot(xx,Id);
  
  
  
  %u=ded_read_hdf(fnu);
  %dAA=(u.y(2)-u.y(1))*(u.z(2)-u.z(1));
  %qv2=squeeze(dAA*sum(sum(u.u.*wr.nwdr)));p
  %plot(u.x,qv1+qv2);
  %qv1=a.qv1;
% $$$ subplot(2,1,1);plot(a.qv1(:));ylabel('qv1');title(fna);
% $$$ subplot(2,1,2);plot(b.qvx(:));ylabel('qvx');title(fnb);
% $$$ drawnow;
% $$$ 
end
