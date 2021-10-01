function ded_aliasing2(nm);
p=ded_read_param(nm);

c=ded_coord(nm);
dx=c.dx;
dz=c.dz;

fns=ded_get_fn(nm,'final',[],'state');
a=ded_read_hdf(fns{end});
pr=ded_read_parity(nm);
wb=ded_read_g(nm,'force','wb',[],'wb');
if c.dim==3
  bx=reshape(permute(a.b,[3 2 1]),c.Nx,c.Ny*c.Nz);
  bz=reshape(a.b,c.Nz,c.Nx*c.Ny);
else
  bz=a.b;
  bx=a.b';
end
N=64;
[fx kx]=pmtm(bx, [], N, 1/dx);  
[fz kz]=pmtm(bz, [], N, 1/dx);  

f=findmax(fx(end,:));

wb=max(wb.wb);
x4=min(c.Jx(wb<max(wb)/100));

figure(1);clf;
subplot(2,1,1);loglog(kx,mean(fx,2),kx,fx(:,f));
title('b x spectrum')
subplot(2,1,2);loglog(kz,mean(fz,2));
title('b z spectrum')
fl=find(c.x<x4);
figure(2);clf;
subplot(3,1,1);plot(c.x(fl),bx(fl,f)-mean(bx(fl,f)),'s-');title('b worst')
subplot(3,1,2);plot(c.x(fl),bx(fl,1)-1,'s-');title('b bottom')
subplot(3,1,3);plot(c.x(end-39:end),bx(end-39:end,1),'s-');title('b bottom')
figure(3);clf;

if c.dim==3
  subplot(2,1,1);
  plot(c.z,mean(sum(abs(diff(a.b,1,3)),3),2));
  title('int |db|')
  xlabel('z');
  subplot(2,1,2);
  plot(c.x,squeeze(mean(sum(abs(diff(a.b,1,1)),1),2)));
  xlabel('x');
else
  subplot(2,1,1);
  plot(c.z,sum(abs(diff(a.b,1,2)),2));
  title('int |db|')
  xlabel('z');
  subplot(2,1,2);
  plot(c.x,sum(abs(diff(a.b,1,1)),1));
  xlabel('x');
end

fb=[];wb=[];
fns=ded_get_fn(nm,'force',[],'wb');
if ~isempty(fns) 
  wb=ded_read_hdf(fns{end});
  [wfx wkx]=pmtm(wb.wb(1,:),[],N,1/dx);  
  figure(4);clf;
  loglog(wkx,wfx);hold('on');
  disp(sum(abs(diff(wb.wb(1,:))))-2);
  title('Buoyancy forcing weight function psd');
  xlabel('kx');
  xlabel('Power');
end

fns=ded_get_fn(nm,'force',[],'fb');
if ~isempty(fns) & false
  fb=ded_read_hdf(fns{end});
  [ffx fkx]=pmtm(fb.fb(1,:),[],N,1/dx);
  [ffz fkz]=pmtm(fb.fb(:,1),[],16,1/dz);
  figure(4);clf;
  subplot(2,1,1);loglog(fkx,ffx);
  title('Buoyancy forcing');
  xlabel('x');
  subplot(2,1,2);loglog(fkz,ffz);
  disp([sum(abs(diff(fb.fb(:,1))))-1 sum(abs(diff(fb.fb(1,:))))-1]);
  xlabel('z');
  keyboard
end

if isfield(a,'db') & ~isempty(wb)
  figure(5);clf;
  db=ichebf2c(a.db,1);db=ichebc2f(db(1:c.Nz),1); 
  imagesc((a.b(:,1:end/2)-db).*wb.wb(:,1:end/2)); 
  disp(max(max(abs((a.b(:,1:end/2)-db).*wb.wb(:,1:end/2)))));
end

%imagesc(log(abs(eps+pr_highpass(a.b,c.dx,2,c.dx))));   

if ~isempty(fb) & ~isempty(wb)
  figure(6);clf;
  e=(a.b-fb.fb).*wb.wb;
  imagesc(e(:,1:end/2));
  disp(max(max(abs(e))));
  f=findmax(abs(e(:,1)));
  figure(6);
  subplot(3,1,1);plot(c.x(1:10),e(end,1:10),'s-');title('b deviations from forcing top');
  subplot(3,1,2);plot(c.x(1:10),e(f,1:10),'s-');title('b deviations from forcing worst');
  subplot(3,1,3);plot(c.x(1:10),e(1,1:10),'s-');title('b deviations from forcing bottom');
end
%ded_cmp_sims('gc/f7/test/*');
