function ded_aliasing(nm)
%ded_aliasing('gc/f6/f/20');
%ded_aliasing('gc/test/20');
%nm='gc/f6/f/20';

c=ded_coord(nm);
p=ded_read_param(nm);
pr=ded_read_parity(nm);
ffn=0;
ff={'fw','fu','fb','wb','wu','final','a'};
ff={'final'};
tol=1e-25;
dx=p.L/p.Nx;
dz=p.H/p.Nz;
dy=p.W/p.Ny;

for j=1:length(ff)
  switch(ff{j})
    case 'final'
      fn=ded_get_fn(nm,'final',[],'state');
      gg={'u','v','w','b','p'};
      gg={'u','v','w','b','p'};
    case {'fw','fu','fb','wb','wu'}
      fn=ded_get_fn(nm,'force',[],ff{j});
      gg={'fu','fv','fw','wu','wv','ww'};
    case 'a'
      gg={'u','v','w','b','p'};
      fn=ded_get_fn(nm,'a');
  end
  if isempty(fn)
    disp(['empty ' ff{j}])
    continue;
  end
  fn=fn{end};
  disp(sprintf('ded_aliasing: %s',fn));
  
  a=ded_read_hdf(fn);
  
  gg=intersect(gg,fieldnames(a));
  
  for k=1:length(gg)
    g=gg{k};
    x=ded_parity(a.(g),flip(pr.(g)));
    ffn=ffn+1;
    figure(ffn);
    clf;
    switch(c.dim)
      case 2
        %x=x(:,100:end-99);
        Nx=size(x,2);
        Nz=size(x,1);
        disp(sprintf('%s %u %u %u',g,Nx,Nz));
        w=ichebintw(Nz);
        
        [ffx kx]=pmtm(permute(x,[2 1]),[],256,1/dx);
        fx=ffx*w';
        ffz=squeeze(ichebf2c(x,1).^2);
        fz=mean(ffz,2);
        kz=linspace(0,1,Nz)';
        fff=find(fx>tol & kx>0);fx=fx(fff); kx=kx(fff);
        fff=find(fz>tol & kz>0);fz=fz(fff); kz=kz(fff);
        kx=kx/max(kx)/dx;
        kz=kz/max(kz)/dz;
        fx=fx/max(fx);
        fz=fz/max(fz);
        subplot(3,1,1);h=loglog(kx,fx,kz,fz);
        legend(h,{'x','z'});
        subplot(3,1,2);semilogy(c.z,ffx(end-5:end,:));xlabel('z');
        subplot(3,1,3);semilogy(c.x,ffz(end-5:end,1:c.Nx));xlabel('x');
      case 3
        x=x(:,:,end/2:end);
        Nx=size(x,3);
        Ny=size(x,2);
        Nz=size(x,1);
        disp(sprintf('%s %u %u %u',g,Nx,Ny,Nz));
        w=ichebintw(Nz);
        
        [fx kx]=pmtm(reshape(permute(x,[3 2 1]),[Nx,Ny*Nz]),[],256,1/dx);
        fx=reshape(fx,[length(kx),Ny,Nz]);
        fx=squeeze(mean(fx,2))*w';
        fy=w*squeeze(mean(abs(fft(x,[],2)).^2,3));
        fz=squeeze(mean(mean(ichebf2c(x,1).^2,3),2));
        ky=imag(fft_modes(Ny,dy)); 
        kz=linspace(0,1,Nz)';
        fff=find(fx>tol & kx>0);fx=fx(fff); kx=kx(fff);
        fff=find(fy>tol & ky>0);fy=fy(fff); ky=ky(fff);
        fff=find(fz>tol & kz>0);fz=fz(fff); kz=kz(fff);
        kx=kx/max(kx)/dx;
        ky=ky/max(ky)/dy;
        kz=kz/max(kz)/dz;
        fx=fx/max(fx);
        fy=fy/max(fy);
        fz=fz/max(fz);
        if sum(sum(sum(abs(diff(x,1,2)))))<1e-6
          h=loglog(kx,fx,kz,fz);
          legend(h,{'x','z'});
        else
          h=loglog(kx,fx,ky,fy,kz,fz);
          legend(h,{'x','y','z'});
        end
    end
    subplot(3,1,1);title(sprintf('%s: %s %s',nm,ff{j},g));
    axis('tight');
    %aa=axis;
    %axis([aa(1)*0.9 aa(2)/0.9 aa(3)*0.9 aa(4)/0.9]);
    drawnow;
    if strcmp(g,'fffu')
      keyboard;
    end
  end
end

return;


L=4;
N=256;
M=32;
w=logspace(-1,1,20);
w=logspace(-0.35,-0.3,20);
w=linspace(1,3,20);
w=linspace(1,5,100);
clf;
for j=1:length(w)
  dx=L/N;
  x=(0:N-1)'*dx;
  y=L/2+(0:M-1)*dx/M*10;
  f=erf((x-y)/(dx*w(j)));
  [fx kx]=pmtm(f,[],64,1/dx);       
  %h(j)=loglog(kx,mean(fx,2));
  ff(j)=mean(fx(end,1));
  hold('on');
end

% $$$ The length scale is sqrt(2*H/(Pe*U)/dx>=2
% $$$ 
% $$$ H=1;U=1;Pe=2000;
% $$$ dx=sqrt(H/(2*Pe*U))
% $$$ dx=0.0158; Should resolve the front
% $$$ 
% $$$ y=linspace(0,dy
% $$$ Nx=linspace(0,L,N);
% $$$ Ny=linspace(
% $$$ x=li
% $$$ 
% $$$     
% $$$ 
% $$$ c=ded_coord('gc/f7/md/00100/27000');
% $$$ a=ded_read_hdf('~/gc/f7/md/00100/27000/state-00031.hdf5');
% $$$ F=ded_read_hdf('~/gc/f7/md/00100/27000/force/fx.hdf5');
% $$$ b=a.b';
% $$$ dx=c.dx;
% $$$ [fx kx]=pmtm(b(end/2:end,:),[],64,1/dx);  
% $$$ 
% $$$ clf;subplot(2,1,1);
% $$$ loglog(kx,mean(fx,2));
% $$$ subplot(2,1,2);
% $$$ f=findmax(fx(end,:));
% $$$ plot(c.x(1:50),b(1:40,f));
% $$$ 
% $$$ 
% $$$ 
% $$$ [fx kx]=pmtm(a.wbl,[],8,1/dx);  
% $$$ loglog(kx,fx);
% $$$ 
% $$$ mpiexec -n 32 ded_gc.py --reset --rfn f7/md/00100/27000 --pfn f7/md/00100/27000 --fbmult True --fbmax False gc/f7/test/01
% $$$ 
% $$$ mkdir -p ~/gc/f7/test
% $$$ rsync -vap hamilton:gc/f7/test/01 ~/gc/f7/test



%mpiexec -n 32 ded_gc.py -rfn md/00100/27000  --pfn md/00100/27000 --reset gc/f7/test/01
%mpiexec -n 32 ded_gc.py --rfn f7/md/00100/27000  --pfn f7/md/00100/27000 --reset gc/f7/test/01
%mpiexec -n 32 ded_gc.py --rfn f7/test/01 --pfn f7/test/01 --db None --reset gc/f7/test/03

ded_aliasing2('gc/f7/test/02');
