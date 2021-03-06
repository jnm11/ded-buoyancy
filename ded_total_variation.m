function ded_total_variation(nm)
%ded_aliasing('gc/f6/f/20');
%ded_aliasing('gc/test/20');
%nm='gc/f6/f/20';

c=ded_coord(nm);
p=ded_read_param(nm);
ffn=0;
ff={'fw','fu','fb','wb','wu','final','a'};
ff={'final'};
tol=1e-25;
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
    x=a.(g);
    ffn=ffn+1;
    figure(ffn);
    clf;
    switch(c.dim)
      case 2
        Nx=size(x,2);
        Nz=size(x,1);
        dx=p.L/Nx;
        dz=p.H/Nz;
        w=ichebintw(Nz);
        
        [fx kx]=pmtm(permute(x,[2 1]),[],256,1/dx);
        fx=fx*w';
        %k=kx/(2*pi);
        %fx=w*abs(fft(x,[],2)).^2;
        fz=mean(ichebf2c(x,1).^2,2);
        %kx=imag(fft_modes(Nx,dx)); fx=fx(kx>0); kx=kx(kx>0);
        kz=imag(fft_modes(Nz,dz)); fz=fz(kz>0); kz=kz(kz>0)/(2*pi);
        h=loglog(kx,fx,kz,fz);
        legend(h,{'x','z'});
      case 3
        Nx=size(x,3);
        Ny=size(x,2);
        Nz=size(x,1);
        dx=p.L/Nx;
        dy=p.W/Ny;
        dz=p.H/Nz;
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
    title(sprintf('%s: %s %s',nm,ff{j},g));
    axis('tight');
    aa=axis;
    axis([aa(1)*0.9 aa(2)/0.9 aa(3)*0.9 aa(4)/0.9]);
    drawnow;
    if strcmp(g,'fffu')
      keyboard;
    end
  end
end



    