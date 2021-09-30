function [psi phi omega div]=ded_helmholtz(a)
a.w(:,[1 end],:)=0;
sz=size(a.w);
sz(end+1:3)=1;
psi=zeros(sz([2 1 3]));
phi=psi;
omega=psi;
div=psi;
dx=a.x(2)-a.x(1);
dz=a.z(2)-a.z(1);
for j=1:size(a.w,3);
  %idiv   =  pr_diff(a.u(:,:,j),dx,2)'+a.wdz(:,:,j)';
  %iomega =  pr_diff(a.w(:,:,j),dx,2)'-a.udz(:,:,j)'; 
  [psi(:,:,j) phi(:,:,j) omega(:,:,j) div(:,:,j)]=helmholtz_fft(a.u(:,:,j)',a.w(:,:,j)',dx,dz);
end
psi   = permute(psi,  [2 1 3]);
phi   = permute(phi,  [2 1 3]);
omega = permute(omega,[2 1 3]);
div   = permute(div,  [2 1 3]);




