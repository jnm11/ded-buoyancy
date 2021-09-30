function psi=ded_streamfunction(L,H,u,v)
%psi=ded_streamfunction(L,H,u,v) find streamfunction for an incompressible velocity field
%System is periodic in first dimension 
% and interior chebyshev in second direction


m=size(u,1);
n=size(u,2);
dx=L/m;
wx = fft_modes(m,dx);

uc=ichebf2c(u,1,H);

w=pr_diff(v,dx,2)-ichebdiffc(uc,1,H); % Calculate vorticity

fw=fft2(w);

% Knockout the top mode which is ambigous in directions with even size
if mod(m,2)==0
  fw(1+m/2,:)=0;
end
% In fourier space we have f_yy+wx^2f=0 
f = ifft2(-fw./wx.^2, 'symmetric');