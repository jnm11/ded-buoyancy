function [psi phi w d sr rr q]=ded_helmholtz3(b,p)
Nx=size(b.u,2);
Nz=size(b.u,1);
u=b.u;
w=b.w;

if strcmp(p.Tx,'SinCos') & 0
  u=[b.u -fliplr(b.u)];
  w=[b.w  fliplr(b.w)];
end
u=[u; flipud(u)];
w=[w;-flipud(w)];

[psi phi w d sr rr q]=helmholtz_fft(w,u,b.dz,b.dx);
psi=psi(1:Nz,1:Nx);
phi=phi(1:Nz,1:Nx);
