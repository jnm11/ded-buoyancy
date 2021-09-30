function [psi omega div]=ded_psi(u,w,Pu,Pw,H,dx,zz)
%ded_helmholtz_slip(a) Find the streamfunction psi, vorticity omega and divergence.
%in x regions where div is zero psi=int(u,dz)
%in x regions where div~=0 integrate in from edges
% u psi_x + w psi_w = 0  
%   psi_x + w ~= 0    since u=psi_w 
% thus psi = - int w dx


if nargin<7
  zz=[];
end


if ndims(u)>2
  psi=0*u;
  omega=psi;
  div=psi;
  for j=1:size(u,3)
     [psi(:,:,j) omega(:,:,j) div(:,:,j)]=ded_psi(u,w,Pu,Pw,H,dx,zz)
  end
  return;
end



dimx=2;
dimz=1;

% Chebychev transform;
uc = ichebf2c(u,dimz);
wc = ichebf2c(w,dimz);

nz=size(uc,dimz);
nx=size(uc,dimx);

psic   = ichebintc(uc,dimz,H,1,true);
Iwc    = pr_int(wc,dx,dimx,1,Pw);
divc   = pr_diff(uc,dx,dimx,1,[],Pu) + ichebdiffc(wc,dimz,H,[],true);
omegac = pr_diff(wc,dx,dimx,1,[],Pw) - ichebdiffc(uc,dimz,H,[],true);
div    = ichebc2f(divc,dimz,zz,H);
omega  = ichebc2f(omegac,dimz,zz,H);

W=ichebintw(nz)*H/2;

sdiv = W*div.^2/H;
tol=max(max(sdiv)/1e3,1e-10);
sdiv=sdiv>tol;
if false
  f=find(diff([0 sdiv 0])~=0);
  f1=f(1:2:end);
  f2=f(2:2:end)-1;
  for j=1:length(f1)
    rg=f1(j):f2(j);
    if     f1(j) == 1;  psic(:,rg)=psic(:,f2(j)+1)-Iwc(:,rg)+Iwc(:,f2(j)+1); % Integrate from the right
    elseif f2(j) == nx; psic(:,rg)=psic(:,f1(j)-1)+Iwc(:,rg)-Iwc(:,f1(j)-1); % Integrate from the left
    else psic(:,rg)=(psic(:,rg(end)+1)-Iwc(:,rg)+Iwc(:,rg(end)+1)+psic(:,rg(  1)-1)+Iwc(:,rg)-Iwc(:,rg(  1)-1))/2;
    end
  end
  psi=ichebc2f(psic,dimz,zz,H);
else
  
  psi=ichebc2f(psic,dimz,zz,H);
  dz=H*sin((1:nz-1)/nz*pi)*sin(pi/nz/2);
  f=min(find(sdiv));
  sdiv(1:f)=1;
  psi = filter_stream(w,u,psi,repmat(sdiv,nz,1),dz,dx);
% $$$   keyboard
% $$$   x=(0.5:nx-0.5)*dx;
% $$$   z=H*(1-cos((0.5:nz-0.5)'/nz*pi))/2;
% $$$   clf;hold('on');
% $$$   for j=1:2
% $$$     switch(j)
% $$$       case(1);
% $$$         [X Y]=ndgrid(x(end),linspace(z(1),z(end),20));
% $$$         XY = stream2(x,z,u,w,X,Y);
% $$$       case (2);
% $$$         [X Y]=ndgrid(x(1),linspace(z(1),z(end),20));
% $$$         XY = stream2(x,z,-u,-w,X,Y);
% $$$     end
% $$$     for j=1:length(XY)
% $$$       if ~isempty(XY{j})
% $$$         plot(XY{j}(:,1),XY{j}(:,2));
% $$$       end
% $$$     end
% $$$   end
  %psi = filter_interpnan2(psi,repmat(sdiv,nz,1));
end


return

nm='gc/f7/ma/00145/25125';
c=ded_coord(nm);
dx=c.dJx;
p=ded_read_param(nm);
H=p.H;
fns=ded_get_fn(nm,'final',[],'state');
a=ded_read_hdf(fns{2});
x=c.Jx';
z=c.Jz;

dc;[psi omega div]=ded_psi(a.u,a.w,a.u_parity(1),a.w_parity(1),H,dx);contour(x,z,psi,50);
