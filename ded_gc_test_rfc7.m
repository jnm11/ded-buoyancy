function ded_gc_test_rfc7(nm)
%nm='gc/test/001';
param=ded_read_param(nm);
%a=ded_read_hdf('force/force_s1.hdf5');
%ded_gc_test_rfc7('gc/f7/g/1400/1000');
%rsync -vap asahi:gc/f7/g/1400/1000 ~/gc/f7/g/1400

a1=ded_read_hdf([ded_dedalus_data_dir '/' nm '/bmax.hdf5']);	
a2=ded_read_hdf([ded_dedalus_data_dir '/' nm '/bmin.hdf5']);
u1=ded_read_hdf([ded_dedalus_data_dir '/' nm '/u1.hdf5']);	
u2=ded_read_hdf([ded_dedalus_data_dir '/' nm '/u2.hdf5']);	
fns=sort(cellstr_ls([ded_dedalus_data_dir '/' nm '/a/*']));
a=ded_read_hdf(fns{end});
a3=ded_read_hdf([ded_dedalus_data_dir '/' nm '/rfu.hdf5']);	

w=ichebintw(param.Nz);
ww=ichebintw(param.Nz*param.AA);
if isempty(a3)
  rfu=NaN*a1.bmaxf;
else
  rfu=a3.rfu;
end
x=a.x';
z=a.z;
xx=a1.x';
zz=a1.z;
dxx=xx(2)-xx(1);
dx=x(2)-x(1);

dimx=2;
dimz=1;
Pu=-1;
Pp=1;

N=3;
exx = [xx xx+1*param.L xx+2*param.L]/param.L;
ex  = [ x  x+1*param.L  x+2*param.L]/param.L;

u=a.u/a.dt;
B=a.b/a.dt;
p=a.p/a.dt;
uu=a.uu/a.dt;
dudx=pr_diff(u,dx,dimx,1,[],Pu);
duudx=pr_diff(uu,dx,dimx,1,[],Pu^2);
dpdx=pr_diff(p,dx,dimx,1,[],Pp);
rfuIx=pr_int(rfu,dx,dimx,1,Pu);
BIz=ichebintf(B,dimz,z,param.H);
M=w*(uu+p);
P=uu/2+p; % The same


if param.conservative
  F=1
else
  F=0.5;
end


for j=1:1
  if j==2;
    w=0*w;w(end)=1;
    ww=0*ww;ww(end)=1;
  end
  
  
  emaxb = ded_parity(max(a1.bmaxf,[],1),   Pp,N,dimx);
  eminb = ded_parity(a2.bminf',   Pp,N,dimx);
  eu1   = ded_parity(u1.u1',      Pu,N,dimx);
  eu2   = ded_parity(u2.u2',      Pu,N,dimx);
  erfu  = ded_parity(ww*rfu,      Pu,N,dimx);
  ecux  = ded_parity(F*2*w*(u.*dudx),Pu,N,dimx);
  ecuu  = ded_parity(F*w*duudx,    Pu,N,dimx);
  ep    = ded_parity(w*p,        Pp,N,dimx);
  eM    = ded_parity(M,          Pp,N,dimx);
  edpdx = ded_parity(w*dpdx,     Pu,N,dimx);

  Irfu  = ded_parity(ww*rfuIx, 1,N,dimx);
  Icux  = ded_parity(F*w*u.^2,   1,N,dimx);
  Icuu  = ded_parity(F*w*uu,     1,N,dimx);
  
  Irfu=Irfu-Irfu(end);
  Icuu=Icuu-Icuu(end);
  Icux=Icux-Icux(end);
  ep=ep-ep(end);
  eM=eM-eM(end);
  
  figure(1)
  subplot(2,1,2*j-1);
  h=plot(exx,erfu,ex,ecux,ex,ecuu,ex,edpdx,ex,ecuu+edpdx);
  legend(h,{'rfu','2ududx','duudx','dpdx','dMdx'},'location','eastoutside');
  subplot(2,1,2*j);
  h=plot(exx,Irfu,ex,Icux,ex,Icuu,ex,ep,ex,eM);
  legend(h,{'rfu','u^2','uu','p','M'},'location','eastoutside');

% $$$   figure(2);
% $$$   subplot(3,1,1);plot(exx,erfu);ylabel('rfu');axis([0 N -4.05 0.05]);
% $$$   subplot(3,1,2);plot(exx,eu1);ylabel('u1');axis([0 N -1.05 0.05]);
% $$$   subplot(3,1,3);plot(exx,eu2);ylabel('u2');axis([0 N -1.05 0.05]);


% $$$   figure(3)
% $$$   subplot(2,1,1);plot(exx,emaxb);ylabel('bmaxf');axis([0 N -0.05 1.05]);
% $$$   subplot(2,1,2);plot(exx,eminb);ylabel('bminf');axis([0 N -0.05 1.05]);
end

%a2=ded_read_hdf([ded_dedalus_data_dir '/' nm '/bmin.hdf5']);	

% $$$ figure(4)
% $$$ mesh(exx,zz,ded_parity(-rfu,Pp,N,dimx));

figure(2);
subplot(2,1,1);
plot(z,u(1:20));
