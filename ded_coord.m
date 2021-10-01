function c=ded_coord(nm)
fnc=[ded_dedalus_data_dir '/' nm '/' 'coord.hdf5'];
fnp=[ded_dedalus_data_dir '/results/' nm '/' 'coord.mat'];
if ~isfile(fnc) & ~isfile(fnp)
  c=[];
  return;
end

if isfile(fnc)
  c=ded_read_hdf(fnc);
  if ~isfield(c,'Ax');
    c=[];
    return;
  end
elseif isfile(fnp)
  c=load(fnp);
  c=c.c;
  return;
end
  
p=ded_read_param(nm);

if isfield(p,'Tx') & ~isfield(c,'x');c=calc(c,'x',p.L,p.Nx,p.Tx,p);end
if isfield(p,'Ty') & ~isfield(c,'y');c=calc(c,'y',p.W,p.Ny,p.Ty,p);end
if isfield(p,'Tz') & ~isfield(c,'z');c=calc(c,'z',p.H,p.Nz,p.Tz,p);end

c.dx  = c.x(2) -c.x(1);
c.dAx = c.Ax(2)-c.Ax(1);
c.dSx = c.Sx(2)-c.Sx(1);
c.dJx = c.Jx(2)-c.Jx(1);
c.Nx  = length(c.x);
c.NAx = length(c.Ax);
c.NSx = length(c.Sx);
c.NJx = length(c.Jx);
c.dz  = c.z(2) -c.z(1);
c.dAz = c.Az(2)-c.Az(1);
c.dSz = c.Sz(2)-c.Sz(1);
c.dJz = c.Jz(2)-c.Jz(1);
c.Nz  = length(c.z);
c.NAz = length(c.Az);
c.NSz = length(c.Sz);

c.NJz = length(c.Jz);
if isfield(c,'y')
  if length(c.y)<=2
    c=rmfieldifexist(c,{'Ny','y','Ty','Jy','Sy','NAy','NJy','NSy'});
  end
end
if isfield(c,'y')
  c.dy  = c.y(2) -c.y(1);
  c.dAy = c.Ay(2)-c.Ay(1);
  c.dSy = c.Sy(2)-c.Sy(1);
  c.dJy = c.Jy(2)-c.Jy(1);
  c.Ny  = length(c.y);
  c.NAy = length(c.Ay);
  c.NSy = length(c.Sy);
  c.NJy = length(c.Jy);
  c.szx  = [1 1 c.Nx]; 
  c.szAx = [1 1 c.NAx]; 
  c.szSx = [1 1 c.NSx]; 
  c.szJx = [1 1 c.NJx]; 
  c.szy  = [1 c.Ny  1]; 
  c.szAy = [1 c.NAy 1]; 
  c.szSy = [1 c.NSy 1]; 
  c.szJy = [1 c.NJy 1]; 
  c.szz  = [c.Nz  1 1]; 
  c.szAz = [c.NAz 1 1]; 
  c.szSz = [c.NSz 1 1]; 
  c.szJz = [c.NJz 1 1]; 
  c.dim  = 3;
  c.dimx = 3;
  c.dimy = 2;
  c.dimz = 1;
  c.dd  = [ c.dx  c.dy  c.dz];
  c.ddA = [c.dAx c.dAy c.dAz];
  c.ddJ = [c.dJx c.dJy c.dJz];
  c.c   = 'xyz';
else
  c.szx  = [1 c.Nx]; 
  c.szAx = [1 c.NAx]; 
  c.szSx = [1 c.NSx]; 
  c.szJx = [1 c.NJx]; 
  c.szz  = [c.Nz  1]; 
  c.szAz = [c.NAz 1]; 
  c.szSz = [c.NSz 1]; 
  c.szJz = [c.NJz 1]; 
  c.dim  = 2;
  c.dimx = 2;
  c.dimz = 1;
  c.dd  = [ c.dx  c.dz];
  c.ddA = [c.dAx c.dAz];
  c.ddJ = [c.dJx c.dJz];
  c.c   = 'xz';
end


function c=calc(c,a,L,N,T,p)
c.(['N' a])  = N;
NA = round(N*p.AA);
NJ = round(N*p.AAJ);
NS = round(N*p.AAS);
L=[0 L];    
switch(T)
  case 'SinCos'
    x  = L(1)+(0.5:N -0.5)/N *(L(2)-L(1));
    Ax = L(1)+(0.5:NA-0.5)/NA*(L(2)-L(1));
    Jx = L(1)+(0.5:NJ-0.5)/NJ*(L(2)-L(1));
    Sx = L(1)+(0.5:NS-0.5)/NS*(L(2)-L(1));
  case 'Fourier'
    if strcmp(p.sType,'pm');L=[-L(2)/2 L(2)/2];end;
    x  = L(1)+(0:N -1)/N *(L(2)-L(1));
    Ax = L(1)+(0:NA-1)/NA*(L(2)-L(1));
    Jx = L(1)+(0:NJ-1)/NJ*(L(2)-L(1));
    Sx = L(1)+(0:NS-1)/NS*(L(2)-L(1));
  case 'Cheb'
    x  = (L(1)+L(2))/2+acos(((N -0.5):-2:(0.5-N ))/N )*(L(2)-L(1))/2;
    Ax = (L(1)+L(2))/2+acos(((NA-0.5):-2:(0.5-NA))/NA)*(L(2)-L(1))/2;
    Jx = (L(1)+L(2))/2+acos(((NJ-0.5):-2:(0.5-NJ))/NJ)*(L(2)-L(1))/2;
    Sx = (L(1)+L(2))/2+acos(((NS-0.5):-2:(0.5-NS))/NS)*(L(2)-L(1))/2;
end

c.(['' a])  = x;
c.(['A' a]) = Ax;
c.(['J' a]) = Jx;
c.(['S' a]) = Sx;
c.(['NA' a]) = NA;
c.(['NJ' a]) = NJ;
c.(['NS' a]) = NS;

