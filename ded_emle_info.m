function a=ded_emle_info(nm,a)
if nargin<2
  a=[];
end
if iscell(nm)
  for j=1:length(nm)
    a(j)=ded_emle_info(nm{j});
  end
  return;
end
a.dimple=0;
a.noise=0.005;
p=ded_read_param(nm);
switch(nm)
  case 'gc/ccle/046'
    a.noise=0.0025;
  case 'gc/emle/001'
    a.noise=0.0025;
  case 'gc/emle/003'
    a.noise=0.01;
  case 'gc/emle/011'
    a.dimple=0.005;
  case 'gc/emle/012'
    a.dimple=0.010;
  case 'gc/emle/013'
    a.dimple=0.020;
  case 'gc/emle/014'
    a.dimple=0.040;
end
a.Wx=p.Wx;
a.dx=p.L/p.Nx;
a.Wn=a.Wx/a.dx;
a.Sc=p.Scb;
a.Nz=p.Nz;
a.W=p.W;
a.L=p.L;
