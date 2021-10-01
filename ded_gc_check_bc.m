function ded_gc_check_bc(nm)
% Check the boundary conditions
% Check it for PIDX>0 and non PIDX

%ded_gc-1.31.py --qgc --lbc -U 3 gc/bc1
if nargin<1
  dd=ded_dedalus_data_dir;
  nms=cellstr_ls([dd '/gc/test/*']);
  nms=cellstrremoveprefix(nms,[dd '/']);
  for j=1:length(nms)
    ded_gc_check_bc(nms{j});
  end
  return;
end

a=ded_read_g(nm,'y');
bl=ded_read_g(nm,'left');
br=ded_read_g(nm,'right');
c=ded_read_g(nm,'xyz');
p=ded_read_param(nm);
s=ded_read_stats(nm);

if isempty(bl)
  disp(sprintf('ded_gc_check_bc: No left boundary information %s', nm));
  return;
end
if isempty(br)
  disp(sprintf('ded_gc_check_bc: No right boundary information %s', nm));
  return;
end
if isempty(a)
  disp(sprintf('ded_gc_check_bc: No y-integrated fields %s', nm));
  return;
end

nd=2+(p.Ny~=1);
if nd==3
  bl.v = +bl.v/p.W;
  br.v = +br.v/p.W;
  a.v  = -a.v/p.W;
  a=ded_zgrid(a,a.nz,{'u','v'},{'u'},{},{},p.H);
else
  a=ded_zgrid(a,a.nz,{'u'},{'u'},{},{},p.H);
end

bl.u = -bl.u/p.W;
br.u = -br.u/p.W;
a.u  = -a.u/p.W;

x=a.x;
z=a.z;

t0=1;
fa=find(a.t  >= t0);
fl=find(bl.t >= t0);
fr=find(br.t >= t0);

a.t  = a.t(fa);
a.lu   = squeeze(a.u(1,:,fa));
a.ru   = squeeze(a.u(end,:,fa));
a.ludz = squeeze(a.udz(1,:,fa));
a.rudz = squeeze(a.udz(end,:,fa));

if nd==3
  a.lv   = squeeze(a.v(1,:,fa));
  a.rv   = squeeze(a.v(end,:,fa));
  bl.v   = bl.v(:,fl);
  br.v   = br.v(:,fl);
  lbv = bl.v;
  lav = a.lv;
  rbv = br.v;
  rav = a.rv;
end

bl.t   = bl.t(fl);
bl.u   = bl.u(:,fl);
bl.udz = bl.udz(:,fl);

br.t   = br.t(fl);
br.u   = br.u(:,fl);
br.udz = br.udz(:,fl);

[ts sl] = ded_interp_stats(s,p,bl.t);
[ts sr] = ded_interp_stats(s,p,br.t);
[ts sa] = ded_interp_stats(s,p,a.t);

% $$$ if p.PIDX>0
% $$$   U=-c.u/(p.U*p.L*p.W*p.H);
% $$$   U=U(:)';
% $$$   Ua=interp1(c.t,U,a.t,'linear')';
% $$$   Ul=interp1(c.t,U,bl.t,'linear')';
% $$$   Ur=interp1(c.t,U,br.t,'linear')';
% $$$ else
% $$$   U=1;
% $$$   Ua=1;
% $$$   Ul=1;
% $$$   Ur=1;
% $$$ end


switch(p.lbc)
  case 'noslip'
    lbu = bl.u-repmat(sl.U',size(bl.u,1),1);
    lau = a.lu-repmat(sa.U',size(a.lu,1),1);
  case 'slip'
    lbu = bl.udz;
    lau = a.ludz;
  otherwise
    error('unknown left boundary condition');
end

switch(p.ubc)
  case 'noslip'
    rbu = br.u-repmat(sr.U',size(br.u,1),1);
    rau = a.ru-repmat(sa.U',size(a.ru,1),1);
  case 'slip'
    rbu = br.udz;
    rau = a.rudz;
  otherwise
    error('unknown right boundary condition');
end

clf;
ah=jsubplot([2,4],[0.05 0.05],[0.05 0.05],[0.01 0.05]);

if nd==3
  dedgcbc3(a.t,bl.t,x,p.lbc,ah(1,:),lbu,lbv,lau,lav);
  dedgcbc3(a.t,br.t,x,p.ubc,ah(2,:),rbu,rbv,rau,rav);
else
  dedgcbc2(a.t,bl.t,x,p.lbc,ah(1,:),lbu,lau);
  dedgcbc2(a.t,br.t,x,p.ubc,ah(2,:),rbu,rau);
end
%set(ah,'yscale','log');


function  dedgcbc3(at,bt,x,titl,aa,bu,bv,au,av)
axes(aa(1));
h=plot(x,sqrt(mean(bu.^2,2)),'r',x,sqrt(mean(bv.^2,2)),'b');
legend(h,{'u','v'},'location','best');
title(titl);
xlabel('x');
axis('tight');
axes(aa(2));
plot(bt,sqrt(mean(bu.^2,1)),bt,sqrt(mean(bv.^2,1)))
xlabel('t');
axis('tight');
axes(aa(3));
h=plot(x,sqrt(mean(au.^2,2)),'r',x,sqrt(mean(av.^2,2)),'b');
legend(h,{'u','v'},'location','best');
xlabel('x');
axis('tight');
axes(aa(4));
plot(at,sqrt(mean(au.^2,1)),at,sqrt(mean(av.^2,1)))
xlabel('t');
axis('tight');
return

function  dedgcbc2(at,bt,x,titl,aa,bu,au)
axes(aa(1));
h=plot(x,sqrt(mean(bu.^2,2)),'r');
legend(h,{'u'},'location','best');
title(titl);
xlabel('x');
axis('tight');
axes(aa(2));
plot(bt,sqrt(mean(bu.^2,1)),'r')
xlabel('t');
axis('tight');
axes(aa(3));
h=plot(x,sqrt(mean(au.^2,2)),'r');
legend(h,{'u'},'location','best');
xlabel('x');
axis('tight');
axes(aa(4));
plot(at,sqrt(mean(au.^2,1)))
xlabel('t');
axis('tight');
return

