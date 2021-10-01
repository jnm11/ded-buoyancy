function [a b]=ded_gc_diff_int(a,p,dx,calc_diff,calc_int,calc_NS)

if nargin<4
  calc_diff=[];
end
if nargin<5
  calc_int=[];
end
if nargin<6
  calc_NS=[];
end

if isempty(calc_diff)
  calc_diff=1;
end

if isempty(calc_int)
  calc_int=1;
end
if isempty(calc_NS)
  calc_NS=1;
end

calc_diff = calc_diff | calc_NS;
calc_int  = calc_int  | calc_NS;

if ~strcmp(p.Tz,'Cheb')
  error('ded_gc_diff_int: expecting Tz is Cheb');
end

if strcmp(p.Tx,'SinCos')
  Pu =-1;
  Pv = 1;
  Pp = 1;
  Pw = 1;
  Pb = 1;
  Pp = 1;
else
  Pu = 0;
  Pv = 1;
  Pp = 0;
  Pw = 0;
  Pb = 0;
  Pp = 0;
end

dimz=1;
dimx=2;

if calc_diff
  if ~isfield(a,'dpdx') &isfield(a,'p');  disp('Calculating dpdx');  a.dpdx  = pr_diff(a.p,dx,dimx,1,[],Pp);     end
  if ~isfield(a,'dbdx') &isfield(a,'b');  disp('Calculating dbdx');  a.dbdx  = pr_diff(a.b,dx,dimx,1,[],Pb);     end
  if ~isfield(a,'dudx') &isfield(a,'u');  disp('Calculating dudx');  a.dudx  = pr_diff(a.u,dx,dimx,1,[],Pu);     end
  if ~isfield(a,'dvdx') &isfield(a,'v');  disp('Calculating dvdx');  a.dvdx  = pr_diff(a.v,dx,dimx,1,[],Pv);     end
  if ~isfield(a,'dwdx') &isfield(a,'w');  disp('Calculating dwdx');  a.dwdx  = pr_diff(a.w,dx,dimx,1,[],Pw);     end
  if ~isfield(a,'dbdxx')&isfield(a,'b');; disp('Calculating dbdxx'); a.dbdxx = pr_diff(a.b, dx,dimx,2,[],Pb);    end
  if ~isfield(a,'dudxx')&isfield(a,'u');; disp('Calculating dudxx'); a.dudxx = pr_diff(a.u, dx,dimx,2,[],Pu);    end
  if ~isfield(a,'dvdxx')&isfield(a,'v');; disp('Calculating dvdxx'); a.dwdxx = pr_diff(a.v, dx,dimx,2,[],Pv);    end
  if ~isfield(a,'dwdxx')&isfield(a,'w');; disp('Calculating dwdxx'); a.dwdxx = pr_diff(a.w, dx,dimx,2,[],Pw);    end
  
  if ~isfield(a,'dbudx')&isfield(a,'bu'); disp('Calculating dbudx'); a.dbudx = pr_diff(a.bu,dx,dimx,1,[],Pb*Pu); end
  if ~isfield(a,'dbvdx')&isfield(a,'bv'); disp('Calculating dbvdx'); a.dbvdx = pr_diff(a.bv,dx,dimx,1,[],Pb*Pv); end
  if ~isfield(a,'dbwdx')&isfield(a,'bw'); disp('Calculating dbwdx'); a.dbwdx = pr_diff(a.bw,dx,dimx,1,[],Pb*Pw); end
  if ~isfield(a,'duudx')&isfield(a,'uu'); disp('Calculating duudx'); a.duudx = pr_diff(a.uu,dx,dimx,1,[],Pu*Pu); end
  if ~isfield(a,'dvvdx')&isfield(a,'vv'); disp('Calculating dvvdx'); a.dvvdx = pr_diff(a.vv,dx,dimx,1,[],Pv*Pv); end
  if ~isfield(a,'dwwdx')&isfield(a,'ww'); disp('Calculating dwwdx'); a.dwwdx = pr_diff(a.ww,dx,dimx,1,[],Pw*Pw); end
  if ~isfield(a,'duvdx')&isfield(a,'uv'); disp('Calculating duvdx'); a.duvdx = pr_diff(a.uv,dx,dimx,1,[],Pu*Pv); end
  if ~isfield(a,'duwdx')&isfield(a,'uw'); disp('Calculating duwdx'); a.duwdx = pr_diff(a.uw,dx,dimx,1,[],Pu*Pw); end
  if ~isfield(a,'dvwdx')&isfield(a,'vw'); disp('Calculating dvwdx'); a.dvwdx = pr_diff(a.vw,dx,dimx,1,[],Pv*Pw); end
end

fn=fieldnames(a);
for j=1:length(fn)
  sz(j,:)=size(a.(fn{j}));
end
f=find(all(sz==repmat(max(sz),size(sz,1),1),2));
nm={fn{f}};

for j=1:length(nm)
  e=nm{j};
  c.(e)=ichebf2c(a.(e),dimz);
end

% First derivatives
if calc_diff
  nmd=intersect(nm,{'b','p','u','v','w','bb','bu','bv','bw','uu','vv','ww','uv','uw','vw'});
  for j=1:length(nmd)
    e=nmd{j};
    f=['d' e 'dz'];
    if ~isfield(c,f)
      c.(f)=ichebdiffc(c.(e),dimz,p.H,1);
    end
  end
  
  nmdd=intersect(nm,{'b','u','v','w'});
  for j=1:length(nmdd)
    e=nmdd{j};
    f=['d' e 'dzz'];
    if ~isfield(c,f)
      c.(f)=ichebdiffc(c.(e),dimz,p.H,2);
    end
  end
end
if calc_diff
  if ~isfield(c,'d');     disp('Calculating d');     c.d     = c.dudx+c.dwdz;    end
  
  if ~isfield(c,'ududx')&isfield(c,'duudx'); disp('Calculating ududx'); c.ududx =  c.duudx/2;end
  if ~isfield(c,'ududz')&isfield(c,'duudz'); disp('Calculating ududz'); c.ududz =  c.duudz/2;end
  if ~isfield(c,'udwdz')&isfield(c,'ududx'); disp('Calculating udwdz'); c.udwdz = -c.ududx;  end
  
  if ~isfield(c,'wdwdx') disp('Calculating wdwdx'); c.wdwdx =  c.dwwdx/2;end
  if ~isfield(c,'wdwdz')&isfield(c,'dwwdx'); disp('Calculating wdwdx'); c.wdwdx =  c.dwwdx/2;end
  if ~isfield(c,'wdwdz')&isfield(c,'dwwdz'); disp('Calculating wdwdz'); c.wdwdz =  c.dwwdz/2;end
  if ~isfield(c,'wdudx')&isfield(c,'wdwdz'); disp('Calculating wdudx'); c.wdudx = -c.wdwdz;  end
  
  if ~isfield(c,'wdudz')&isfield(c,'duwdz');; disp('Calculating wdudz'); c.wdudz  =  c.duwdz - c.udwdz;end
  if ~isfield(c,'udwdx')&isfield(c,'duwdx');; disp('Calculating uduwx'); c.udwdx  =  c.duwdx - c.wdudx;end
  
  if ~isfield(c,'bd')&all(isfield(c,{'bdudx','bdwdz'}));disp('Calculating bd'); c.bd = c.bdudx + c.bdwdz; end
  if ~isfield(c,'ud')&all(isfield(c,{'ududx','udwdz'}));disp('Calculating ud'); c.ud = c.ududx + c.udwdz; end
  if ~isfield(c,'wd')&all(isfield(c,{'wdudx','wdwdz'}));disp('Calculating wd'); c.wd = c.wdudx + c.wdwdz; end
end

if calc_NS
  c.P       = c.p+c.uu/2+c.ww/2;
  c.dPdz    = c.dpdz+c.ududz/2+c.dwwdz/2;
  c.dPdx    = c.dpdx+c.ududx/2+c.dwwdx/2;
end

if calc_int
  nmI=intersect(fieldnames(c),{'b','u','v','w','ud','vw','wd','uu','vv','ww','uv','uw','vw','uw','p','dudx','dvdx','dwdx','dpdx','dbdx','dPdx',...
                      'duudx','dvvdx','dwwdx','duvdx','duwdx','dvwdx','wdudz','udwdz','udwdx','wdudx','P','dudxx','dwdxx','ud','wd'});
  for j=1:length(nmI)
    e=nmI{j};
    f=[e 'Iz'];
    if ~isfield(c,f)
      c.(f)=ichebintc(c.(e),dimz,p.H);
    end
  end
end

cfn=fieldnames(c);
c.z=linspace(0,p.H,2*max(sz(:,1)));
for j=1:length(cfn)
  f=cfn{j};
  c.(f)=ichebc2f(c.(f),dimz,c.z,p.H);
end

a=c;

if calc_int
  a.bIx     = pr_int(a.b,     dx, dimx,[],Pb);
  a.dudzIx  = pr_int(a.dudz,  dx, dimx,[],Pu);
  a.dwdzIx  = pr_int(a.dwdz,  dx, dimx,[],Pw);
  a.dpdzIx  = pr_int(a.dpdz,  dx, dimx,[],Pp);
  a.dPdzIx  = pr_int(a.dPdz,  dx, dimx,[],Pp);
  a.dudzzIx = pr_int(a.dudzz, dx, dimx,[],Pu);
  a.dwdzzIx = pr_int(a.dwdzz, dx, dimx,[],Pw);
  a.wdudzIx = pr_int(a.wdudz, dx, dimx,[],Pu*Pw);
  a.udwdxIx = pr_int(a.udwdx, dx, dimx,[],Pu*Pw);
  
  a.dwwdzIx = pr_int(a.dwwdz, dx, dimx,[],Pw*Pw);
  a.duudzIx = pr_int(a.duudz, dx, dimx,[],Pu*Pu);
  a.duwdzIx = pr_int(a.duwdz, dx, dimx,[],Pu*Pw);
end

if calc_NS
  nu=1/p.Re;
  
  a.wy=a.dudz-a.dwdx;
  
  b.cm    = a.dudx+a.dwdz-a.d;
  b.cpx   = a.duudx  + a.duwdz - a.ud + a.dpdx-nu*(a.dudxx+a.dudzz);
  b.cpz   = a.duwdx  + a.dwwdz - a.wd + a.dpdz-nu*(a.dwdxx+a.dwdzz)+p.g*a.b;
  b.npx   = a.ududx  + a.wdudz        + a.dpdx-nu*(a.dudxx+a.dudzz);
  b.npz   = a.udwdx  + a.wdwdz        + a.dpdz-nu*(a.dwdxx+a.dwdzz)+p.g*a.b;
  
  
  b.cmIx  = a.u+a.dwdzIx;
  b.cmIz  = a.dudxIz+a.w;
  b.cpxIx = a.uu  + a.duwdzIx + a.p      - nu*(a.dudx+   a.dudzzIx);
  b.cpzIx = a.dwwdzIx  + a.uw + a.dpdzIx - nu*(a.dwdx+   a.dwdzzIx)+p.g*a.bIx;
  b.cpxIz = a.duudxIz  + a.uw + a.dpdxIz - nu*(a.dudxxIz+a.dudz);
  b.cpzIz = a.ww  + a.duwdxIz + a.p      - nu*(a.dwdxxIz+a.dwdz)+p.g*a.bIz;
  
  b.npxIx = a.wdudzIx - a.ww/2      + a.P      - nu*(a.dudx    + a.dudzzIx);
  b.npzIx = a.udwdxIx - a.duudzIx/2 + a.dPdzIx - nu*(a.dwdx    + a.dwdzzIx)+p.g*a.bIx;
  b.npxIz = a.wdudzIz - a.dwwdxIz/2 + a.dPdxIz - nu*(a.dudxxIz + a.dudz);
  b.npzIz = a.udwdxIz - a.uu/2      + a.P      - nu*(a.dwdxxIz + a.dwdz)+p.g*a.bIz;
end
