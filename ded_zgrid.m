function [b c]=ded_zgrid(a,nz,nm,nmd,nmdd,nmi,H,dim)
%[b c]=ded_zgrid(a,nz,nm,nmd,nmdd,nmi,H) Grid interior chebychev data and calculate
%derivatives if requested
% c is structure of chebychev coefficents
% b is structure of evaluated functions
% H is height of domain for scaling derivative
% nmd is fieldnames or flag for 1st derivatives
% nmdd is fieldnames or flag for 2nd derivatives
if isempty(a)
  b=[];
  c=[];
  return;
end

if nargin<2
  nz=[];
end
if nargin<3
  nm=[];
end
if nargin<4
  nmd={};
end
if nargin<5
  nmdd={};
end
if nargin<6
  nmi={};
end
if nargin<7
  H=[];
end
if nargin<8
  dim=[];
end
if isempty(H)
  H=1;
end

if isempty(nz)
  nz=2*length(a.z);
end
if isempty(dim)
  dim=1;
end

z=linspace(-1,1,nz);
b.z=linspace(0,H,nz);
if isfield(a,'t')
  b.t=a.t;
end
if isfield(a,'y')
  b.y=a.y;
end
if isfield(a,'i')
  b.i=a.i;
end
if isfield(a,'x'); b.x=a.x; end
if isfield(a,'y'); b.y=a.y; end

a=rmfield(a,intersect({'Tz','t','x','y','z','nt','nx','dt','nz','sz','i','constant','iteration','kx','sim_time','timestep','wall_time','world_time','write_number','i1','i2','t1','t2','n','wn','write_num','param','nm','typ','b_parity','bz_parity','p_parity','u_parity','uz_parity','w_parity','wz_parity'},fieldnames(a)));

if isempty(nm)
  nm=setdiff(fieldnames(a),{'z','t','gsize','interval','type'});
end
if ~iscell(nmd)
  if nmd ~=0
    nmd=nm;
  end
end
if ~iscell(nmdd)
  if nmdd ~=0
    nmdd=nm;
  end
end
if ~iscell(nmi)
  if nmi ~=0
    nmi=nm;
  end
end

fna=fieldnames(a);


nm   = intersect(fna,nm);
nmd  = intersect(fna,nmd);
nmdd = intersect(fna,nmdd);

% Function values
for j=1:length(nm)
  e=nm{j};
  c.(e)=ichebf2c(a.(e),dim);
  b.(e)=ichebc2f(c.(e),dim,z);
end

% First derivatives
for j=1:length(nmd)
  e=nmd{j};
  f=['d' e 'dz'];
  if ~isfield(b,f)
    c.(f)=ichebdiffc(c.(e),dim,H);
    b.(f)=ichebc2f(c.(f),dim,z);
  end
end


% Second derivatives
for j=1:length(nmdd)
  e=nmdd{j};
  f=['d' e 'dzz'];
  if ~isfield(b,f)
    c.(f)=ichebdiffc(c.(e),dim,H,2);
    b.(f)=ichebc2f(c.(f),dim,z);
  end
end

% Integrals
for j=1:length(nmi)
  e=nmi{j};
  f=[e 'Iz'];
  if ~isfield(b,f)
    if isfield(c,e)
      c.(f)=ichebintc(c.(e),dim,H);
      b.(f)=ichebc2f(c.(f),dim,z);
    end
  end
end





