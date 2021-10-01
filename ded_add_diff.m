function a=ded_add_diff(a,p,parity,c,nm)

if nargin<5
  nm=fieldnames(a);
end
if isfield(a,'c')
  ccc=a.c;
else
  if c.dim==3
    ccc='zyx';
  else
    ccc='zx';
  end
end
% If have derivatives directly from the simulation rename them
if isfield(a,'bz'); a.dbdz=a.bz; a=rmfield(a,'bz');end; 
if isfield(a,'sz'); a.dsdz=a.sz; a=rmfield(a,'sz');end; 
if isfield(a,'uz'); a.dudz=a.uz; a=rmfield(a,'uz');end; 
if isfield(a,'vz'); a.dvdz=a.vz; a=rmfield(a,'vz');end; 
if isfield(a,'wz'); a.dwdz=a.wz; a=rmfield(a,'wz');end; 

T=flip(p.T);
X=flip(p.X);
dd=flip(c.dd);
cc=flip(c.c);
for j=1:length(nm)
  e=nm{j};
  for k=1:length(ccc)
    f=['d' e 'd' ccc(k)];
    pk=find(cc==ccc(k));
    if isfield(a,f)
      continue;
    end
    switch(T{pk})
      case 'Cheb'
        a.(f)=ichebdifff(a.(e),k,X(pk)); % Assume spacing is Chebychev points
      otherwise
        a.(f)=pr_diff(a.(e),dd(pk),k,1,[],parity.(e)(c.dim+1-pk));
    end
  end
end
