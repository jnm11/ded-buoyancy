function b=ded_diff(a,p,parity,c,nm)

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

T=flip(p.T);
X=flip(p.X);
dd=flip(c.dd);
cc=flip(c.c);
for j=1:length(nm)
  e=nm{j};
  for k=1:length(ccc)
    f=['d' e 'd' ccc(k)];
    pk=find(cc==ccc(k));
    switch(T{pk})
      case 'Cheb'
        b.(f)=ichebdifff(a.(e),k,X(pk)); % Assume spacing is Chebychev points
      otherwise
        b.(f)=pr_diff(a.(e),dd(pk),k,1,[],parity.(e)(c.dim+1-pk));
    end
  end
end

if nargin<5
if c.dim==3
  if isfield(a,'ududx') & isfield(a,'vdudy') & isfield(a,'wdudz')
    b.cdu   = b.duudx +b.duvdy +b.duwdz;
    b.cdv   = b.duvdx +b.dvvdy +b.dvwdz;
    b.cdw   = b.duwdx +b.dvwdy +b.dwwdz;
  end
else
  if isfield(b,'duudx') & isfield(b,'duwdz') & isfield(b,'dwwdz') & isfield(b,'udwdz') & isfield(b,'wdudz')
    b.cdu   = b.duudx  +b.duwdz;
    b.cdw   = b.duwdx  +b.dwwdz;
    b.wdudx = a.wd     -b.dwwdz/2; %b.wdudx=b.w(dudx+dwdz) - b.wdwdz  
    b.udwdz = a.ud     -b.duudx/2; %b.udwdz=b.u(dwdz+dudx) - b.ududx
    b.wdudz = b.duwdz - b.udwdz;
    b.udwdx = b.duwdx - b.wdudx;
  end
end
nm1={'d','s','b','p','u','v','w'};
for j=1:length(nm1)
  for k=1:c.dim
    w1=[nm1{j} 'd' nm1{j} 'd' c.c(k)];
    w2=['d' nm1{j} nm1{j} 'd' c.c(k)];
    if isfield(b,w2)
      b.(w1)=b.(w2)/2;
    end
  end
end
end
