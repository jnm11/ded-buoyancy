function a=ded_add_int(a,p,parity,c,nm)

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
    f=['I' e 'd' ccc(k)];
    pk=find(cc==ccc(k));
    if isfield(a,f)
      continue;
    end
    switch(T{pk})
      case 'Cheb'
        a.(f)=ichebintf(a.(e),k,[],X(pk),[],true); 
      otherwise
        a.(f)=pr_int(a.(e),dd(pk),k,1,-parity.(e)(c.dim+1-pk));
    end
  end
end
