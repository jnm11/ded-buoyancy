function b=ded_gc_find_front_mom(nm,T)

if nargin<2
  T=[];
end

if isstruct(nm)
  a=nm;
  b.T=T;
else
  a=ded_read_mom(nm);
  b.T=ded_convergence_T(nm);
end
if isempty(a)
  b.X=[];
  b.t=[];
  b.Xm=NaN;
else
  if isfinite(b.T)
    T=b.T;
  else
    T=a.t(end);
  end
  b.X = a.mx+sqrt(3)*a.sx;
  f=find(a.t>=min(T,a.t));
  mb=sum(a.b(f));
  mx=sum(a.bx(f))/mb;
  sx=sqrt(sum(a.bxx(f))/mb-mx^2);
  b.Xm = mx+sqrt(3)*sx;
end
