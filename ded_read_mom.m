function a=ded_read_mom(n)

a1=ded_read_g(n,'mom');
a2=ded_read_g(n,'momb');
if isempty(a1) & isempty(a2)
  a=[];
  return;
elseif ~isempty(a1) & ~isempty(a2)
  a1=lower_fieldnames(a1);
  a.t   = [      a1.t;     a2.t];
  a.b   = [      a1.b;     a2.b];
  a.bx  = [     a1.bx;    a2.bx];
  a.bxx = [    a1.bxx;   a2.bxx];
  a.bxz = [    a1.bxz;   a2.bxz];
  a.bz  = [     a1.bz;    a2.bz];
  a.bzz = [    a1.bzz;   a2.bzz];
  
  a.bxy = [    NaN*a1.t;   a2.bxy];
  a.byy = [    NaN*a1.t;   a2.byy];
  a.byz = [    NaN*a1.t;   a2.byz];
  a.nt  =length(a.t);
elseif ~isempty(a1) 
  a=lower_fieldnames(a1);
else
  a=a2;
end

if isempty(a)
  return;
end
a.mx  = a.bx ./a.b;
a.mxx = a.bxx./a.b;
a.mz  = a.bz ./a.b;
a.mzz = a.bzz./a.b;
a.mxz = a.bxz./a.b;
a.sx  = sqrt(max(0,a.mxx-a.mx.^2));
a.sz  = sqrt(max(0,a.mzz-a.mz.^2));
a.sxz = sqrt(max(0,a.mxz-a.mx.*a.mz));
f=findstr('dedalus',n);
if isempty(f)
  a.nm=n;
else
  a.nm=n(f+8:end);
end



