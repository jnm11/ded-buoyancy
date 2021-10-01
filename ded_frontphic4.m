function xm=ded_frontphic4(x,s,display)
%EM_FRONTPHIC4(X,s,DISPLAY)
% Find the front location
% Find the inflection point in the maximum of f
% and then the maximum concentration
% xm is x location

cs(length(s):-1:1) = cummax(s(end:-1:1));
cs=s;
dds=diff(cs,2);
ds =diff(cs,1);
ms = 0.5*(cs(2:end)+cs(1:end-1));
tol=cs(end)+0.1*(cs(1)-cs(end));
f=find(dds(1:end-1)<0 & dds(2:end)>0 & ds(2:end-1)<0 & ms(2:end-1)>tol);
if ~isempty(f)
  mf = max(f);
  xm = (x(mf+1)*dds(mf+1)-x(mf+2)*dds(mf))/(dds(mf+1)-dds(mf));
else
  xm=NaN;
end

return

nm='gc/qle1';
nm='gc/ccle/2d/01'
a=ded_read_g(nm,'yz');
X=0*a.t;
for j=1:a.nt
  X(j)=ded_frontphic4(a.x,a.b(:,j));
end
plot(a.t,X);

mpiexec -n 4 ded_gc.py --preset ccle --Ny 1 --dtyz 0.01 ccle/2d/01

  