%function ded_gc_left(nm)
%function ded_gc_left(nm) Look at left conditions
nm='gc/test/f7/05';
p=ded_read_param(nm);
c=ded_coord(nm);
fns=ded_get_fn(nm,'final',[],'state');
a=ded_read_hdf(fns{end});
w=ichebintw(c.NSz);

fx=ded_read_hdf([ded_dedalus_data_dir '/' nm '/force/fx.hdf5']);
fd=ded_read_hdf([ded_dedalus_data_dir '/' nm '/force/fd.hdf5']);

U=p.V*pr_int(fd.fd,c.dJx,c.dimx,1,1);

dudz=ichebdifff(a.u,c.dimz,p.L);

fu=dudz./(w*dudz);

x=c.Sx;

qu=w*max(0,a.u); qU=w*max(0,U); 
hu=w*(a.u>0);    hU=w*(U>0);    
uu=qu./hu;       uU=qU./hU;         

%(1-s)*a.u(1:end-1,:)+s*a.u(2:end,:)=0
s=-a.u(1:end-1,:)./diff(a.u);
s(s<0|s>1)=NaN;
s(:,sum(isfinite(s))~=1)=NaN;
[ss sf]=min(s);
hu=(1-ss).*c.Jz(sf)'+ss.*c.Jz(sf+1)';

for j=1:size(a.u,2);
  pp(j,:)=fit_erf(c.Sz,a.u(:,j));
end
eU2=pp(:,1)+pp(:,2);
eU1=pp(:,1)-pp(:,2);
ehu=pp(:,4);
ewu=1./pp(:,3);

subplot(3,1,1);
plot(x,qu,x,qU);
ylabel('qu');
axis('tight');

subplot(3,1,2);
plot(x,hu,x,ehu,x,hU);
ylabel('hu');
axis([0 p.L 0 p.H/2]);

subplot(3,1,3);
plot(x,uu,x,uU);
ylabel('qu./hu');
xlabel('h');
axis('tight');

return;

%hu 0.65
%U1 1.2

%ded_gc.py --rfn test/f7/02  --pfn test/f7/02 --hu 0.65 --U1 1.2 --reset gc/test/f7/05
%--Wz 0.7