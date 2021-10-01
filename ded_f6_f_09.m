%rsync -uavp hamilton:gc/f6/f/09 ~/gc/f6/f --exclude "check*" --exclude "final*" 
%rsync -uavp hamilton:gc/f6/f/09/b/b-00075.hdf5 ~/gc/f6/f/b

clear('all');

famat='~/gc-f6-f-09a.mat';
fnmat='~/gc-f6-f-09.mat';
nm='gc/f6/f/09';

if ~isfile(famat)
  [a p]=ded_gc_read_javrg(nm,[0 inf]);
  save(famat,'a','p');
else
  load(famat);
end
if ~isfile(fnmat)
  [a p]=ded_gca_3(a,p,nm);
  save(fnmat,'a','p');
else
  load(fnmat);
end

sz=[6 3];
fnt=12;
nx=length(a.x);
nz=length(a.z);

x0=39.8;

u=-a.u(1,:);
f=1+max(find(u<0)):length(u);
plot(a.x(f),u(f))

uu=u(f);
xx=a.x(f);
xx=xx-xx(1);
xx=xx(:);
uu=uu(:);

ff = @(p) (p(1)*sqrt(xx)+p(2)*xx)./(1 + p(2)*xx);

fff = @(p) ff(p)-uu;
pp=lsqnonlin(fff,[0.5124    5.0810],[0 0]);
plot(xx,uu,xx,ff(pp));
axis([0 2 0 1]);






b=ded_read_hdf('~/gc/f6/f/09/b/b-00012.hdf5');
b.b=squeeze(b.b(1,:,:));




f=find(a.x>5);
x=a.x'-x0;
h=max(0,a.bIz(end,:));
f1=findmin(abs(a.x-5));
f2=findmax(h.*(a.x>30));
Fr=p.U./sqrt(p.g*h);
plot(a.x(f1:f2)-x0,Fr(f1:f2));




