function a=ded_pm_r(nm,numi,numo)
% Grid 3d time averaged averaged data radially
%ded_pm_r(nm,num)
%a=ded_pm_r('pm/025',4);
%nm='pm/025';numi=10;numo=10;
dd=ded_dedalus_data_dir;

fnr=sprintf('%s/%s/a/a-r-%05u.mat',dd,nm,numo);
if isfile(fnr) 
  load(fnr)
  if all(a.num==numi)
    return;
  end
  clear('a');
end

p=ded_read_param(nm);

for k=1:length(numi)
  fna=sprintf('%s/%s/a/a-%05u.hdf5',dd,nm,numi(k));
  a=ded_read_hdf(fna);
  if k==1
    dr=sqrt((a.y(2)-a.y(1))*(a.z(2)-a.z(1))/2);
    [yy zz]=ndgrid(a.y,a.z);
    r=sqrt(yy.^2+zz.^2);
    nlz=0;%Normalise
    r=r(:)'/dr;
    nr=length(r);
    rmin=0;
    rmax=ceil(max(r));
    fnm=fieldnames(a);
    sz=[length(a.z) length(a.y) length(a.x)]; 
    for j=1:length(fnm)
      if prod(size(a.(fnm{j})))~=prod(sz)
        if k==1
          b.(fnm{j})=a.(fnm{j});
        end
        continue
      end
      y=reshape(a.(fnm{j}),[nr sz(3)])';
      [pa fa rr]=jgrid1(r,y,rmin,rmax,'cubic');
      if rr(2)~=0
        error('Cannot reflect r');
      end
      fa(:,3)=fa(:,3)+fa(:,1);
      if k==1
        b.(fnm{j})=fa(:,2:end-1);
      else
        b.(fnm{j})=b.(fnm{j})+fa(:,2:end-1);
      end
    end
    b.dt(k) = a.dt;
    b.t1(k) = a.t1;
    b.t2(k) = a.t2;
    b.n(k)  = a.n;
  end
end
b.num = numi;
b.dt = sum(b.dt);
b.n  = sum(b.n);
b.t1 = min(a.t1);
b.t2 = max(a.t2);

pa(3)=pa(3)+pa(1);
for j=1:length(fnm)
  if prod(size(a.(fnm{j})))==prod(sz)
    b.(fnm{j})=b.(fnm{j})./(b.dt*pa(2:end-1));
  end
end


b.rw=pa(2:end-1);
b.r=rr(2:end-1)*dr;
a=b;
save(fnr,'a');


return;
nm='pm/025';
numi=10;
numo=10;

nm='pm/029';
numi=4;
numo=4;

nm='pm/test/01';
numi=20:21;
numo=2021



a=ded_pm_r(nm,numi,numo);


opt=optimset('display','none');

f= @(p,r) p(1)*exp(-p(2)*r.^2)+polyval(p(3:end),r);

param=ded_read_param(nm);

kmax=1+isfield(param,'B')+isfield(param,'S');
pp=zeros([length(a.x) 5,kmax]);
for k=1:kmax
  switch(k)
    case(1)
      p=[2 1 0 0 0];np=5;
      lp=[0 0 -inf -inf -inf];
      x=a.u/param.U;
    case(2)
      p=[0 1 0];np=3;
      lp=[0 0 -inf];
      x=a.b/param.B;
    case(3)
      p=[0 1 0];np=3;
      lp=[0 0 -inf];
      x=a.s/param.S;
   end
  minx=min(x(:));
  maxx=max(x(:));
  for j=1:length(a.x)
    ff= @(p) f(p,a.r)-x(j,:);
    p=lsqnonlin(ff,p,lp,[],opt);
    plot(a.r,x(j,:),a.r,f(p,a.r));
    pp(j,1:np,k)=p;
    axis([0 a.r(end) minx, maxx]);
    title(sprintf('x=%6.2f',a.x(j)));
    drawnow;
  end
end

for j=1:length(a.x)
  plot(a.r,a.d(j,:));
  drawnow;
end

  

RS=zeros(length(a.x),kmax);
RS(:,1)=sqrt(max(0,a.u*a.r.^3')./max(0,a.u*a.r'));
RS(:,2)=sqrt(max(0,a.b*a.r.^3')./max(0,a.b*a.r'));
if isfield(a,'s')
  RS(:,3)=sqrt(max(0,a.s*a.r.^3')./max(0,a.s*a.r'));
end

R=squeeze(1./sqrt(pp(:,1,:)/2));
R(a.x<param.x6,:)=1;
RS(a.x<param.x6,:)=1;
h=plot(a.x,R,a.x,RS);
legend(h,{'u','b','s','u','b','s'});
axis([0 inf 0 10]);



b=a.b;

nx=size(b,1);
R=zeros(nx,1);
tol=0.01;
for j=1:nx
  %  tol=max(b(j,:))/20;
  f1=min(find((b(j,1:end-1)-tol).*(b(j,2:end)-tol)<0));
  if isempty(f1)
    continue;
  end
  b1=b(j,f1);
  b2=b(j,f1+1);
  r1=a.r(f1);
  r2=a.r(f1+1);
  
  % (R-r1)*b2+(r2-R)*b1=tol*(r2-r1)
  % R*(b2-b1)+r2*b1-r1*b2=tol
  R(j) = (r1*b2-r2*b1+tol*(r2-r1))/(b2-b1);
end

