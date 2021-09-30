function a=ded_s_stats(nm,typ)
d=[ded_dedalus_data_dir '/' nm '/' typ '-stats.mat'];
if isfile(d) 
  load(d);
  return;
end
NN=1000;
MM=100;
fn=ded_get_fn(nm,typ);
n=length(fn);
a.n=zeros(NN,1);
a.x=zeros(MM,n);
a.max=zeros(1,n);
a.min=zeros(1,n);
a.t=zeros(1,n);
a.nm=nm;
a.typ=typ;
h=linspace(0,10,NN);
for j=1:n
  disp(fn{j});
  x=ded_read_hdf(fn{j});
  s=sort(x.(typ)(:));
  n=histc(s,h);
  a.n=a.n+n;
  a.max(j)=s(end);
  a.min(j)=s(1);
  a.x(end-MM+1:end,j)=s(end-MM+1:end);
end
a.h=h;
save(d,'a');

return;
a=ded_s_stats('pm/f7/e/25','s');
a=ded_s_stats('pm/f7/e/24','s');

a=ded_s_stats('pm/f7/e/26','s');
a=ded_s_stats('pm/f7/e/27','s');
a=ded_s_stats('pm/f7/e/28','s');

rg=2:max(find(a.n>0));

semilogy(a.h,a.n);
a=ded_s_stats(nm,'s');



p=ded_read_param(nm);

b=ded_s_stats(nm,'b');

uz=squeeze(sum(sum(abs(fft(b.u,[],1)).^2,2),3));
uy=squeeze(sum(sum(abs(fft(b.u,[],2)).^2,1),3));
ux=squeeze(sum(sum(abs(fft(b.u,[],3)).^2,1),2));
nx=length(ux);fx=1:nx/2;ux=ux(fx)
ny=length(uy);fy=1:ny/2;uy=uy(fy)
nz=length(uz);fz=1:nz/2;uz=uz(fz)


p=polyfit(log(fx'),log(ux),1);

clf;loglog(fx,ux,fy,uy,fz,uz);
hold('on');
plot(fx,1e7*fx.^(-5/3));
axis([0    2.4594   -0.2851    8.4614]);


c=ded_coord(nm);
[Pxx,F] = pwelch(reshape(b.u,[c.Nz,c.Ny*c.Nx]),[],[],[],1/c.dz);
loglog(F,mean(Pxx,2),F,1e-4*F.^(-5/3));
us=squeeze(sum(sum(abs(fft(b.s,[],1)).^2,2),3));
ub=squeeze(sum(sum(abs(fft(b.b,[],1)).^2,2),3));
clf;h=loglog(uz(1:end/2));
hold('on');h(2)=loglog(ub(1:end/2));
legend('h',{'s','b'});
