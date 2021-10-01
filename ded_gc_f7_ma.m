%for a in $(seq 1 1000) ; do rsync -vap --progress asahi:gc/f7/ma/00100/27100 ~/gc/f7/ma/00100 --exclude a --exclude ayz --exclude avar --exclude sev  --exclude final --exclude force --bwlimit 250 ; done
%for a in $(seq 1 1000) ; do rsync -vap --progress asahi:gc/f7/mb/01000/20986 ~/gc/f7/mb/01000 --exclude a --exclude ayz --exclude avar --exclude sev  --exclude final --exclude force  ; done

nm='gc/f7/ma/00100/27100';
nm='gc/f7/mb/01000/20986';
p=ded_read_param(nm);
a=ded_read_javrg(nm,'a');
c=ded_coord(nm);
b=ded_zgrid(a,2*c.Nz,[],[],[],[],p.H,1);

[a b]=ded_diff_int(a,p,parity,dx,calc_diff,calc_int,calc_NS,fd)

k=size(a.b,1);
w=ichebintw(k);
b0=w*a.b;
b1=1./(w*(a.b.*(1-a.b)))/sqrt(2*pi);
b1(b0<1e-2)=0;



f= @(p) (1-erf(p(2)*(b.z(:)-p(1))))/2;
dz=b.z(2)-b.z(1);
opt=optimset('display','none');
q=[1 10];
clf;preprint([1 2],8);colour_lines
for j=1:10:840
  ff = @(p) f(p)-b.b(:,j);
  P(j,1:2)=lsqnonlin(ff,[b0(j),b1(j)],[-1 0],[2 100],opt);
  rb(j)=sqrt(mean(ff(P(j,1:2)).^2));
  plot(f(P(j,1:2)),b.z,b.b(:,j),b.z);
  axis([-0.05 1.05 0 2]);
  drawnow;
  nm=sprintf('~/pics/sdns/b-%05u.eps',j);
  print('-depsc2',nm);
end



rg=min(find(c.Jx>=p.x6)):max(find(b0>=4e-2));
hb=P(:,1);
wb=1./P(:,2);
fb=find(rb<2e-3 & b0>1e-2);
fb=min(fb):max(fb);
fb=rg;

figure(1);clf;
subplot(3,1,1);
plot(c.Jx(fb),hb(fb));
axis([-inf inf 0 max(hb)*1.05]);
ylabel('hb');

subplot(3,1,2);
plot(c.Jx(fb),wb(fb));
axis([-inf inf 0 inf]);
ylabel('wb');

subplot(3,1,3);
plot(c.Jx(fb),rb(fb));
axis([-inf inf 0 inf]);
ylabel('err');

m0=w*a.u/p.H;
m1=((c.Jz.^1)'.*w)*(a.u-m0)/p.H;
m2=((c.Jz.^2)'.*w)*(a.u-m0)/p.H;
m3=((c.Jz.^3)'.*w)*(a.u-m0)/p.H;
dU=b.u(end,:)-b.u(1,:);
hu=3/2*m2./m1-p.H;
wu=2*m1.^2./max(1e-2,8*m3.*m1-9*m2.^2+6*m2.*m1*p.H-4*p.H^2*m1.^2);
ru=0*hu;
hu=zeros(c.Nx,1);
wu=zeros(c.Nx,1);
clf;preprint([1 2],8);colour_lines

for j=1:10:c.Nx
  f = @(q) ded_uerf(b.z(:),q(1),q(2),q(3),m0(j),p.H); % ded_uerf(y,h,w,dU,U,H,n)
  ff = @(q) f(q)-b.u(:,j); % ded_uerf(y,h,w,dU,U,H,n)
  q=[hu(j) wu(j) b.u(end,j)-b.u(1,j)];
  q=lsqnonlin(ff,q,[-1 0 -2*p.V*p.H],[2 100 0],opt);
  ru(j)=sqrt(mean(ff(q).^2));
  plot(f(q),b.z,b.u(:,j),b.z);
  hu(j)=q(1);
  wu(j)=1/q(2);
  dU(j)=q(3);
  axis([-2.1 .25 0 2]);
  drawnow;
  nm=sprintf('~/pics/sdns/u-%05u.eps',j);
  print('-depsc2',nm);
end

fu=find((ru<p.V/10) & (-dU>p.V/5) & (hu>p.H/20));
fu=min(fu):max(fu);
fu=rg;

figure(2);clf;
subplot(4,1,1);
plot(c.Jx(fu),hu(fu));
axis([-inf inf 0 max(hu)*1.05]);
ylabel('hu');

subplot(4,1,2);
plot(c.Jx(fu),wu(fu));
axis([-inf inf 0 inf]);
ylabel('wu');

subplot(4,1,3);
plot(c.Jx(fu),-dU(fu));
axis([-inf inf 0 inf]);
ylabel('dU');

subplot(4,1,4);
plot(c.Jx(fu),ru(fu));
axis([-inf inf 0 inf]);
ylabel('err');

figure(3);clf;
subplot(3,1,1);
h=plot(c.Jx(rg),hu(rg),c.Jx(rg),hb(rg));
axis([-inf inf 0 1.05]);
legend(h,{'u','b'});
ylabel('h');

subplot(3,1,2);
h=plot(c.Jx(rg),wu(rg),c.Jx(rg),wb(rg));
axis([-inf inf 0 max(wu(rg))*1.05]);
legend(h,{'u','b'});
ylabel('w');

subplot(3,1,3);
h=plot(c.Jx(rg),(hb(rg)-hu(rg))./wu(rg),c.Jx(rg),(hb(rg)-hu(rg))./wb(rg));
axis([-inf inf 0 inf]);
legend(h,{'u','b'});
ylabel('hb-hu');



subplot(2,1,1);
wf=find(c.Jx>1&c.Jx<8);
f=@(p)p(1)*sqrt(p(2)-c.Jx(wf))-wu(wf);
p=lsqnonlin(f,[max(wu(wf)),max(c.Jx(wf))]);
pp=polyfit(c.Jx(wf),wu(wf),2);
plot(c.Jx(wf),wu(wf),c.Jx(wf),p(1)*sqrt(max(0,p(2)-c.Jx(wf))),c.Jx(wf),polyval(pp,c.Jx(wf)));


subplot(2,1,2);
wf=find(c.Jx>1&c.Jx<8);
f=@(p)p(1)*sqrt(p(2)-c.Jx(wf))-wb(wf);
p=lsqnonlin(f,[max(wb(wf)),max(c.Jx(wf))]);
pp=polyfit(c.Jx(wf),wb(wf),1);
plot(c.Jx(wf),wb(wf),c.Jx(wf),0*p(1)*sqrt(max(0,p(2)-c.Jx(wf))),c.Jx(wf),polyval(pp,c.Jx(wf)));

           

parity=ded_read_hdf([ded_dedalus_data_dir '/' nm '/' 'parity.hdf5']);
a=ded_read_javrg(nm,'a',trg,'combine');
nms=intersect(fieldnames(a),nms);
da=ded_diff(a,p,parity,c,nms);
Ia=ded_int(a,p,parity,c,nms);
