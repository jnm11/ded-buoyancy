function  a=ded_gc_fit_erf_profiles_a(a,p,c,parity)


z=c.Jz(:);
x=c.Jx(:);
Nx=c.NJx;
Nz=c.NJz;
H=p.H;
V=p.V;

k=size(a.b,1);
w=ichebintw(k);

ws=sqrt(w(:));
b0=w*a.b;
bz=(z(:)'.*w)*a.b;
a.c='zx';
ad=ded_diff(a,p,parity,c,{'b','u','w'});
dbs=w*(ad.dbdx.^2+ad.dbdz.^2);
dus=w*(ad.dudx.^2+ad.dwdz.^2 + (ad.dudz+ad.dwdx).^2/2);
uu2=w*a.u.^2;
b1=1./(w*(a.b.*(1-a.b)))/sqrt(2*pi);
b1(b0<1e-2)=0;


opt=optimset('display','none');
q=[1 10];

figs=isfield(qqq,'sz');
if figs;
  clf;preprint(qqq.sz,qqq.fnt);colour_lines
end

f1=max(find(b0>0.05*max(b0)));
x1=x(f1);x2=x(f1+1);
ba=b0(f1);bb=b0(f1+1);
X=((x2-x1)*0.05+x1*bb-x2*ba)/(bb-ba);

r=struct;
f= @(p,z) p(3)*serf(z,p(1),p(2));

for j=1:Nx
  ff = @(p) (f(p,z)-a.b(:,j)).*ws;
  [P p.rb(j)]=lsqnonlin(ff,[b0(j),b0(j)/2,1],[0 0 0],[p.H 100 1],opt);
  r.bh(j)=P(1);
  r.bw(j)=P(2);
  r.B1(j)=f(P,0);
  r.B2(j)=f(P,H);
  if mod(j,10)==0 & figs
    plot(a.b(:,j),z,f(P,z),z,'--');
    axis([-0.10 1.1 0 p.H]);
    text(1,p.H*0.95,sprintf('x=%5.2f',X-x(j)),'horizontalalignment','right','verticalalignment','top');
    set(gca,'ygrid','on','xgrid','on','box','on');
    drawnow;
    fnm=sprintf('%s-b-%05u.eps',nn,j);
    print('-depsc2',fnm);
  end
end
minu=min(a.u(:));
maxu=max(a.u(:));
rgu=[minu maxu]+(maxu-minu)/10*[-1 1];

for j=1:Nx
  u=a.u(:,j);
  f= @(p,z) w*u/H+p(3)*(serf(z,p(1),p(2))-1/H);
  ff = @(p) (f(p,z)-u).*ws;
  [P r.ru(j)]=lsqnonlin(ff,[1 1 1],[0 0 0],[],opt);
  r.uh(j)=P(1);
  r.uw(j)=P(2);
  r.U1(j)=f(P,0);
  r.U2(j)=f(P,H);
  if mod(j,10)==0 & figs
    plot(u,z,f(P,z),z,'--');
    axis([rgu 0 p.H]);
    text(rgu(2)*0.95+rgu(1)*0.05,p.H*0.95,sprintf('x=%5.2f',X-x(j)),'horizontalalignment','right','verticalalignment','top');
    set(gca,'ygrid','on','xgrid','on','box','on');
    drawnow;
    fnm=sprintf('%s-u-%05u.eps',nn,j);
    print('-depsc2',fnm);
  end
end

r.Tuu=w*(a.u.*a.u);
r.Tuw=w*(a.u.*a.w);
r.Tww=w*(a.w.*a.w);
r.Tbu=w*(a.b.*a.u);
r.Tbw=w*(a.b.*a.w);
r.Tbb=w*(a.b.*a.b);
r.Fuu=w*(a.uu-a.u.*a.u);
r.Fuw=w*(a.uw-a.u.*a.w);
r.Fww=w*(a.ww-a.w.*a.w);
r.Fbu=w*(a.bu-a.b.*a.u);
r.Fbw=w*(a.bw-a.b.*a.w);
r.Fbb=w*(a.bb-a.b.*a.b);

if false
P=zeros(Nx,3);
for k=1:3
  switch(k)
    case 1
      uw=a.uw-a.u.*a.w;
      mm='uw';
    case 2
      uw=a.uu-a.u.*a.u;
      mm='uu';
    case 3
      uw=a.ww-a.w.*a.w;
      mm='ww';
  end
  minu=min(uw(:));
  maxu=max(uw(:));
  rgu=[minu maxu]+(maxu-minu)/10*[-1 1];
  P=repmat(NaN,Nx,3);
  for j=1:Nx
    u=uw(:,j);
    p0=[r.uh(j) 1/r.uw(j) u(findmax(abs(u)))];
    if any(~isfinite(p0))
      continue;
    end
    ff = @(p) (f(p)-u).*ws;
    P(j,:)=lsqnonlin(ff,p0,[],[],opt);
    rb(j)=sqrt(mean(ff(P(j,:)).^2));
    if mod(j,10)==0 & figs
      plot(uw(:,j),z,f(P(j,:)),z,'--');
      axis([rgu 0 p.H]);
      text(rgu(2)*0.95+rgu(1)*0.05,p.H*0.95,sprintf('x=%5.2f',X-x(j)),'horizontalalignment','right','verticalalignment','top');
      set(gca,'ygrid','on','xgrid','on','box','on');
      drawnow;
      fnm=sprintf('%s-%s-%05u.eps',nn,mm,j);
      print('-depsc2',fnm);
    end
  end
  P(x>=X,:)=NaN;
  
  switch(k)
    case 1
      r.Fuw=w*uw;
      r.uwh=P(:,1)';
      r.uww=1./P(:,2)';
      r.ruw=rb(:)';
      r.Auw=P(:,3)';
    case 2
      r.Fuu=w*uw;
      r.uuh=P(:,1)';
      r.uuw=1./P(:,2)';
      r.ruu=rb(:)';
      r.Auu=P(:,3)';
    case 3
      r.Fww=w*uw;
      r.wwh=P(:,1)';
      r.www=1./P(:,2)';
      r.rww=rb(:)';
      r.Aww=P(:,3)';
  end
end
end

[w1 w2]=ichebendsw(c.Nz);
f=fieldnames(a);
for j=1:length(f)
  if all(size(a.(f{j}))>1)
    r.([f{j} '1']) = w1*a.(f{j});
    r.([f{j} '2']) = w2*a.(f{j});
    r.([f{j} 'I']) = w*a.(f{j});
  end
end
r.x=x(:)';
r.nm=nm;
r.X=X;

r.u1=r.U1;
r.u2=r.U2;
r.b1=r.B1;
r.b2=r.B2;

r.U1=w1*a.u;
r.U2=w2*a.u;
r.B1=w1*a.b;
r.B2=w2*a.b;

r.bz=bz;
r.dbs=dbs;
r.dus=dus;
r.uu2=uu2;

a=r;

[dd fn]=fileparts(fnmat);
if ~isdir(dd); mkdir(dd); end



save(fnmat,'a');




return;

plot(x,r.uh,x,r.bh);axis([0 p.L 0 p.H]);



fb=find(rb<2e-3 & b0>1e-2);
fb=min(fb):max(fb);
fb=rg;

figure(1);clf;
subplot(3,1,1);
plot(c.Jx(fb),bh(fb));
axis([-inf inf 0 max(bh).*1.05]);
ylabel('bh');

subplot(3,1,2);
plot(c.Jx(fb),bw(fb));
axis([-inf inf 0 inf]);
ylabel('bw');

subplot(3,1,3);
plot(c.Jx(fb),rb(fb));
axis([-inf inf 0 inf]);
ylabel('err');

m0=w.*a.u./p.H;
m1=((c.Jz.^1)'.*w).*(a.u-m0)./p.H;
m2=((c.Jz.^2)'.*w).*(a.u-m0)./p.H;
m3=((c.Jz.^3)'.*w).*(a.u-m0)./p.H;
dU=b.u(end,:)-b.u(1,:);
uh=3./2.*m2./m1-p.H;
uw=2.*m1.^2./max(1e-2,8.*m3.*m1-9.*m2.^2+6.*m2.*m1.*p.H-4.*p.H.^2.*m1.^2);
ru=0.*uh;
uh=zeros(c.Nx,1);
uw=zeros(c.Nx,1);
clf;preprint([1 2],8);colour_lines

for j=1:c.Nx
  f = @(q) ded_uerf(b.z(:),q(1),q(2),q(3),m0(j),p.H); % ded_uerf(y,h,w,dU,U,H,n)
  ff = @(q) f(q)-b.u(:,j); % ded_uerf(y,h,w,dU,U,H,n)
  q=[uh(j) uw(j) b.u(end,j)-b.u(1,j)];
  q=lsqnonlin(ff,q,[-1 0 -2.*p.V.*p.H],[2 100 0],opt);
  ru(j)=sqrt(mean(ff(q).^2));
  plot(f(q),b.z,b.u(:,j),b.z);
  uh(j)=q(1);
  uw(j)=1./q(2);
  dU(j)=q(3);
  axis([-2.1 .25 0 2]);
  drawnow;
  nm=sprintf('~./pics./sdns./%s-u-%05u.eps',nn,j);
  print('-depsc2',nm);
end



fu=find((ru<p.V./10) & (-dU>p.V./5) & (uh>p.H./20));
fu=min(fu):max(fu);
fu=rg;

figure(2);clf;
subplot(4,1,1);
plot(c.Jx(fu),uh(fu));
axis([-inf inf 0 max(uh).*1.05]);
ylabel('uh');

subplot(4,1,2);
plot(c.Jx(fu),uw(fu));
axis([-inf inf 0 inf]);
ylabel('uw');

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
h=plot(c.Jx(rg),uh(rg),c.Jx(rg),bh(rg));
axis([-inf inf 0 1.05]);
legend(h,{'u','b'});
ylabel('h');

subplot(3,1,2);
h=plot(c.Jx(rg),uw(rg),c.Jx(rg),bw(rg));
axis([-inf inf 0 max(uw(rg)).*1.05]);
legend(h,{'u','b'});
ylabel('w');

subplot(3,1,3);
h=plot(c.Jx(rg),(bh(rg)-uh(rg))./uw(rg),c.Jx(rg),(bh(rg)-uh(rg))./bw(rg));
axis([-inf inf 0 inf]);
legend(h,{'u','b'});
ylabel('bh-uh');



subplot(2,1,1);
wf=find(c.Jx>1&c.Jx<8);
f=@(p)p(1).*sqrt(p(2)-c.Jx(wf))-uw(wf);
p=lsqnonlin(f,[max(uw(wf)),max(c.Jx(wf))]);
pp=polyfit(c.Jx(wf),uw(wf),2);
plot(c.Jx(wf),uw(wf),c.Jx(wf),p(1).*sqrt(max(0,p(2)-c.Jx(wf))),c.Jx(wf),polyval(pp,c.Jx(wf)));


subplot(2,1,2);
wf=find(c.Jx>1&c.Jx<8);
f=@(p)p(1).*sqrt(p(2)-c.Jx(wf))-bw(wf);
p=lsqnonlin(f,[max(bw(wf)),max(c.Jx(wf))]);
pp=polyfit(c.Jx(wf),bw(wf),1);
plot(c.Jx(wf),bw(wf),c.Jx(wf),0.*p(1).*sqrt(max(0,p(2)-c.Jx(wf))),c.Jx(wf),polyval(pp,c.Jx(wf)));

           

parity=ded_read_hdf([ded_dedalus_data_dir './' nm './' 'parity.hdf5']);
a=ded_read_javrg(nm,'a',trg,'combine');
nms=intersect(fieldnames(a),nms);
da=ded_diff(a,p,parity,c,nms);
Ia=ded_int(a,p,parity,c,nms);




fns=cellstr_ls('~/data/dedalus/results/gc/f6/g/*/*/*/[23]*/profile.mat');

for j=1:length(fns)
  load(fns{j})
  if isfield(a,'u0')
    disp([fns{j} ': u0']);
  end
end
rm -rf data/dedalus/results/gc/f6/g/1000/0400/2304/29400 data/dedalus/results/gc/f6/g/1000/1600/2304/23100 data/dedalus/results/gc/f6/g/2000/0400/2304/28600 data/dedalus/results/gc/f6/g/4000/2000/4608/21500 data/dedalus/results/gc/f6/g/8000/2000/4608/21300 