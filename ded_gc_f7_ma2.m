%for a in $(seq 1 1000) ; do rsync -vap --progress asahi:gc/f7/ma/00100/27100 ~/gc/f7/ma/00100 --exclude a --exclude ayz --exclude avar --exclude sev  --exclude final --exclude force --bwlimit 250 ; done
%for a in $(seq 1 1000) ; do rsync -vap --progress asahi:gc/f7/mb/01000/20986 ~/gc/f7/mb/01000 --exclude a --exclude ayz --exclude avar --exclude sev  --exclude final --exclude force  ; done
%rsync -vap hamilton:gc/f7/mb/01000/20986 ~/gc/f7/mb/01000
%rsync -vap hamilton:gc/f7/ma/00100/32700 ~/gc/f7/ma/00100
nm='gc/f7/ma/00100/32700';

T=ded_convergence_T(nm);
if ~(T<0.8*b.t(end)) T = 0.8*b.t(end);end;
T=400;
trg=[T inf];
p=ded_read_param(nm);
c=ded_coord(nm);
pr=ded_read_hdf([ded_dedalus_data_dir '/' nm '/' 'parity.hdf5']);
a=ded_read_javrg(nm,'ay',trg,'combine');


fn={'b','d','p','u','w','uu','ww','uw','bu','bw'};
da=ded_diff(a,p,pr,c,fn);
Ia=ded_int(a,p,pr,c,{'b','d','p','u','w','uu','ww','uw','bu','bw'});
% $$$ dudx = pr_diff(a.u,c.dJx(1),2,1,[],-1);
% $$$ disp(sum(sum(dudx.*a.dudx))/sum(a.dudx(:).^2));
% $$$ disp(sum(sum(dudx.*da.dudx))/sum(a.dudx(:).^2));
da.ddbdxx   = pr_diff(a.b,c.dJx(1),2,2,[], 1);
da.ddudxx   = pr_diff(a.u,c.dJx(1),2,2,[],-1);
da.ddwdxx   = pr_diff(a.w,c.dJx(1),2,2,[], 1);
da.ddbdzz   = ichebdifff(a.b,1,p.H,2); % Assume spacing is Chebychev points
da.ddudzz   = ichebdifff(a.u,1,p.H,2); % Assume spacing is Chebychev points
da.ddwdzz   = ichebdifff(a.w,1,p.H,2); % Assume spacing is Chebychev points
Ia.dbwdzIx  = pr_int(da.dbwdz,c.dJx(1),2,[], +1);
Ia.duwdzIx  = pr_int(da.duwdz,c.dJx(1),2,[], -1);
Ia.wdudzIx  = pr_int( a.wdudz,c.dJx(1),2,[], -1);
Ia.wdwdzIx  = pr_int( a.wdwdz,c.dJx(1),2,[], +1);
Ia.ddbdzzIx = pr_int(da.ddbdzz,c.dJx(1),2,[],+1);
Ia.ddudzzIx = pr_int(da.ddudzz,c.dJx(1),2,[],-1);
Ia.ddwdzzIx = pr_int(da.ddwdzz,c.dJx(1),2,[],+1);
Ia.ddudxxIz = ichebintf(da.ddudxx,1,c.Jz,p.H);
Ia.ddwdxxIz = ichebintf(da.ddwdxx,1,c.Jz,p.H);
Ia.duwdxIz  = ichebintf(da.duwdx, 1,c.Jz,p.H);
w=ichebintw(c.NJz);
rg=find(c.Jx>2*p.x5 & c.Jx<p.L-2*p.x5);
lbl={'tot','uu','uw','p','uxx','uzz'};
wwuu=(w*a.ww)./(w*a.uu);
ded_eq_rms(rg, cat(3, a.ududx, a.wdudz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl,-1,0.01,2);
ded_eq_rms(rg ,cat(3,da.duudx,da.duwdz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl,-1,0.01,2);

ded_eq_rms(rg, cat(3,a.uu/2-a.uu(:,rg(end))/2,Ia.wdudzIx,a.p-a.p(:,rg(end)),-a.dudx/p.Re,-Ia.ddudzzIx/p.Re),lbl);
rg2=find(c.Jx>2 & c.Jx<3);
ded_eq_rms(rg2,cat(3, a.ududx, a.wdudz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl,-1,0.01,2);
ded_eq_rms(rg2,cat(3,da.duudx,da.duwdz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl,-1,0.01,2);

ded_eq_rms(rg, cat(3,a.ududx,a.wdudz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl,-1,c.dJx,2);
aa=cat(3,a.uu/2-a.uu(:,rg(end))/2,Ia.wdudzIx,a.p-a.p(:,rg(end)),-a.dudx/p.Re,-Ia.ddudzzIx/p.Re);aa=aa-aa(:,rg(end),:);
ded_eq_rms(rg,aa,lbl);



if false
for j=1:size(a.u,1)
  rg = max(find(a.u(j,:)>0));
  if isempty(rg)
    break;
  end
  rg = rg+(round(-0.5/c.dJx):round(0.5/c.dJx));
  fu = @(p) p(1)+p(2)*max(0,c.Jx(rg)'-p(3)).^(1/3);
  ff = @(p) fu(p)-a.u(j,rg);
  pp(j,:) = lsqnonlin(ff,[max(a.u(j,rg)) min(a.u(j,rg))-max(a.u(j,rg)) mean(c.Jx(rg))]);
  plot(c.Jx(rg),a.u(j,rg),c.Jx(rg),fu(pp(j,:)));
  axis('tight');
  drawnow;
end
end


f=rg(findmax(a.u(1,rg)));

fx=f+(-20:20);
fz=1:51;
u=a.u(fz,fx);
p=a.p(fz,fx);
w=a.w(fz,fx);
b=a.b(fz,fx);
psi1 =  Ia.uIz(fz,fx); psi1=psi1-min(psi1(:));
psi2 = -Ia.wIx(fz,fx); psi2=psi2-min(psi2(:));
z=c.Jz(fz);x=c.Jx(fx)';

P=a.p(fz,fx)+u.^2/2+w.^2/2+p.g*Ia.bIz(fz,fx);P=P-min(P(:));
contour(x,z,P);


[cu nu]=poly_fit2(z,x,u,3,3);
[cw nw]=poly_fit2(z,x,w,3,3);


x = linspace(-2,2,1e3);
clear(r);
for j=1:length(x)
  r(j,1:3)=roots([1 0 -x(j) 1e-2]);
end
r(imag(r)~=0)=NaN;
plot(x,r);

[xx zz]=ndgrid(x,z);

contour(x,z,psi1)
contour(x,z,psi2)

v=a.w(1:21,f-10:f+10);




plot(c.Jx,(w*a.ww)./(w*a.uu));


plot(w*a.b);

x=c.Jx(rg);
Pc=a.uu/2+a.ww/2+a.p;
Pd=-(a.dudx+Ia.ddudzzIx)/p.Re;
P=Pc+Pd;
figure(1);clf;
subplot(4,1,1);h=plot(x,P([1 end],rg)-P([1 end],rg(end)));legend(h,{'bottom','top'});ylabel('total')
title('Integrated x momentum');
subplot(4,1,2);h=plot(x,Pc([1 end],rg)-Pc([1 end],rg(end)));legend(h,{'bottom','top'});ylabel('P')
subplot(4,1,3);h=plot(x,Pd([1 end],rg)-Pd([1 end],rg(end)));legend(h,{'bottom','top'});ylabel('visc')
subplot(4,1,4);h=plot(x,a.SR([1 end],rg));legend(h,{'bottom','top'});ylabel('disp')
xlabel('x');

figure(2);clf;
Pc=da.duudx/2+da.dwwdx/2+da.dpdx;
Pd=-(da.ddudxx+da.ddudzz)/p.Re;


P=Pc+Pd;
subplot(4,1,1);h=plot(x,P([1 end],rg)-P([1 end],rg(end)));legend(h,{'bottom','top'});ylabel('total')
title('x momentum');
subplot(4,1,2);h=plot(x,Pc([1 end],rg)-Pc([1 end],rg(end)));legend(h,{'bottom','top'});ylabel('P')
subplot(4,1,3);h=plot(x,Pd([1 end],rg)-Pd([1 end],rg(end)));legend(h,{'bottom','top'});ylabel('visc')
subplot(4,1,4);h=plot(x,a.SR([1 end],rg));legend(h,{'bottom','top'});ylabel('disp')
xlabel('x');

ded_eq_rms(rg,cat(3,a.ududx,a.wdudz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl,-1,0.01,2);



w=ichebintw(c.NJz);
rms   = @(p) sqrt(mean(p(:,rg).^2));
rmsw  = @(p) sqrt(mean(w*p(:,rg).^2));
rmstb = @(p) sqrt(mean(p([1 end],rg).^2,2));
p.Peb=p.Re*p.Scb;
Db  = da.dbudx   + da.dbwdz                          - da.ddbdxx/p.Peb - da.ddbdzz/p.Peb;
Ib  = a.bu       + Ia.dbwdzIx                        - da.dbdx    /p.Peb - Ia.ddbdzzIx/p.Peb;
Dx  = da.duudx   + da.duwdz   + da.dpdx              - da.ddudxx  /p.Re   - da.ddudzz  /p.Re;
Ix  = a.uu       + Ia.duwdzIx + a.p                  -    a.dudx  /p.Re   - Ia.ddudzzIx/p.Re;Ix=Ix-Ix(:,rg(end));
Dx  = a.ududx    + a.wdudz    + da.dpdx              - da.ddudxx  /p.Re   - da.ddudzz  /p.Re;
Ix  = a.uu/2     + Ia.wdudzIx + a.p                  -    a.dudx  /p.Re   - Ia.ddudzzIx/p.Re;Ix=Ix-Ix(:,rg(end));
Dz  = da.duwdx   + da.dwwdz   + da.dpdz + p.g*a.b    - da.ddwdxx  /p.Re   - da.ddwdzz  /p.Re;
Iz  = Ia.duwdxIz + a.ww       + a.p     + p.g*Ia.bIz - Ia.ddwdxxIz/p.Re   -  a.dwdz    /p.Re;
rDb = [rmsw(da.dbudx)     rmsw(da.dbwdz)    rmsw(da.ddbdxx/p.Peb)  rmsw(da.ddbdzz/p.Peb)];
rIb = [rmsw(a.bu)         rmsw(Ia.dbwdzIx)  rmsw(da.dbdx/p.Peb)    rmsw(Ia.ddbdzzIx/p.Peb)];
rDx = [rmsw(da.ududx)     rmsw(a.wdudz)     rmsw(da.dpdx)     rmsw(da.ddudxx/p.Re)   rmsw(da.ddudzz/p.Re)];
wDb = [rms(w*da.dbudx)    rms(w*da.dbwdz)   rms(w*da.ddbdxx/p.Peb) rms(w*da.ddbdzz/p.Scb)];
wIb = [rms(w*a.bu)        rms(w*Ia.dbwdzIx) rms(w*da.dbdx/p.Peb)   rms(w*Ia.ddbdzzIx/p.Peb)];
wDx = [ rms(w*da.duudx/2) rms(w*da.dpdx)    rms(w*da.ddudxx/p.Re)  rms(w*da.ddudzz/p.Re)];
wIx = [ rms(w*a.uu/2)     rms(w*a.p)        rms(w*da.dudx/p.Re)    rms(w*Ia.ddudzzIx/p.Re)];
aa=cat(3,a.uu,Ia.duwdzIx,a.p,-a.dudx/p.Re,-Ia.dduxxIx/p.Re,-Ia.ddudzzIx/p.Re);
aa=cat(3,da.duudx,da.duwdz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re);

rIx1  = [rmsw(a.uu)       rmsw(Ia.duwdzIx) rmsw(a.p)     rmsw(a.dudx/p.Re)    rmsw(Ia.ddudzzIx/p.Re)];
rIx2  = [rmsw(a.uu/2)     rmsw(Ia.wdudzIx) rmsw(a.p)     rmsw(a.dudx/p.Re)    rmsw(Ia.ddudzzIx/p.Re)];
rDx1  = [rmsw(da.duudx)   rmsw(da.duwdz)   rmsw(da.dpdx) rmsw(da.ddudxx/p.Re) rmsw(da.ddudzz  /p.Re)];
rDx2  = [rmsw(da.duudx/2) rmsw(a.wdudz)    rmsw(da.dpdx) rmsw(da.ddudxx/p.Re) rmsw(da.ddudzz  /p.Re)];


ded_eq_rms(rg,cat(3,da.duudx,da.duwdz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl);

ded_eq_rms(rg,cat(3,a.ududx,a.wdudz,da.dpdx,-da.ddudxx/p.Re,-da.ddudzz/p.Re),lbl);

plot(x,Ix(:,rg))
 
% D0 d(bu)/dx and kappa*d^2(b)/dx^2 are important
% D1 d(bu)/dx and d(bw)/dz          are important
% I0   bu     and kappa*d(b)/dx     are important
% I1   bu     and int(d(bw)/dz,dx)  are important
z=c.Jz'-mean(c.Jz);
w=ichebintw(c.NJz);
D0=[rms(w*da.dbudx)  rms(w*da.dbwdz)   [rms(w*da.ddbdxx) rms(w*da.ddbdzz)]/(p.Re*p.Scb) ];
I0=[rms(w*a.bu)    rms(w*Ia.dbwdzIx)   [rms(w*da.dbdx) rms(w*Ia.ddbdzzIx)]/(p.Re*p.Scb) ];
w=w.*c.Jz';
D1=[rms(w*da.dbudx)  rms(w*da.dbwdz)   [rms(w*da.ddbdxx) rms(w*da.ddbdzz)]/(p.Re*p.Scb) ];
I1=[rms(w*a.bu)    rms(w*Ia.dbwdzIx)   [rms(w*da.dbdx) rms(w*Ia.ddbdzzIx)]/(p.Re*p.Scb) ];
w=w.*c.Jz';
D2=[rms(w*da.dbudx)  rms(w*da.dbwdz)   [rms(w*da.ddbdxx) rms(w*da.ddbdzz)]/(p.Re*p.Scb) ];
I2=[rms(w*a.bu)    rms(w*Ia.dbwdzIx)   [rms(w*da.dbdx) rms(w*Ia.ddbdzzIx)]/(p.Re*p.Scb) ];
figure(3);clf;
subplot(4,1,1);plot(x,w*Db(:,rg));
subplot(4,1,2);plot(x,w*Ib(:,rg));


z=c.Jz'-mean(c.Jz);
Dx0= [rms(w*da.duudx/2) rms(w*da.dpdx) rms(w*da.ddudxx/p.Re) rms(w*da.ddudzz/p.Re)];
Ix0= [rms(w*a.uu/2)     rms(w*a.p)     rms(w*da.dudx/p.Re)   rms(w*Ia.ddudzzIx/p.Re)];
w=w.*c.Jz';
Dx1= [rms(w*da.duudx/2) rms(w*da.dpdx) rms(w*da.ddudxx/p.Re) rms(w*da.ddudzz/p.Re)];
Ix1= [rms(w*a.uu/2)     rms(w*a.p)     rms(w*da.dudx/p.Re)   rms(w*Ia.ddudzzIx/p.Re)];
w=w.*c.Jz';
Dx2= [rms(w*da.duudx/2) rms(w*da.dpdx) rms(w*da.ddudxx/p.Re) rms(w*da.ddudzz/p.Re)];
Ix2= [rms(w*a.uu/2)     rms(w*a.p)     rms(w*da.dudx/p.Re)   rms(w*Ia.ddudzzIx/p.Re)];


w=ichebintw(c.NJz);
plot(x,a.p(1,rg)-a.p(end,rg)-w*a.b(:,rg))

Dz = da.duwdx + da.dwwdz + da.dpdz + p.g*a.b - da.ddwdxx/p.Re - da.ddwdzz/p.Re;
Iz = Ia.duwdxIz + a.ww   +  a.p    + p.g*Ia.bIz - Ia.ddwdxxIz/p.Re - a.dwdz/p.Re;
Iz=Iz-Iz(end,:);
subplot(4,1,1);plot(c.Jz,Dz(:,rg));
subplot(4,1,2);plot(c.Jz,Iz(:,rg));
subplot(4,1,3);plot(x,w*Dz(:,rg));
subplot(4,1,4);plot(x,w*Iz(:,rg));



subplot(2,1,1);plot(c.Jx,dudx(end,:),c.Jx,a.dudx(end,:));
subplot(2,1,2);plot(c.Jx,dudx(end,:)-a.dudx(end,:));

plot(c.Jx,dudx,c.Jx,a.dudx)
a.dudxx = pr_diff(da.dudx,c.dd(1),2,1,[],1);
b=ded_diff(da,p,pr,c,{'dudx','dudz','dwdx','dudz'});
%a.duudz=da.

rg=find(v.Jx>p.x5 & v.Jx<p.L-p.x5);
rg=find(c.Jx>p.x5 & c.Jx<8.2);
x=c.Jx(rg);
subplot(4,1,1);
h=plot(x,da.duudx(1,rg),x,da.duwdz(1,rg),x,da.dpdx(1,rg));legend(h,{'uu','uw','p'});
subplot(4,1,2);
plot(x,da.duudx(1,rg)+da.duwdz(1,rg)+da.dpdx(1,rg));
subplot(4,1,3);
h=plot(x,da.ududx(1,rg),x,da.wdudz(1,rg),x,da.dpdx(1,rg));legend(h,{'uu','uw','p'});

subplot(2,1,1);
plot(x,da.duudx(1,rg)/2+da.dpdx(1,rg));


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
for j=1:c.Nx
  ff = @(p) f(p)-b.b(:,j);
  P(j,1:2)=lsqnonlin(ff,[b0(j),b1(j)],[-1 0],[2 100],opt);
  rb(j)=sqrt(mean(ff(P(j,1:2)).^2));
  plot(b.z,f(P(j,1:2)),b.z,b.b(:,j));
  drawnow;
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

for j=1:c.Nx
  f = @(q) ded_uerf(b.z(:),q(1),q(2),q(3),m0(j),p.H); % ded_uerf(y,h,w,dU,U,H,n)
  ff = @(q) f(q)-b.u(:,j); % ded_uerf(y,h,w,dU,U,H,n)
  q=[hu(j) wu(j) b.u(end,j)-b.u(1,j)];
  q=lsqnonlin(ff,q,[-1 0 -2*p.V*p.H],[2 100 0],opt);
  ru(j)=sqrt(mean(ff(q).^2));
  plot(b.z,f(q),b.z,b.u(:,j));
  hu(j)=q(1);
  wu(j)=1/q(2);
  dU(j)=q(3);
  drawnow;
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

           
g=p.g;
w=ichebintw(c.NJz);
[w1 w2]=ichebendsw(c.NJz);
dw=w2-w1;
subplot(2,1,1);
plot(c.Jx,g*w*Ia.bIz,c.Jx,-dw*a.p,c.Jx,dw*a.dwdz/p.Re,c.Jx,dw*Ia.ddwdxxIz/p.Re,c.Jx,dw*Ia.duwdxIz);
subplot(2,1,2);
dp=a.p(1,:)-a.p(end,:);
plot(c.Jx,dw*a.p+dw*Ia.duwdxIz);
plot(c.Jx,-dp+dw*Ia.duwdxIz);


f= find(sum(a.d.^2)<1e-3);f=f(1)+3:f(end)-3;
x=c.Jx(f);
subplot(5,1,1);
plot(x,w2*da.dpdx(:,f)+w2*da.duudx(:,f)/2-w2*da.ddudxx(:,f)/p.Re-w2*da.ddudzz(:,f)/p.Re);
hold('on');
plot(x,w2*da.dpdx(:,f)+w2*da.duudx(:,f)/2);
title('x momentum top');
subplot(5,1,2);
plot(x,w1*da.dpdx(:,f)+w1*da.duudx(:,f)/2-w1*da.ddudxx(:,f)/p.Re-w1*da.ddudzz(:,f)/p.Re);
title('x momentum bottom');
subplot(5,1,3);
plot(x,w*da.dpdx(:,f)+w*da.duudx(:,f));
title('x momentum integrated');


Idpdx   = ichebintf(da.dpdx, 1,c.Jz,p.H);
Iduudx  = ichebintf(da.duudx,1,c.Jz,p.H);
Iduwdz  =            a.uw;
Iddudxx = ichebintf(da.ddudxx,1,c.Jz,p.H);
Iddudzz =            a.dudz;

IIdpdx   = ichebintf(Idpdx, 1,c.Jz,p.H);
IIduudx  = ichebintf(Iduudx,1,c.Jz,p.H);
IIduwdz  = ichebintf(Iduwdz,1,c.Jz,p.H);
IIddudxx = ichebintf(Iddudxx,1,c.Jz,p.H);
IIddudzz =            a.u;

subplot(5,1,4);
plot(x,w*(Idpdx(:,f)+Iduudx(:,f)+Iduwdz(:,f)+ Iddudxx(:,f)/p.Re+Iddudzz(:,f)/p.Re));

Ex   = da.duudx + da.duwdz + da.dpdx - (da.ddudxx+da.ddudzz)/p.Re;
Ex   = Ex(:,f);
IEx  = ichebintf(Ex,1,c.Jz,p.H);
IIEx = ichebintf(IEx,1,c.Jz,p.H);
subplot(5,1,1);plot(x,w2*Ex);title('x momentum top');
subplot(5,1,2);plot(x,w1*Ex);title('x momentum bottom');
subplot(5,1,3);plot(x, w*Ex);title('x momentum I');
subplot(5,1,4);plot(x, w*IEx);title('x momentum II');
subplot(5,1,5);plot(x, w*IIEx);title('x momentum III');

Ey   = da.duwdx + da.dwwdz + da.dpdz - (da.ddwdxx+da.ddwdzz)/p.Re +p.g*a.b;
Ey   = Ey(:,f);
IEy  = ichebintf(Ey,1,c.Jz,p.H);
IIEy = ichebintf(IEy,1,c.Jz,p.H);
subplot(5,1,1);plot(x,w2*Ey);title('y momentum top');
subplot(5,1,2);plot(x,w1*Ey);title('y momentum bottom');
subplot(5,1,3);plot(x, w*Ey);title('y momentum I');
subplot(5,1,4);plot(x, w*IEy);title('y momentum II');
subplot(5,1,5);plot(x, w*IIEy);title('y momentum III');

Eb  = da.dbudx + da.dbwdz  - (da.ddbdxx+da.ddbdzz)/p.Re/p.Scb;
Eb  = Eb(:,f);
Eb1 = ichebintf(Eb,1,c.Jz,p.H);
Eb2 = ichebintf(Eb1,1,c.Jz,p.H);
IEb   = w*Eb;
IIEb  = w*Eb1;
IIIEb = w*Eb2;

IEb =  ichebintf(da.dbudx-da.ddbdxx/p.Re/p.Scb,1,c.Jz,p.H);IEb=IEb(:,f);
subplot(5,1,1);plot(x,w2*Eb);title('b conservation top');
subplot(5,1,2);plot(x,w1*Eb);title('b conservation bottom');
subplot(5,1,3);plot(x,IEb);title('b conservation I');
subplot(5,1,4);plot(x,IIEb);title('b conservation II');
subplot(5,1,5);plot(x,IIIEb);title('b conservation III');
