function ded_gc_shw(nm,cons,tol,rgx)
% Test shallow water equations
%a=data/dedalus/results/gc/f6/mg/2000/0200/2304/41105/ ; rsync -uvap hamilton:$a ~/$a --delete
%a=data/dedalus/results/gc/f6/mg/2000/0800/2304/23645/ ; rsync -uvap hamilton:$a ~/$a --delete
%ded_gc_shw('gc/f6/mg/2000/0200/2304/41105',true,1e-3,[10 20])
%ded_gc_shw('gc/f6/mg/2000/0800/2304/23645',true,1e-3,[10 20])
%nm='gc/f6/mg/2000/0200/2304/41105';
%nm='gc/f6/mg/2000/0800/2304/23645';
fn=[ded_dedalus_data_dir '/results/' nm];

pm=load([fn '/param.mat']);pm=pm.p;
pr=load([fn '/profile-ray.mat']);pr=pr.a;
pa=load([fn '/parity.mat']);pa=pa.p;
c=load([fn '/coord.mat']);c=c.c;
a=load([fn '/ay.mat']);a=a.a;

pa.duudx=pa.ududx;
pa.duwdx=pa.udwdx;
pa.dwwdx=pa.wdwdx;
pa.duudz=pa.ududz;
pa.duwdz=pa.udwdz;
pa.dwwdz=pa.wdwdz;
pa.ddudxdx=pa.u;
pa.ddudzdz=pa.u;
pa.ddwdxdx=pa.w;
pa.ddwdzdz=pa.w;

g=pm.g;
V=pm.V;
H=pm.H;
F=V./sqrt(H*g*pr.b1);
Re=pm.Re;
Pe=pm.Re*pm.Scb;

w  = ichebintw(c.Nz);w=w(:);z=c.z;
w0 = w.*(z.^0)/H^1;
w1 = w.*(z.^1)/H^2;
w2 = w.*(z.^2)/H^3;
w3 = w.*(z.^3)/H^4;

% $$$ IEX1=U*H - (H/2 - uh).*U.^2 - uh./F.^2;
% $$$ plot(x,IEX1(f));
% $$$ plot(pr.x,(U*H - (H/2 - uh).*U.^2)./uh);

a=ded_add_diff(a,pm,pa,c,{'b','u','w','uw','p','uu','ww','bu','bw'});
a=ded_add_diff(a,pm,pa,c,{'dudx','dudz','dwdx','dwdz','dbdx','dbdz'});
a.ududx=a.duudx/2;
a.wdwdz=a.dwwdz/2;
a.udwdx=a.duwdx+a.wdwdz;
a.wdudz=a.duwdz+a.ududx;

a=ded_add_int(a,pm,pa,c,{'dbdx','b','wdudz','udwdx','duwdx','duwdz','ddudxdx','ddudzdz','ddwdxdx','ddwdzdz'});

a.Iddudzz=pr_int(a.ddudzdz,c.dx,2,1,1);


%ub0=(-w*H + H*uh - uh.^2).*U - H*uh;

f=find(pr.x>=rgx(1)&pr.x<=rgx(2));
x=pr.x(f);
z=c.z;
u1 = pr.u1(f);
b1 = pr.b1(f);
uh = pr.uh(f);
bh = pr.bh(f);
uw = pr.uw(f);
bw = pr.bw(f);

P1=a.p(  1,:)+a.u(  1,:).^2/2-a.dudx(  1,:)/Re-a.Iddudzz(  1,:)/Re; % Bernouii along top and bottom stream lines
P2=a.p(end,:)+a.u(end,:).^2/2-a.dudx(end,:)/Re-a.Iddudzz(end,:)/Re; 

if isempty(f)
  d=sqrt(sum(w0.*(a.dudx+a.dwdz).^2));
  dtol=1e-4;f=find(d<=dtol);f=min(f):max(f);
end
x=c.x(f);
z=c.z;

% Collect the terms for the Navier Stokes equations 
% 1 z direction
% 2 x direction
% 3 individual terms

nmb={'bu','bw','bxx','bzz'};

En=cat(3,a.ududx,a.udwdx,a.wdudz,a.wdwdz);En=En(:,f,:);
Ec=cat(3,a.duudx,a.duwdx,a.duwdz,a.dwwdz);Ec=Ec(:,f,:);
Eb=cat(3, +a.dbudx, +a.dbwdz, -a.ddbdxdx/Pe, -a.ddbdzdz/Pe            );
Ib=cat(3, +a.dbudx, +a.dbwdz, -a.ddbdxdx/Pe, -a.ddbdzdz/Pe            );

if cons  % Conservative
  Ix=cat(3,+a.p,       +a.uu, +a.Iduwdzdx,       -a.dudx/Re, -a.Iddudzdzdx/Re        );
  Iz=cat(3,+a.p, +a.Iduwdxdz,       +a.ww, -a.Iddwdxdxdz/Re,       -a.dwdz/Re, +g*a.Ibdz);
  %dpdx=cat(3, +a.Idduwdxdxdz,       +a.dwwdx,  +a.ddudxdx/Re, +g*a.Idbdxdz);  %-a.Idddwdxdxdxdz/Re  % Need to add on bottom
  Ex=cat(3, +a.dpdx, +a.duudx,   +a.duwdz, -a.ddudxdx/Re, -a.ddudzdz/Re            );
  %Ex=cat(3, +g*a.Idbdxdz, +a.duudx,   +a.duwdz, -a.ddudxdx/Re, -a.ddudzdz/Re            );
  Ez=cat(3, +a.dpdz, +a.duwdx,   +a.dwwdz, -a.ddwdxdx/Re, -a.ddwdzdz/Re,+g*a.b);
  nmx={'dpdx','duudx','duwdz','dudxx','dudzz'    };
  nmz={'dpdz','duwdx','dwwdz','dwdxx','dwdzz','b'};
else
  Ix=cat(3,+a.p,     +a.uu/2, +a.Iwdudzdx,       -a.dudx/Re, -a.Iddudzdzdx/Re        );
  Iz=cat(3,+a.p, +a.Iudwdxdz,     +a.ww/2, -a.Iddwdxdxdz/Re,       -a.dwdz/Re, +g*a.Ibdz);
  Ex=cat(3, +a.dpdx, +a.ududx, +a.wdudz, -a.ddudxdx/Re, -a.ddudzdz/Re            );
  Ez=cat(3, +a.dpdz, +a.udwdx, +a.wdwdz, -a.ddwdxdx/Re, -a.ddwdzdz/Re,+g*a.b);
  nmx={'dpdx','ududx','wdudz','dudxx','dudzz'    };
  nmz={'dpdz','udwdx','wdwdz','dwdxx','dwdzz','b'};
end
Ex=Ex(:,f,:);
Ez=Ez(:,f,:);
Eb=Eb(:,f,:);
Ix=Ix(:,f,:);
Iz=Iz(:,f,:);
Ib=Ib(:,f,:);
Ix=Ix-Ix(:,1,:);
Iz=Iz-Iz(1,:,:);

% $$$ figure(4);clf;
% $$$ subplot(2,1,1);h=semilogy(x,sqrt(squeeze(sum(w0.*En.^2,1))));legend(h,{'ududx','udwdx','wdudz','wdwdz'});axis('tight');set(h([1 3]),'color',[0 0 1]);
% $$$ subplot(2,1,2);h=semilogy(x,sqrt(squeeze(sum(w0.*Ec.^2,1))));legend(h,{'duudx','duwdx','duwdz','dwwdz'});axis('tight');set(h([1 3]),'color',[0 0 1]);



figure(1);clf;
subplot(3,1,1);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ex,3).^2,1))));ylabel('Nx');xlabel('x');axis('tight');
subplot(3,1,2);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ez,3).^2,1))));ylabel('Nz');xlabel('x');axis('tight');
subplot(3,1,3);semilogy(x,sqrt(squeeze(sum(w0.*sum(Eb,3).^2,1))));ylabel('Nb');xlabel('x');axis('tight');
subplot(3,1,3);title('Navier stokes rms residuals');

figure(5);clf;
subplot(3,1,1);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ix,3).^2,1))));ylabel('Ix');xlabel('x');axis('tight');
subplot(3,1,2);semilogy(x,sqrt(squeeze(sum(w0.*sum(Iz,3).^2,1))));ylabel('Iz');xlabel('x');axis('tight');
subplot(3,1,3);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ib,3).^2,1))));ylabel('Ib');xlabel('x');axis('tight');
subplot(3,1,3);title('Navier stokes rms residuals');

if false
  % Which terms are important in the NS equations
  figure(2);clf;
  subplot(3,1,1);h=semilogy(x,sqrt(squeeze(sum(w0.*Ex.^2,1))));legend(h,nmx);axis('tight');ylabel('Nx');
  subplot(3,1,2);h=semilogy(x,sqrt(squeeze(sum(w0.*Ez.^2,1))));legend(h,nmz);axis('tight');ylabel('Nz');
  subplot(3,1,3);h=semilogy(x,sqrt(squeeze(sum(w0.*Eb.^2,1))));legend(h,nmb);axis('tight');ylabel('Nb');
  subplot(3,1,1);title('Navier stokes terms rms');
end


ff=find(x>rgx(1)&x<rgx(2));

[px rx sx]=ded_find_significant(Ex(:,ff,:),w0);fx=min([length(rx) find(rx<=tol)]);ppx=px{fx};rrx=rx(fx);tx=cellsprintf('%s ',{nmx{ppx}});tx=[tx{:}];
[pz rz sz]=ded_find_significant(Ez(:,ff,:),w0);fz=min([length(rz) find(rz<=tol)]);ppz=pz{fz};rrz=rz(fz);tz=cellsprintf('%s ',{nmz{ppz}});tz=[tz{:}];
[pb rb sb]=ded_find_significant(Eb(:,ff,:),w0);fb=min([length(rb) find(rb<=tol)]);ppb=pb{fb};rrb=rb(fb);tb=cellsprintf('%s ',{nmb{ppb}});tb=[tb{:}];
disp(rx);disp(rz);disp(rb);
figure(3);clf;
subplot(3,1,1);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ex(:,:,ppx),3).^2,1))));axis('tight');ylabel('Nx');title(sprintf('%5.2e %s ',rrx,tx));
subplot(3,1,2);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ez(:,:,ppz),3).^2,1))));axis('tight');ylabel('Nz');title(sprintf('%5.2e %s ',rrz,tz));
subplot(3,1,3);semilogy(x,sqrt(squeeze(sum(w0.*sum(Eb(:,:,ppb),3).^2,1))));axis('tight');ylabel('Nb');title(sprintf('%5.2e %s ',rrb,tb));

[px rx sx]=ded_find_significant(Ix(:,ff,:),w0);fx=min([length(rx) find(rx<=tol)]);ppx=px{fx};rrx=rx(fx);tx=cellsprintf('%s ',{nmx{ppx}});tx=[tx{:}];
[pz rz sz]=ded_find_significant(Iz(:,ff,:),w0);fz=min([length(rz) find(rz<=tol)]);ppz=pz{fz};rrz=rz(fz);tz=cellsprintf('%s ',{nmz{ppz}});tz=[tz{:}];
[pb rb sb]=ded_find_significant(Ib(:,ff,:),w0);fb=min([length(rb) find(rb<=tol)]);ppb=pb{fb};rrb=rb(fb);tb=cellsprintf('%s ',{nmb{ppb}});tb=[tb{:}];
disp(rx);disp(rz);disp(rb);
figure(6);clf;
subplot(3,1,1);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ix(:,:,ppx),3).^2,1))));axis('tight');ylabel('Nx');title(sprintf('%5.2e %s ',rrx,tx));
subplot(3,1,2);semilogy(x,sqrt(squeeze(sum(w0.*sum(Iz(:,:,ppz),3).^2,1))));axis('tight');ylabel('Nz');title(sprintf('%5.2e %s ',rrz,tz));
subplot(3,1,3);semilogy(x,sqrt(squeeze(sum(w0.*sum(Ib(:,:,ppb),3).^2,1))));axis('tight');ylabel('Nb');title(sprintf('%5.2e %s ',rrb,tb));


spi=sqrt(pi);
s2=sqrt(2);
expu=exp(-uh.^2./uw.^2);
erfu=erf(uh./uw);
expb=exp(-bh.^2./bw.^2);
erfb=erf(bh./bw);
h=(uh+bh)/2;
U = H*(V + u1)./(erfu*H - uh);
B = b1./erf(bh./bw);
%uh=h;
%bh=h;



ww=sqrt(bw.^2 + uw.^2);
hd=(+bh - uh)./ww;
hm=(+bh + uh)./ww;
um0 = -V.*H.*uh.^0;
um1 = -H.^2.*V./2 - (-U.*expu.*uh.*uw + spi.*(-(uh.^2 + uw.^2./2).*U.*erfu + H.*U.*uh))./(2.*spi);
um2 = ((2.*uh.^3 + (-2.*H.^2 + 3.*uw.^2).*uh).*U)./6 - H.^3.*V./3;
um3 = -H.^4.*V./4 - U.*((-uh.^3.*uw - 5./2.*uh.*uw.^3).*expu + spi.*((-uh.^4 - 3.*uh.^2.*uw.^2 - 3./4.*uw.^4).*erfu + H.^3.*uh))./(4*spi);

bm0 = bh.*B;
bm1 = B.*(bh.*bw.*expb + erfb.*spi.*(bh.^2 + bw.^2./2))./(2.*spi);
bm2 = bh.*B.*(2.*bh.^2 + 3.*bw.^2)./6;
bm3 = B.*(((4.*bh.^3.*bw + 10.*bh.*bw.^3).*expb)./4 + erfb.*spi.*(bh.^4 + 3.*bh.^2.*bw.^2 + 3./4.*bw.^4))./(4.*spi);

bum0 = -B.*(V*bh + U.*bh.*uh./H) + B.*U.*ww.*(hm.*erf(hm) + exp(-hm.^2)./spi-hd.*erf(hd) - exp(-hd.^2)./spi)./2;
bum0 = -B.*(V*bh + U.*bh.*uh./H) + B.*U.*ww.*(hm -(1+hd.^2)/spi)./2;
uum0 =  U.^2.*s2.*uw.*exp(-2.*uh.^2./uw.^2)./(2.*spi) + U.^2.*uh.*erf(uh.*s2./uw) - U.^2.*s2.*uw./(2.*spi) - U.^2.*uh.^2./H + H.*V.^2;
uum0 =  U.*U.*(uh - s2.*uw./(2.*spi) - uh.^2./H) + H.*V.*V;
p=a.p(1,f);




pm0 = H.^1/1*p +g.*B.*(2.*bh.^2 + bw.^2).*erfb./4 + bh.*bw.*B.*g.*expb./(2.*spi) - g.*B.*H.*bh;
pm2 = H.^3/3*p +g.*B.*(4.*bh.^4 + 12.*bh.^2.*bw.^2 + 3.*bw.^4).*erfb./48 + g.*B.*bh.*bw.*(2.*bh.^2 + 5.*bw.^2).*expb./(24.*spi) - g.*B.*H.^3.*bh./3;
pm1 = H.^2/2*p +g.*B.*bh.*( -6.*H.^2 +2.*bh.^2  +3.*bw.^2)/12;
pm3 = H.^4/4*p +g.*B.*bh.*(-20.*H.^4 +4.*bh.^4 +15.*bw.^4 +20.*bh.^2.*bw.^2)./80;
pm0 = g*B*(2*bh^2 + bw^2)*erf(bh/bw)/4 + g*B*bh*bw*exp(-bh^2/bw^2)/(2*sqrt(Pi)) - H*(B*bh*g - P)


pm0 = g*B*(2*bh^2 + bw^2)*erf(bh/bw)/4 + g*B*bh*bw*exp(-bh^2/bw^2)/(2*sqrt(Pi)) - H*(B*bh*g - P)
pm1 = B*g*bh^3/6 - g*B*(H^2 - bw^2/2)*bh/2 + P*H^2/2;
pm2 = g*B*(4*bh^4 + 12*bh^2*bw^2 + 3*bw^4)*erf(bh/bw)/48 + g*B*bh*bw*(2*bh^2 + 5*bw^2)*exp(-bh^2/bw^2)/(24*sqrt(Pi)) - H^3*(B*bh*g - P)/3

pm3 = B*g*bh^5/20 + B*g*bh^3*bw^2/4 - ((H^4 - (3*bw^4)/4)*g*B*bh)/4 + P*H^4/4




figure
subplot(4,3, 1);plot(x,pr.um0(f)-um0);ylabel('u0');subplot(4,3, 2);plot(x,pr.bm0(f)-bm0);ylabel('b0');subplot(4,3, 3);plot(x,pr.pm0(f)-pm0);;ylabel('p0');
subplot(4,3, 4);plot(x,pr.um1(f)-um1);ylabel('u1');subplot(4,3, 5);plot(x,pr.bm1(f)-bm1);ylabel('b1');subplot(4,3, 6);plot(x,pr.pm1(f)-pm1);;ylabel('p1');
subplot(4,3, 7);plot(x,pr.um2(f)-um2);ylabel('u2');subplot(4,3, 8);plot(x,pr.bm2(f)-bm2);ylabel('b2');subplot(4,3, 9);plot(x,pr.pm2(f)-pm2);;ylabel('p2');
subplot(4,3,10);plot(x,pr.um3(f)-um3);ylabel('u3');subplot(4,3,11);plot(x,pr.bm3(f)-bm3);ylabel('b3');subplot(4,3,12);plot(x,pr.pm3(f)-pm3);;ylabel('p3');

keyboard
plot(x,pr.bum0(f),x,bum0);
plot(x,pr.uum0(f)-uum0);

Exm0 = uum0
plot(x,u1,x,a.u(1,f));


U(X).*H - (H./2 - h(X)).*U(X).^2 - h(X)./F.^2

F=V./sqrt(H.*g.*pr.b1);


P1=pr.u1.^2./2+pr.p1;
dwdz=a.dwdz(end,f)-a.dwdz(1,f);

wwm0=w'*a.ww(:,f);
plot(x,pr.pm0(f),x,-g*pr.bm1(f)+p + wwm0/2);

a.pm0=w'*(a.p-a.p(1,:));
pmm0b=-g*(H*bm0-bm1)-w'*a.ww;
bb=-g*(W0*a.b(:,f));
plot(x,a.pm0(f),x,bb);


