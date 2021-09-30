function ded_serf_moments(a,p,b)
% Check serf moment calculations

g  = p.g;
V  = p.V;
H  = p.H;
Re = p.Re;
Pe = p.Re*p.Scb;

if false
  w  = ichebintw(c.Nz);w=w(:);z=c.z;
  w0 = w.*(z.^0)/H^1;
  w1 = w.*(z.^1)/H^2;
  w2 = w.*(z.^2)/H^3;
  w3 = w.*(z.^3)/H^4;
end

u1 = a.u1(f);
b1 = a.b1(f);
uh = a.uh(f)/H;
bh = a.bh(f)/H;
uw = a.uw(f)/H;
bw = a.bw(f)/H;

U = H*(V + u1)./(erfu*H - uh)/V;
B = b1./erf(bh./bw);

ww=sqrt(bw.^2 + uw.^2);
hd=(+bh - uh)./ww;
hm=(+bh + uh)./ww;

P1=a.p(  1,:)+a.u(  1,:).^2/2-a.dudx(  1,:)/Re-a.Iddudzz(  1,:)/Re; % Bernouii along top and bottom stream lines
P2=a.p(end,:)+a.u(end,:).^2/2-a.dudx(end,:)/Re-a.Iddudzz(end,:)/Re; 


spi=sqrt(pi);
s2=sqrt(2);
expu  = exp(-uh.^2./uw.^2)/spi.*uw;
expb  = exp(-bh.^2./bw.^2)/spi.*bw;
expm  = exp(-hm.^2)/spi;
expd  = exp(-hd.^2)/spi;
expu2 = exp(-2*uh.^2./uw.^2)/spi;
expb2 = exp(-2*bh.^2./bw.^2)/spi;
expm2 = exp(-2*hm.^2)/spi;
expd2 = exp(-2*hd.^2)/spi;
erfu  = erf(uh./uw);
erfb  = erf(bh./bw);
erfm  = erf(hm);
erfd  = erf(hd);
erfu2 = erf(spi*uh./uw);
erfb2 = erf(spi*bh./bw);
erfm2 = erf(spi*hm);
erfd2 = erf(spi*hd);


um0 = -V+V*U*uh.*(-1); 
um1 = -V+V*U*uh.*(-1 + erfu.*(uh.^2 +     uw.^2/2) + expu) ;
um2 = -V+V*U*uh.*(-1 +        uh.^2 + 3/2*uw.^2          );
um3 = -V+V*U*uh.*(-1 + erfu.*(uh.^4 + 3*uh.^2.*uw.^2 + 3/4*uw.^4) + expu.*(uh.^2 + 5/2*uw.^3));

pm0 = P + g*B.*bh.*(-1 + erfb.*(bh.^2/2 + bw.^2/4) + expb/2);
pm1 = P + g*B.*bh.*(-1 +        bh.^2/3 + bw.^2/2          );
pm2 = P + g*B.*bh.*(-1 + erfb.*(bh.^4 + 3*bh.^2*bw.^2 + 3/4*bw.^4)/4 + expb.*(bh.^2/4 + 5/8*bw.^2));
pm3 = P + g*B.*bh.*(-1 +        bh.^4/5 + bh.^2*bw.^2 + 3/4*bw.^4                                 );



bm0 = B.*bh.*(1);
bm1 = B.*bh.*(erfb.*(bh.^2 +     bw.^2/2) + expb); 
bm2 = B.*bh.*(       bh.^2 + 3/2*bw.^2          ); 
bm3 = B.*bh.*(erfb.*(bh.^4 +   3*bh.^2.*bw.^2 + 3/4*bw.^4) + expb.*(bh.^2 + 5/2*bw.^2);


bum0 = -B.*(V*bh + U.*bh.*uh./H) + B.*U.*ww.*(hm.*erfm + expm./spi-hd.*erfd- expd./spi)./2;
bum0 = -B.*(V*bh + U.*bh.*uh./H) + B.*U.*ww.*(hm -(1+hd.^2)/spi)./2;
uum0 =  U.^2.*s2.*uw.*expu2./(2.*spi) + U.^2.*uh.*erfu2 - U.^2.*s2.*uw./(2.*spi) - U.^2.*uh.^2./H + H.*V.^2;
uum0 =  U.*U.*(uh - s2.*uw./(2.*spi) - uh.^2./H) + H.*V.*V;
p=a.p(1,f);




pm0 = H.^1/1*p +g.*B.*(2.*bh.^2 + bw.^2).*erfb./4 + bh.*bw.*B.*g.*expb./(2.*spi) - g.*B.*H.*bh;
pm2 = H.^3/3*p +g.*B.*(4.*bh.^4 + 12.*bh.^2.*bw.^2 + 3.*bw.^4).*erfb./48 + g.*B.*bh.*bw.*(2.*bh.^2 + 5.*bw.^2).*expb./(24.*spi) - g.*B.*H.^3.*bh./3;
pm1 = H.^2/2*p +g.*B.*bh.*( -6.*H.^2 +2.*bh.^2  +3.*bw.^2)/12;
pm3 = H.^4/4*p +g.*B.*bh.*(-20.*H.^4 +4.*bh.^4 +15.*bw.^4 +20.*bh.^2.*bw.^2)./80;pm0 := g*B*(2*bh^2 + bw^2)*erf(bh/bw)/4 + g*B*bh*bw*expb/(2*sqrt(Pi)) - H*(B*bh*g - P)


pm0 := g*B*(2*bh^2 + bw^2)*erf(bh/bw)/4 + g*B*bh*bw*expb/(2*sqrt(Pi)) - H*(B*bh*g - P)
pm1 := B*g*bh^3/6 - g*B*(H^2 - bw^2/2)*bh/2 + P*H^2/2;
pm2 := g*B*(4*bh^4 + 12*bh^2*bw^2 + 3*bw^4)*erf(bh/bw)/48 + g*B*bh*bw*(2*bh^2 + 5*bw^2)*expb/(24*spi) - H^3*(B*bh*g - P)/3

pm3 := B*g*bh^5/20 + B*g*bh^3*bw^2/4 - ((H^4 - (3*bw^4)/4)*g*B*bh)/4 + P*H^4/4




figure
subplot(4,3, 1);plot(x,a.um0(f)-um0);ylabel('u0');subplot(4,3, 2);plot(x,a.bm0(f)-bm0);ylabel('b0');subplot(4,3, 3);plot(x,a.pm0(f)-pm0);;ylabel('p0');
subplot(4,3, 4);plot(x,a.um1(f)-um1);ylabel('u1');subplot(4,3, 5);plot(x,a.bm1(f)-bm1);ylabel('b1');subplot(4,3, 6);plot(x,a.pm1(f)-pm1);;ylabel('p1');
subplot(4,3, 7);plot(x,a.um2(f)-um2);ylabel('u2');subplot(4,3, 8);plot(x,a.bm2(f)-bm2);ylabel('b2');subplot(4,3, 9);plot(x,a.pm2(f)-pm2);;ylabel('p2');
subplot(4,3,10);plot(x,a.um3(f)-um3);ylabel('u3');subplot(4,3,11);plot(x,a.bm3(f)-bm3);ylabel('b3');subplot(4,3,12);plot(x,a.pm3(f)-pm3);;ylabel('p3');

keyboard
plot(x,a.bum0(f),x,bum0);
plot(x,a.uum0(f)-uum0);

Exm0 = uum0
plot(x,u1,x,a.u(1,f));


U(X).*H - (H./2 - h(X)).*U(X).^2 - h(X)./F.^2

F=V./sqrt(H.*g.*a.b1);


P1=a.u1.^2./2+a.p1;
dwdz=a.dwdz(end,f)-a.dwdz(1,f);

wwm0=w'*a.ww(:,f);
plot(x,a.pm0(f),x,-g*a.bm1(f)+p + wwm0/2);

a.pm0=w'*(a.p-a.p(1,:));
pmm0b=-g*(H*bm0-bm1)-w'*a.ww;
bb=-g*(W0*a.b(:,f));
plot(x,a.pm0(f),x,bb);


