nm='gc/f6/f/4032/22900';
fd='~/Dropbox/GFD-2019-Gravity-Currents/JFM-Paper/ofigs/4032/22900';
fsz=[4 2];
fnt=10;

%q.trg=[-inf inf];
%qq.sz=[2 4];
%qq.fnt=10;
%a=ded_gc_fit_erf_profiles(nm,q,'ay');

%
%ay=ded_read_javrg(nm,'ay');
%save('~/gc/f6/f/4032/22900/ay.mat','ay');
%rsync -vap hamilton:gc/f6/f/4032/22900/*.hdf5 ~/gc/f6/f/4032/22900
%rsync -vap hamilton:gc/f6/f/4032/22900/ay.mat ~/gc/f6/f/4032/22900
%rsync -vap hamilton:gc/f6/f/4032/22900/*.mat ~/gc/f6/f/4032/22900

load('~/gc/f6/f/4032/22900/ay.mat');
pr=load('~/gc/f6/f/4032/22900/profile.mat');pr=pr.a;
p=ded_read_param(nm);
c=ded_coord(nm);
w=ichebintw(c.Nz);
[w0 w1]=ichebendsw(c.Nz);
dw=w1-w0;

%plot(c.x,w*ay.b);

xrg=[7 42];


rg=find(c.x>xrg(1) & c.x<xrg(2));
X=max(c.x(w*ay.b>0.01));
x=c.x(rg)-X;
xrg=xrg-X;
P=w*(ay.uu+ay.p);P=P(rg);P=P-P(end);
ub=w*ay.bu;ub=ub(rg);
b=w*ay.b;b=b(rg);

dc;

figure;clf;preprint(fsz,fnt);colour_lines;
plot(x,ub); % This shows lack of mass balance
xlabel('$x$','interpreter','latex');
ylabel('$ub$','interpreter','latex');
line(x,0*x,'color',0.7*[1 1 1]);
axis('tight');
print('-depsc2',[fd '/ubx.eps']);

figure;clf;preprint(fsz,fnt);colour_lines;
plot(x,P); % This show momenyum balance
xlabel('$x$','interpreter','latex');
ylabel('$M$','interpreter','latex');
line(x,0*x,'color',0.7*[1 1 1]);
axis('tight');
print('-depsc2',[fd '/Mx.eps']);



figure;clf;preprint(fsz,fnt);colour_lines;
Kxx=w*(ay.uu-ay.u.*ay.u);Kxx=Kxx(rg);
Kzz=w*(ay.ww-ay.w.*ay.w);Kzz=Kzz(rg);
Kxz=w*(ay.uw-ay.w.*ay.u);Kxz=Kxz(rg);
h=plot(x,Kxx,x,Kxz,x,Kzz,'g');
legend(h,{'Kxx','Kxz','Kzz'},'location','best');
xlabel('$x$','interpreter','latex');
ylabel('$K$','interpreter','latex');
axis('tight');
print('-depsc2',[fd '/Kx.eps']);

figure;clf;preprint(fsz,fnt);colour_lines;
Uxx=w*((p.V+ay.u).*(p.V+ay.u));Uxx=Uxx(rg);
Uzz=w*(ay.w.*ay.w);Uzz=Uzz(rg);
Uxz=w*(ay.w.*ay.u);Uxz=Uxz(rg);
h=plot(x,Uxx,x,Uxz,x,Uzz,'g');
xlabel('$x$','interpreter','latex');
ylabel('$R$','interpreter','latex');
axis('tight');
legend(h,{'Rxx','Rxz','Rzz'},'location','west');
print('-depsc2',[fd '/Rx.eps']);



figure;clf;preprint(fsz,fnt);colour_lines;
h=plot(x,pr.hb(rg),x,pr.hu(rg),x,b);
axis([xrg 0 inf]);
xlabel('$x$','interpreter','latex');
ylabel('$h$','interpreter','latex');
legend(h,{'$h_b$','$h_u$','$m$'},'location','southwest','interpreter','latex');
print('-depsc2',[fd '/hx.eps']);


wu=pr.wu(rg);wu(find(wu==min(wu)):end)=NaN;
wb=pr.wb(rg);wb(find(wb==min(wb)):end)=NaN;
figure;clf;preprint(fsz,fnt);colour_lines;
h=plot(x,wb,x,wu);
axis([xrg 0 0.7]);
xlabel('$x$','interpreter','latex');
ylabel('$w$','interpreter','latex');
legend(h,{'$w_b$','$w_u$'},'location','southwest','interpreter','latex');
print('-depsc2',[fd '/wx.eps']);


p1=w1*ay.p+(w1*ay.u).^2/2;p1=p1-p1(end);p1=p1(rg);
p0=w0*ay.p+(w0*ay.u).^2/2;p0=p0-p0(end);p0=p0(rg);
dp=dw*ay.p+w*ay.b*p.g;dp=dp(rg);
figure;clf;preprint(fsz,fnt);colour_lines;
h=plot(x,p0,x,p1,x,dp);
axis([xrg -0.2 0.2]);
legend(h,{'$P_0$','$P_1$','$p_1-p_0+g\int\rho\,dz$'},'location','southwest','interpreter','latex');

