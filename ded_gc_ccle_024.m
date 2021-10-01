dmat='~/mat-ccle-024';
param=load([dmat '/024-param.mat']);param=param.p;
stats=load([dmat '/024-stats.mat']);stats=stats.s
traj=load([dmat '/024-traj.mat']);
sa=load([dmat '/024-steady.mat']);
ta=load([dmat '/024-time.mat']);
sa=load([dmat '/024-steady.mat']);sa=sa.a;
ta=load([dmat '/024-time.mat']);ta=ta.a;
dfig='~/Dropbox/Jim-Claudia-GC/ofigs';
mv=mean(polyval(dp,[sa.t1 sa.t2]));
n=20;
mt = (t(1:end+1-n)+t(n:end))/2;
dt = (t(1:end+1-n)-t(n:end));
dX = (X(1:end+1-n)-X(n:end));
v=dX./dt;
dp=p(1:2).*[2 1];
t1=t(1);
t2=t(end);
A=param.L*param.H;
W=param.W;
b=ded_read_hdf('b-00000.hdf5');
bb=b.b(:,:,b.x<4.95);
mb=mean(bb(:));
sb=std(bb(:));
brg=[min(bb(:)) max(bb(:))];
disp(sprintf(' b mean %6.3f, range [%6.3f %6.3f], std %6.3f',mb,brg(1),brg(2),sb));
preprint([4 2],8);colour_lines;clf;
hist(bb(:),1000);
xlabel('$\phi$','interpreter','latex');
ylabel('frequency');
set(gca,'ytick',[],'xlim',[0.99 1.01]);
clear('b','bb');
print('-depsc2',[dfig '/024-b-pdf.eps']);



preprint([4 2],8);colour_lines;clf;
plot(t,X,t,polyval(p,t));
axis([0 33.5 5 21]);
xlabel('$t$','interpreter','latex');
ylabel('$X$','interpreter','latex');
set(gca,'position',[0.11 0.21 0.88    0.78]);
a2=axes;
plot(t,X-polyval(p,t),[t1 t2],[0 0]);
set(a2,'position',[0.2 0.6 0.35    0.3],'box','on','xtick',[0 10 20],'ytick',[-0.1 0 0.1]);
set(a2,'xlim',[0 33.5],'ylim',[-0.1 0.1]);
print('-depsc2',[dfig '/024-Xt.eps']);

preprint([4 2],8);colour_lines;clf;
plot(mt,v,mt,polyval(dp,mt),[sa.t1 sa.t2],[mv mv]);
axis([0 33.5 0 0.5]);
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
set(gca,'position',[0.11 0.21 0.88 0.77]);
a2=axes;
plot(mt,v,mt,polyval(dp,mt),[t1 t2],[0 0],[sa.t1 sa.t2],[mv mv]);
text(25,.425,sprintf('v=%6.3f',mv));
set(a2,'position',[0.3 0.35 0.55 0.3],'box','on','xtick',[0 10 20 30],'ytick',[0.40 0.45]);
set(a2,'xlim',[0 33.5],'ylim',[0.39 0.45]);
print('-depsc2',[dfig '/024-vt.eps']);

E0=2.5;
preprint([4 2],8);colour_lines;clf;
ke=(ta.uu+ta.vv+ta.ww)/(2*W*E0);
pe=ta.bz/(W*E0);
h=plot(ta.t,ke,ta.t,pe,ta.t,ke+pe);
xlabel('$t$','interpreter','latex');
ylabel('$E$','interpreter','latex');
legend(h,{'ke','pe','E'},'location','west');
axis([0 33.5 0 1])
print('-depsc2',[dfig '/024-Et.eps']);

dt = (ta.t(1:end+1-n)-ta.t(n:end));
mt = (ta.t(1:end+1-n)+ta.t(n:end))/2;
dke = (ke(1:end+1-n)-ke(n:end))./dt;
dpe = (pe(1:end+1-n)-pe(n:end))./dt;
h=plot(mt,dke,mt,-dpe,mt,-dke-dpe,mt,0*mt);
xlabel('$t$','interpreter','latex');
ylabel('$E$','interpreter','latex');
legend(h,{'dke','-dpe','-dE'},'location','northeast');
axis([0 33.5 -0.02 0.06])
print('-depsc2',[dfig '/024-dEt.eps']);


E0=2.5;
E0=1;
texz=ta.Euw/(W*E0);
tezz=ta.Eww/(W*E0);
teyy=ta.vv/(W*E0);

preprint([4 2],8);colour_lines;clf;
h=plot(ta.t,teyy,ta.t,tezz,ta.t,texy,ta.t,0*ta.t);
xlabel('$t$','interpreter','latex');
ylabel('$E$','interpreter','latex');
legend(h,{'te y','te z','R xz'},'location','northeast');
axis([0 33.5 -0.011 0.11])
print('-depsc2',[dfig '/024-tE.eps']);


