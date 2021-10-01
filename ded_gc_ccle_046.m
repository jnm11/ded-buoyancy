
nm='ccle/046';
nm='emle/019';

dmat=['~/Dropbox/Jim-Claudia-GC/mat/' nm];

param=load([dmat '/param.mat']);param=param.p;
stats=load([dmat '/stats.mat']);stats=stats.s;
traj=load([dmat  '/traj.mat']);traj=traj.a;
ta=load([  dmat  '/time.mat']);ta=ta.a;
%sa=load([  dmat  '-steady.mat']);sa=sa.a;
sa.t1=10;
sa.t2=33.5;
p=traj.p;
dp=traj.dp;

nmf=nm;
nmf(nmf=='/')='-';
dfig='~/Dropbox/Jim-Claudia-GC/ofigs';
mv=mean(traj.fu([sa.t1 sa.t2]));
t=stats.t;
X=stats.X;
n=20;
mt = (t(1:end+1-n)+t(n:end))/2;
dt = (t(1:end+1-n)-t(n:end));
dX = (X(1:end+1-n)-X(n:end));
v=dX./dt;
t1=t(1);
t2=t(end);
A=param.L*param.H;
W=param.W;
fnb=[ded_dedalus_data_dir '/' nm '/b/b-00000.hdf5'];
if isfile(fnb)
  b=ded_read_hdf(fnb);
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
  print('-depsc2',[dfig '/' nmf '-b-pdf.eps']);
end



preprint([4 2],8);colour_lines;clf;
plot(t,X,t,traj.fx(t));
axis([0 33.5 5 21]);
xlabel('$t$','interpreter','latex');
ylabel('$X$','interpreter','latex');
set(gca,'position',[0.11 0.21 0.88    0.78]);
a2=axes;
plot(t,X-traj.fx(t),[t1 t2],[0 0]);
set(a2,'position',[0.2 0.6 0.35    0.3],'box','on','xtick',[0 10 20],'ytick',[-0.1 0 0.1]);
set(a2,'xlim',[0 33.5],'ylim',[-0.1 0.1]);
print('-depsc2',[dfig '/' nmf '-Xt.eps']);

preprint([4 2],8);colour_lines;clf;
plot(mt,v,mt,traj.fu(mt),[sa.t1 sa.t2],[mv mv]);
axis([0 33.5 0 0.5]);
xlabel('$t$','interpreter','latex');
ylabel('$v$','interpreter','latex');
set(gca,'position',[0.11 0.21 0.88 0.77]);
a2=axes;
plot(mt,v,mt,traj.fu(mt),[t1 t2],[0 0],[sa.t1 sa.t2],[mv mv]);
text(25,.425,sprintf('v=%6.3f',mv));
set(a2,'position',[0.3 0.35 0.55 0.3],'box','on','xtick',[0 10 20 30],'ytick',[0.40 0.45]);
set(a2,'xlim',[0 33.5],'ylim',[0.39 0.45]);
print('-depsc2',[dfig '/' nmf '-vt.eps']);

E0=2.5;
preprint([4 2],8);colour_lines;clf;
ke=(ta.uu+ta.vv+ta.ww)/(2*W*E0);
pe=ta.bz/(W*E0);
h=plot(ta.t,ke,ta.t,pe,ta.t,ke+pe);
xlabel('$t$','interpreter','latex');
ylabel('$E$','interpreter','latex');
legend(h,{'ke','pe','E'},'location','west');
axis([0 33.5 0 1])
print('-depsc2',[dfig '/' nmf '-Et.eps']);

dt = (ta.t(1:end+1-n)-ta.t(n:end));
mt = (ta.t(1:end+1-n)+ta.t(n:end))/2;
dke = (ke(1:end+1-n)-ke(n:end))./dt;
dpe = (pe(1:end+1-n)-pe(n:end))./dt;
h=plot(mt,dke,mt,-dpe,mt,-dke-dpe,mt,0*mt);
xlabel('$t$','interpreter','latex');
ylabel('$E$','interpreter','latex');
legend(h,{'dke','-dpe','-dE'},'location','northeast');
axis([0 33.5 -0.02 0.06])
print('-depsc2',[dfig '/' nmf '-dEt.eps']);


E0=2.5;
E0=1;
%texz=ta.Euw/(W*E0);
%tezz=ta.Eww/(W*E0);
texx=ta.uu/(W*E0);
texy=ta.uv/(W*E0);
texz=ta.uw/(W*E0);
teyy=ta.vv/(W*E0);
teyz=ta.vw/(W*E0);
tezz=ta.ww/(W*E0);

preprint([4 2],8);colour_lines;clf;
h=plot(ta.t,teyy,ta.t,tezz,ta.t,texy,ta.t,0*ta.t);
xlabel('$t$','interpreter','latex');
ylabel('$E$','interpreter','latex');
legend(h,{'te y','te z','R xz'},'location','northeast');
axis([0 33.5 -0.011 0.11])
print('-depsc2',[dfig '/' nmf '-tE.eps']);







LL=[-20 5];
trg=[16 33.5];
typ='y';
a=ded_mavrg2(['gc/' nm],trg,traj.fx,traj.fu,typ,LL);
save([dfig '/' nmf '-steady.mat'],'a');

b=ded_zgrid(a);
contour(b.x,b.z,b.b,[0.01 0.02 0.05 0.1 0.5 0.9 0.95 0.98 0.99]);
[w,w2]=ichebintw(length(a.z),param.H)
B1=w*a.b;
B2=(2*(w.*a.z')*a.b)./max(1e-3,B1);
B3=B1.^2./max(1e-3,w*a.b.^2);
plot(a.x,B1,a.x,B2,a.x,B3);


fns=ded_get_fn(nm,'y');
for k=1:length(fns)
  disp(sprintf('%u/%u',k,length(fns)));
  a=ded_read_hdf(fns{k});
  dx=a.x(2)-a.x(1);
  b.u(k)=dx*sum(w*a.u);
  b.w(k)=dx*sum(w*a.w);
  b.v(k)=dx*sum(w*a.v);

  b.uu(k)=dx*sum(w*a.uu);
  b.ww(k)=dx*sum(w*a.ww);
  b.vv(k)=dx*sum(w*a.vv);
  b.uv(k)=dx*sum(w*a.uv);
  b.uw(k)=dx*sum(w*a.uw);
  b.vw(k)=dx*sum(w*a.vw);
  
  b.Auu(k)=dx*sum(w*(a.u.*a.u));
  b.Avv(k)=dx*sum(w*(a.v.*a.v));
  b.Aww(k)=dx*sum(w*(a.w.*a.w));
  b.Auv(k)=dx*sum(w*(a.u.*a.v));
  b.Auw(k)=dx*sum(w*(a.u.*a.w));
  b.Avw(k)=dx*sum(w*(a.v.*a.w));
 
  b.Ruu(k)=dx*sum(w*(a.uu-a.u.*a.u).^2);
  b.Rvv(k)=dx*sum(w*(a.vv-a.v.*a.v).^2);
  b.Rww(k)=dx*sum(w*(a.ww-a.w.*a.w).^2);
  b.Ruv(k)=dx*sum(w*(a.uv-a.u.*a.v).^2);
  b.Ruw(k)=dx*sum(w*(a.uw-a.u.*a.w).^2);
  b.Rvw(k)=dx*sum(w*(a.vw-a.v.*a.w).^2);

  b.Euu(k)=dx*sum(w*(a.uu-a.u.*a.u));
  b.Evv(k)=dx*sum(w*(a.vv-a.v.*a.v));
  b.Eww(k)=dx*sum(w*(a.ww-a.w.*a.w));
  b.Euv(k)=dx*sum(w*(a.uv-a.u.*a.v));
  b.Euw(k)=dx*sum(w*(a.uw-a.u.*a.w));
  b.Evw(k)=dx*sum(w*(a.vw-a.v.*a.w));
 
  b.buu(k)=dx*sum(w*(a.b.*a.uu));
  b.bww(k)=dx*sum(w*(a.b.*a.ww));
  b.bvv(k)=dx*sum(w*(a.b.*a.vv));
  b.buu(k)=dx*sum(w*(a.b.*a.u));
  b.bww(k)=dx*sum(w*(a.b.*a.w));
  b.bvv(k)=dx*sum(w*(a.b.*a.v));
  b.buv(k)=dx*sum(w*(a.b.*a.uv));
  b.buw(k)=dx*sum(w*(a.b.*a.uw));
  b.bvw(k)=dx*sum(w*(a.b.*a.vw));

  b.b(k)=dx*sum(w*a.b);
  b.bz(k)=dx*sum((a.z'.*w)*a.b);
  b.t(k)=a.t1;
  [dfig fn]=fileparts(fns{k});
  b.nm{k}=fn;
end
a=b;
save([dfig '/' nmf '-time.mat'],'a');
