function ded_gc_ccle_plot(nm,tm)
if nargin==0
  dnm=cellstr_ls('~/Dropbox/Jim-Claudia-GC/mat/*le*/*');
  for j=1:length(dnm)
    [dd nm]=fileparts(dnm{j});
    nm=nm(5:end);
    switch(nm)
      case {'ccle/017','ccle/024','ccle/025'}
        tm='-30';
      case {'ccle/020','ccle/021','ccle/022'}
        tm='-18';
      case {'ccle/023'}
        tm='-15';
      case {'ccle/046','emle/001','emle/002','emle/003','emle/019'}
        tm='-24';
    end
    ded_gc_ccle_plot(nm,tm);
  end
  return;
end
%ded_gc_ccle_plot('ccle/046');
%ded_gc_ccle_plot('emle/019);
dc;
dmat=['~/Dropbox/Jim-Claudia-GC/mat/' nm];
dfig='~/Dropbox/Jim-Claudia-GC/ofigs';


fnsteady = [dmat '/' 'steady.mat'];
fntime   = [dmat '/'   'time.mat'];


[traj param stats]=ded_gc_traj(['gc/' nm],[],[5 33.5],1);
if isfield(traj,'t1')
  t1=traj.t1;
  t2=traj.t2;
else
  t1=traj.t(1);
  t2=traj.t(end);
end  
if isfile(fnsteady)
  sa=load(fnsteady);
  sa=sa.a;
  t1=sa.t1;
  t2=sa.t2;
end

if isfile(fntime)
  ta=load(fntime);ta=ta.a;
else
  ta=[];
end

fnb=[dmat '/b-00000.hdf5'];
if isfile(fnb)
  b=ded_read_hdf(fnb);
  bb=b.b(:,:,b.x<4.95);
  mb=mean(bb(:));
  sb=std(bb(:));
  brg=[min(bb(:)) max(bb(:))];
  disp(sprintf(' b mean %6.3f, range [%6.3f %6.3f], std %6.3f',mb,brg(1),brg(2),sb));
  figure;preprint([4 2],8);colour_lines;clf;
  hist(bb(:),1000);
  xlabel('$\phi$','interpreter','latex');
  ylabel('frequency');
  set(gca,'ytick',[],'xlim',[0.99 1.01]);
  clear('b','bb');
  print('-depsc2',[dfig '/' nm '-b-pdf.eps']);
end

t=traj.t;
mt=traj.mt;
X=traj.X;
V=traj.V;
Xf=traj.fx(traj.t);
Vf=traj.fu(mt);
mv=(traj.fu(t1)+traj.fu(t2))/2;
maxt=max(t(X<param.L-1));
minV=min(V(mt>5&mt<maxt));
maxV=max(V(mt>5&mt<maxt));
if isempty(minV); minV=min(V(:)); end;
if isempty(maxV); maxV=max(V(:)); end;


figure;preprint([4 2],8);colour_lines;clf;
plot(t,X,t,Xf);
axis([0 t(end) min(X) max(X)]);
xlabel('$t$','interpreter','latex');
ylabel('$X$','interpreter','latex');
set(gca,'position',[0.11 0.21 0.88    0.78]);
a2=axes;
plot(t,X-Xf,[t1 t2],[0 0]);
set(a2,'position',[0.21 0.64 0.35    0.3],'box','on','xtick',[0 10 20],'ytick',[-0.05 0 0.05]);
set(a2,'xlim',[t(1) t(end)],'ylim',[-0.05 0.05]);
print('-depsc2',[dfig '/' nm '-Xt.eps']);

if length(V)>1
  figure;preprint([4 2],8);colour_lines;clf;
  plot(mt,V,mt,Vf,[t1 t2],[mv mv]);
  axis([0 t(end) 0 0.5]);
  xlabel('$t$','interpreter','latex');
  ylabel('$v$','interpreter','latex');
  set(gca,'position',[0.11 0.21 0.88 0.77]);
  a2=axes;
  plot(mt,V,mt,Vf,[t1 t2],[0 0],[t1 t2],[mv mv]);
  text(25,.425,sprintf('v=%6.3f',mv));
  set(a2,'position',[0.3 0.35 0.55 0.3],'box','on','xtick',[0:10:t(end)],'ytick',[0:0.02:.5]);
  set(a2,'xlim',[0 t(end)],'ylim',[minV maxV]);
  print('-depsc2',[dfig '/' nm '-vt.eps']);
end


if ~isempty(ta)
  E0=ta.bz(1);
  figure;preprint([4 2],8);colour_lines;clf;
  A=param.L*param.H;
  ke=(ta.uu+ta.vv+ta.ww)/2/A;
  pe=ta.bz/A;
  h=plot(ta.t,ke,ta.t,pe,ta.t,ke+pe);
  xlabel('$t$','interpreter','latex');
  ylabel('$E$','interpreter','latex');
  legend(h,{'ke','pe','E'},'location','best');
  axis([0 t(end) 0 E0/A])
  print('-depsc2',[dfig '/' nm '-Et.eps']);
  
  
  figure;preprint([4 2],8);colour_lines;clf;
  n=round(0.2/median(diff(unique(ta.t))));
  dt = (ta.t(1:end+1-n)-ta.t(n:end));
  mt = (ta.t(1:end+1-n)+ta.t(n:end))/2;
  dke = (ke(1:end+1-n)-ke(n:end))./dt;
  dpe = -(pe(1:end+1-n)-pe(n:end))./dt;
  f=find(mt>1 & mt<maxt);
  if isempty(f);f=1:length(dke);end
  maxDE=max([max(dke(f)) max(dpe(f))]); 
  minDE=min([min(dke(f)) min(dpe(f))]); 
  h=plot(mt,dke,mt,dpe,mt,-dke+dpe,mt,0*mt);h(end)=[];
  xlabel('$t$','interpreter','latex');
  ylabel('$E$','interpreter','latex');
  legend(h,{'dke','-dpe','-dE'},'location','best');
  axis([0 t(end) minDE maxDE])
  print('-depsc2',[dfig '/' nm '-dEt.eps']);
  
  texx=(ta.uu-ta.Auu)/A;
  texz=(ta.uw-ta.Auw)/A;
  tezz=(ta.ww-ta.Aww)/A;
  teyy=(ta.vv)/A;
  
  maxE=max(texx(ta.t>1));
  figure;preprint([4 2],8);colour_lines;clf;
  h=plot(ta.t,texx,ta.t,teyy,ta.t,tezz,ta.t,texz,ta.t,0*ta.t);h(end)=[];
  xlabel('$t$','interpreter','latex');
  ylabel('$E$','interpreter','latex');
  legend(h,{'Rxx','Ryy','Rzz','Rxz'},'location','best');
  axis([0 t(end) 0 3e-3])
  print('-depsc2',[dfig '/' nm '-TE.eps']);
  
  T=20;
  f=find(ta.t<T);
  maxE=max(texx(ta.t>1 & ta.t<T));
  maxE=0.1;
  figure;preprint([4 2],8);colour_lines;clf;
  h=plot(ta.t(f),texx(f),ta.t(f),teyy(f),ta.t(f),tezz(f),ta.t(f),texz(f),ta.t(f),0*ta.t(f));h(end)=[];
  xlabel('$t$','interpreter','latex');
  ylabel('$E$','interpreter','latex');
  legend(h,{'Rxx','Ryy','Rzz','Rxz'},'location','best');
  axis([0 T 0 3e-3])
  print('-depsc2',[dfig '/' nm '-TE20.eps']);
end


