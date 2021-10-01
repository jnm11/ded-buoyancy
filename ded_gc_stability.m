function s=ded_gc_stability(nm)

if iscell(nm)
  for j=1:length(nm)
    ss=ded_gc_stability(nm{j});
    if ~isempty(ss)
      s(j)=ss;
    end
  end
  return;
end

dd=ded_dedalus_data_dir;
fna = sprintf('%s/results/%s/profile-ray.mat',dd,nm);
fnp = sprintf('%s/results/%s/param.mat',dd,nm);
fns = sprintf('%s/results/%s/stability.mat',dd,nm);

if ~isfile(fna)
  s=[];
  return;
end

if file_nt(fns,fna) & isfile(fns) & false
  disp(sprintf('ded_gc_stability: loading %s',fns));
  s=load(fns);
  return;
end

disp(sprintf('ded_gc_stability: making %s',fns));

xmin=1;
hdisp=false;
nc=20;


a=load(fna);a=a.a;
p=load(fnp);p=p.p;
h=ded_gc_find_head(a.x,a.bm0,xmin,hdisp);
a.g=p.g;
x=linspace(xmin,h.Xt,nc);
a=ded_gc_erf_RJ(a);
uw = interp1(a.x,a.uw,x,'linear');
bw = interp1(a.x,a.bw,x,'linear');
uh = interp1(a.x,a.uh,x,'linear');
bh = interp1(a.x,a.bh,x,'linear');
u1 = interp1(a.x,a.u1,x,'linear');
u2 = interp1(a.x,a.u2,x,'linear');
b1 = interp1(a.x,a.b1,x,'linear');
b2 = interp1(a.x,a.b2,x,'linear');
J  = interp1(a.x,a.J, x,'linear');
z1=0;
z2=p.H;
lbc='slip';
ubc='slip';

m.type='fd';
m.test=false;
m.sparse=false;
m.order=9;
m.n=500;

alpha = zeros(nc,1);
c     = zeros(nc,1);
x     = x(:);

mina=0.01;
maxa=1.2;

for j=1:nc
  pr(j)=ur_profiles('erf',uw(j),bw(j),J(j),lbc,ubc,uh(j),bh(j),z1,z2,u1(j),u2(j),b1(j),b2(j),p.g);
  a=linspace(mina,maxa,5)';
  s=TGOS(p.Re,a,pr(j),m);
  for k=1:5
    f=findmax(imag(a.*s));
    f=max(2,min(length(a)-1,f));
    a1=(a(f)+a(f-1))/2;a2=(a(f)+a(f+1))/2;
    s1=TGOS(p.Re,a1,pr(j),m);
    s2=TGOS(p.Re,a2,pr(j),m);
    a=[a(1:f-1);a1;a(f);a2;a(f+1:end)];
    s=[s(1:f-1);s1;s(f);s2;s(f+1:end)];
  end
  f=findmax(imag(a.*s));
  f=max(2,min(length(a)-1,f));
  aa{j}=a;
  cc{j}=s;
  alpha(j)=a(f);
  c(j)=s(f);
end
plot(a,imag(s).*a,'s-');title(nm);drawnow;

save(fns,'h','x','c','alpha','pr','m','aa','cc');

clear('s');
s.h=h;
s.x=x;
s.c=c;
s.alpha=alpha;
s.pr=pr;
s.m=m;
s.aa=aa;
s.cc=cc;





return;

clear('a');


fns=cellstr_ls('~/data/dedalus/results/gc/f6/i/8000/*/4608/*/stability.mat');
fns=cellstr_ls('~/data/dedalus/results/gc/f6/mg/4000/*/2304/*/stability.mat');
dd=[ded_dedalus_data_dir '/results/'];
fns=jfm_mgi_nms(1000)
 
for j=1:length(fns);
  fnmat=[dd fns{j} '/stability.mat'];
  a(j)=load(fnmat);
  typ{j}=ded_gc_f7_g_classification(fns{j});
end

figure(1);clf;
H=[];
for j=1:length(a)
  x=a(j).h.Xt-a(j).x;
  f=find(x>=2*a(j).h.ht);
  g=imag(a(j).alpha.*a(j).c);
  c=real(a(j).c);
  J=[a(j).pr.J];
  R=[a(j).pr.uw]./[a(j).pr.bw];
  subplot(4,1,1);  h    = plot(x(f),g(f));hold('on');ylabel('gr');
  subplot(4,1,2);  h(2) = plot(x(f),c(f)); hold('on');ylabel('c');
  subplot(4,1,3);  h(3) = plot(x(f),J(f)); hold('on');ylabel('J');
  subplot(4,1,4);  h(4) = plot(x(f),R(f)); hold('on');ylabel('R');
  switch(typ{j})
    case('S')
      H(1)=h(1);
      set(h,'color',[1 0 0]);
    case('K')
      H(2)=h(1);
      set(h,'color',[0 0 1]);
     case('H')
      H(3)=h(1);
      set(h,'color',[0 1 0]);
   end
end
legend(H,{'S','K','H'},'location','best');


for k=1:length(a)
  figure(1+k);clf;
  for j=1:length(a(k).aa)
    plot(a(k).aa{j},imag(a(k).cc{j}).*a(k).aa{j},'-s');
    hold('on');
  end
  title(fns{k});
end
