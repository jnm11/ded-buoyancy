ded_make_param(cellstr_ls('gc/f6/[hg]/*/*/*/*/status',[],'dir'));
ded_make_avg(cellstr_ls('gc/f6/g/*/*/*/[23]*/status',[],'dir'),'ay');
ded_gc_fit_erf_profiles;

%dc;figure;[c f pp]=ded_gc_fit_g_Re('gc/f6/g/8000/*/*/*',8,'rat34');

%rsync -vap hamilton:data/dedalus/results ~/data/dedalus --exclude ay.mat

if false
  cd('~/');
  fns=cellstr_ls('gc/f6/g/*/*/*/[23]*/status',[],'dir');
  qqq.trg=[-inf inf];
  for j=1:length(fns);
    nm=fns{j};
    fnp=[ded_dedalus_data_dir '/results/' nm '/param.mat'];
    fnc=[ded_dedalus_data_dir '/results/' nm '/coord.mat'];
    [dd fn]=fileparts(fnp);
    if ~isdir(dd); mkdir(dd); end
    p=ded_read_param(nm);
    c=ded_coord(nm);
    save(fnp,'p');
    save(fnc,'c');
  end
end
  
dd='~/Dropbox/GFD-2019-Gravity-Currents/JFM-Paper/ofigs/g';
fsz=[4 2];
fnt=10;

% $$$ nm='1000/0500/2304/27800';
% $$$ pr=load([nm '/profile.mat']);pr=pr.a;
% $$$ p=load([nm '/param.mat']);p=p.p;
% $$$ c=load([nm '/coord.mat']);c=c.c;

wrg=[0 0.7];
hrg=[0 0.9];
cd('~/data/dedalus/results/gc/f6/g');
fns=cellstr_ls('*/*/*/*/profile.mat',[],'dir');
xmin=1;
for j=1:length(fns)
  nm=fns{j};
  pr=load([nm '/profile.mat']);pr=pr.a;
  p=load([nm '/param.mat']);p=p.p;
  c=load([nm '/coord.mat']);c=c.c;
  fr=ded_gc_find_head(pr.x,pr.b,xmin);
  %fru=ded_gc_find_head(pr.x,max(0,pr.uu-p.H),xmin);
  
  f=find(c.x<=fr.Xt & c.x>xmin);
  x=(fr.Xt-c.x)/fr.ht;
  R=pr.wu(:)./pr.wb(:);
  du=pr.u0(:)-pr.u1(:);
  db=pr.b0(:)-pr.b1(:);
  dudz=du./pr.wu/sqrt(pi);
  dbdz=db./pr.wb/sqrt(pi);
  Re=p.Re*du.*pr.hb(:);
  Rib=p.g*db.*pr.hu(:)./du.^2;
  J=dbdz./dudz.^2*p.g;
  Fr=sqrt(max(0,(pr.uu(:)-p.H)./(p.g*pr.bz(:))));
  aa.Fbw=pr.Fbw(f);
  aa.Tbw=pr.Tbw(f);
  aa.x=x(f);
  aa.x=x(f);
  aa.R=R(f);
  aa.Re=Re(f);
  aa.Rib=Rib(f);
  aa.J=J(f);
  aa.Fr=J(f);
  aa.nRe=p.Re; 
  aa.nPe=round(p.Scb*p.Re);
  aa.g = p.g;
  aa.nm=nm;
  a(j)=aa;
end

save([gfd_data_dir '/dns.mat'],'a');


figure(1);clf;hold('on');
for j=1:length(a);plot(a(j).R,a(j).J);end
set(gca,'xscale','log','yscale','log');
xlabel('R');ylabel('J');

figure(2);clf;hold('on');
for j=1:length(a);plot(a(j).R,a(j).Re);end
set(gca,'xscale','log','yscale','log');
xlabel('R');ylabel('Re');

figure(3);clf;hold('on');
for j=1:length(a);plot(a(j).R,a(j).Rib);end
set(gca,'xscale','log','yscale','log');
xlabel('R');ylabel('Rib');

figure(4);clf;hold('on');
for j=1:length(a);plot(a(j).x,a(j).J);end
set(gca,'xscale','linear','yscale','log');
xlabel('x');ylabel('J');


X=[0 25 50 75 100];
col=[[1 0 0];[0 0 1];[1 1 0];[1 0 1];[0 0 1];[0 0 0]];
X=linspace(0,100,200);
col=jet(length(X));
mk='s<>^vpdo*';
k=50;plot(a(k).x,a(k).R)
figure(1);gfd_plot_along_x(a,X,'R','J',col,mk);

figure(2);gfd_plot_along_x(a,X,'R','Re',col,mk);

figure(3);gfd_plot_along_x(a,X,'R','Rib',col,mk);

figure(4);gfd_plot_along_x(a,X,'J','Fbw',col,mk);set(gca,'yscale','linear');


figure(5);gfd_plot_along_x(a,X,'R','Fbw',col,mk);set(gca,'yscale','linear');



figure(2);gfd_plot_along_x(a,X,'Re','Fr',col,mk);set(gca,'yscale','linear');



f=find([a.nPe]==8000);% & [a.nRe]<=1000);
                      %gfd_plot_along_x(a(f),X,'R','J',col,mk);

clf;
for k=1:length(f)
  j=f(k);
  subplot(3,1,1);  plot(a(j).x,a(j).J);hold('on');ylabel('J');
  subplot(3,1,2);  plot(a(j).x,a(j).R);hold('on');ylabel('R');
  subplot(3,1,3);  plot(a(j).R,a(j).J);hold('on'); xlabel('R');ylabel('J');
end
