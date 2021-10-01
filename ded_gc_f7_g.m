
dc;figure;[c f pp]=ded_gc_fit_g_Re('gc/f6/g/8000/*/*/*',8,'rat34');

%rsync -vap tttuck:data/dedalus/results ~/data/dedalus 

if false
  cd('~/');
  fns=cellstr_ls('gc/f6/*/*/*/*/*/status',[],'dir');
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

nm='1000/0500/2304/27800';
pr=load([nm '/profile.mat']);pr=pr.a;
p=load([nm '/param.mat']);p=p.p;
c=load([nm '/coord.mat']);c=c.c;

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
  
  x=fr.Xh-c.x;
  
  X=pr.X;%X=max(c.x(pr.hb>0.01));
  if isempty(c);disp(nm);continue;end
  nm(nm=='/')='-';
  fd=[dd '/' nm];
  if ~isdir(fd);mkdir(fd);end;
  w=ichebintw(c.Nz);
  [w0 w1]=ichebendsw(c.Nz);
  dw=w1-w0;
  rg=find(pr.hu>0&pr.hu<2&pr.hb>0&pr.hb<2&pr.wu>0&pr.wu<2&pr.wb>0&pr.wb<2);
  rg=find(c.x>1&c.x<X&pr.wu>0&pr.wb>0&pr.wu<1&pr.wb<1);
  x=c.x(rg)-X;
  xrg=x([1 end])'+[0 1];
  clf;preprint(fsz,fnt);colour_lines;
  h=plot(x,pr.hb(rg),x,pr.hu(rg));
  axis([xrg hrg]);
  xlabel('$x$','interpreter','latex');
  ylabel('$h$','interpreter','latex');
  legend(h,{'$h_b$','$h_u$'},'location','southwest','interpreter','latex');
  print('-depsc2',[fd '/hx.eps']);
    
  wu=pr.wu(rg);%wu(find(wu==min(wu)):end)=NaN;
  wb=pr.wb(rg);%wb(find(wb==min(wb)):end)=NaN;
  wu(max(find(diff(wu)<0)):end)=NaN;
  wb(max(find(diff(wb)<0)):end)=NaN;

  clf;preprint(fsz,fnt);colour_lines;
  h=plot(x,wb,x,wu);
  axis([xrg wrg]);
  xlabel('$x$','interpreter','latex');
  ylabel('$w$','interpreter','latex');
  legend(h,{'$w_b$','$w_u$'},'location','southwest','interpreter','latex');
  print('-depsc2',[fd '/wx.eps']);
end

