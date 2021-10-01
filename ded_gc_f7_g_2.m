
wrg=[0 0.7];
hrg=[0 0.9];
cd('~/data/dedalus/results/gc/f6/g');
fns=cellstr_ls('*/*/*/*/profile.mat',[],'dir');
for j=1:length(fns)
  nm=fns{j};
  pr=load([nm '/profile.mat']);pr=pr.a;
  p=load([nm '/param.mat']);p=p.p;
  c=load([nm '/coord.mat']);c=c.c;
  X=pr.X;%X=max(c.x(pr.hb>0.01));
  
  Re(j)=p.Re(j);
  Sc(j)=p.Scb(j)*p.Re(j);
  
  X=20;
  u0(j)=interp1(pr.x,pr.u0,X,'linear');
  u1(j)=interp1(pr.x,pr.u1,X,'linear');
  b0(j)=interp1(pr.x,pr.b0,X,'linear');
  b1(j)=interp1(pr.x,pr.b1,X,'linear');
  hu(j)=interp1(pr.x,pr.hu,X,'linear');
  hb(j)=interp1(pr.x,pr.hb,X,'linear');
  wu(j)=interp1(pr.x,pr.wu,X,'linear');
  wb(j)=interp1(pr.x,pr.wb,X,'linear');
  g(j)=p.g;

  
  dudz(j)=
  
  
  findmin(c.x-20);
  
  
  Re(j)=p.Re.
  
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

