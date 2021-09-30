if false
  cd('~/');fns=cellstr_ls('gc/f7/mg/*/*/2304/xxxx[ab]');
  cd('~/');fns=cellstr_ls('gc/f6/mg/*/*/4608/xxxxi');
   
  ctp.typ={'X','g','wu','U1'};
  ctp.tol=[1e-3 1e-4 1e-3 1e-3];
  ctp.T=0;
  
  [m s]=ded_stats_avrg(fns,ctp);
  p=ded_read_param(fns);
  
  Re=[p.Re];
  Pe=round([p.Re].*[p.Scb]);
  dT=[p.dT];
  wu=[m.wu];
  U1=[m.U1];
  g =[m.g];
  
  f=find(Re<600);
  f1 = @(p,R) 1./max(eps,p(1)+p(2)*R/1000);
  fi =@(p) f1(p,Re(f))-wu(f);
  f3 = @(p,R) (1+erf(p(1)*(R/1000-p(2))))/2;
  fwu = @(p,R) f1(p(1:2),R)+f3(p(3:4),R).*(p(5)-f1(p(1:2),R));
  fi =@(p) fwu(p,Re)-wu;
  pwu = lsqnonlin(fi,[1 1 1 1 0.4 0]);
  
  fg = @(p,R) polyval(p(1:3),R/1000)./polyval([p(4:5) 1],R/1000);
  fi =@(p) fg(p,Re)-g;
  pg = lsqnonlin(fi,[1.4197    1.2041    0.4363    0.6308    0.6806]*100);
  
  %fU1 = @(p,R) p(1)+p(5)*R/1000+p(2)*erf(p(3)*(R/1000-p(4)));
  %fi =@(p) fU1(p,Re)-U1;
  %U1 = lsqnonlin(fi,[0.8 0.4 1 1 0]);
  fU1 = @(p,R) polyval(p(1:3),R/1000)./polyval([p(4:5) 1],R/1000)
  fi =@(p) fU1(p,Re)-U1;
  pU1 = lsqnonlin(fi,[0 0 0.4 0 0]);
  
  
  RR =linspace(min(Re),max(Re),1e3);
  figure(1);clf;subplot(4,1,1);groupplot(Re,wu,            Pe);hold('on');plot(RR,fwu(pwu,RR));ylabel('wu');
  figure(2);clf;subplot(4,1,1);groupplot(Re,wu-fwu(pwu,Re),Pe);                   ylabel('dwu');
  
  figure(1);subplot(4,1,2);groupplot(Re,U1,            Pe);hold('on');plot(RR,fU1(pU1,RR));ylabel('U1');
  figure(2);subplot(4,1,2);groupplot(Re,U1-fU1(pU1,Re),Pe);                   ylabel('dU1');
  
  figure(1);subplot(4,1,3);groupplot(Re,g,          Pe);hold('on'); plot(RR,fg(pg,RR));ylabel('g');
  figure(2);subplot(4,1,3);groupplot(Re,g-fg(pg,Re),Pe);                    ylabel('dg');
  
  figure(1);subplot(4,1,4);groupplot(Re,[s.t1],Pe);ylabel('t1');xlabel('Re');
  figure(2);subplot(4,1,4);groupplot(Re,[s.dt],Pe);ylabel('dt');xlabel('Re');
  
  cU1=fU1(pU1,Re);
  cwu=fwu(pwu,Re);
  cg=fg(pg,Re);
  for j=1:length(fns)
    disp(sprintf('--wu %6.4f --U1 %6.4f --g %6.4f %s/%5.0f',cwu(j),cU1(j),cg(j),p(j).name,cg(j)*1e4));
  end
  
  for j=1:length(fns)
    disp(sprintf('--wu %6.4f --U1 %6.4f --g %6.4f %s/%5.0f',cwu(j),cU1(j),cg(j),p(j).name,cg(j)*1e4));
  end
end




if false
  cd('~/');
  fns=cellstr_ls('gc/f6/mg*/*/*/*/xxxxi');
  ded_plot_stats(fns,{'U1','wu'},{'Re'});
  showb=true;
  ctp.typ={'X','g'};
  ctp.tol=[1e-3,1e-4];
  %ctp.T=100;
  ded_make_avg(fns,'ay',ctp);
  q1.nbw=2;q1.nbh=1;q1.nuw=2;q1.nuh=1;q1.nu1=1;
  dc;
  b=ded_gc_extrap_uw(fns,'ray',[ 10 20],[0.1 35],showb,q1,'uw')
  
  p=ded_read_param({b.nm});
  Re=[p.Re];
  Pe=round([p.Re].*[p.Scb]);
  
  f=find(Re <= 1400 & Pe== 1000);
  ded_gc_fit_Re_hw1(Re(f),[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name});
  
  f=find(Re<= 1000 & Pe== 2000);
  ded_gc_fit_Re_hw1(Re(f),[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name});
  
  f=find(Re<= 1000 & Pe== 4000);
  ff=ded_gc_fit_Re_hw1(Re(f),[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name});
  
  ff.fuw([100 150 200 250 300:100:2000])
end

clear all
showb=true;
ctp.typ={'X'};  
ctp.tol=[1e-3];
cd('~/');
fn1=cellstr_ls('gc/f6/mg/*/*/*/[0-9]*');
fn2=cellstr_ls('gc/f6/i/*/*/*/[0-9]*');
fns=cat(1,fn1,fn2);
fns=cellstr_ls('gc/f6/mg/*/*/*/*h');
p=ded_read_param(fns);
ded_make_param(fns);
ded_make_avg(fns,'ay',ctp);
q.nbw=2;q.nbh=1;q.nuw=2;q.nuh=1;q.nu1=1;
dc;b=ded_gc_extrap_uw(fns,'ray',[5 30],[0.1 35],showb,q,'wh30');
p=ded_read_param({b.nm});
Re=[p.Re];
Pe=round([p.Re].*[p.Scb]);
uw  = [p.wu];
U1  = [p.U1];
uh  = [p.hu];
Nx  = [p.Nx];
g   = [p.g];
guw = [b.guw];
gu1 = [b.gu1];
guh = [b.guh];


dc;
p2.mfc=[];

figure;clf;[h lh uc]=groupplot(Re,guw,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('uw');groupplot(Re,uw,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,gu1,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('U1');groupplot(Re,U1,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,guh,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('Uh');groupplot(Re,uh,Pe,p2);


if false
  
  %f=find(Re<1000);a=ded_gc_fit_Re_hw1(Re(f),guh(f),guw(f),gu1(f),g(f),{fns{f}},'mid',Pe(f));
  
  Re=unique(Re);
  Pe=unique(Pe);
  [Re Pe]=ndgrid(Re,Pe);
  Re=Re(:)';Pe=Pe(:)';
  disp(sprintf('--pfn  -rfn --reset --wu %6.4f --U1 %6.4f %04.0f/%04.0f/2304/xxxxh\n',[a.fuw(Re,Pe);a.fu1(Re,Pe);Pe;Re]));
  
  
  
% $$$ disp(sprintf('--pfn  -rfn --reset --wu %6.4f --U1 %6.4f %04.0f/%04.0f/%04.0f/xxxxh\n',[a.fuw(Re,Pe);a.fu1(Re,Pe);Pe;Re;Nx]));
% $$$ 
% $$$ 
% $$$ f=find(Pe==4000 & Nx==2304);a=ded_gc_fit_Re_hw1(Re(f),guh(f),guw(f),gu1(f),g(f),{fns{f}},'low',Pe(f));
% $$$ 
% $$$ 
% $$$ f=find(Pe==4000 & Nx==2304 & Re<=600);a=ded_gc_fit_Re_hw1(Re(f),guh(f),guw(f),gu1(f),g(f),{fns{f}},{'poly2','poly2','rat11','poly2'},Pe(f));
end


























