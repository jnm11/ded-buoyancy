cd('~/');fns=cellstr_ls('gc/f6/mg/1000/*/4608/xxxxi');
ded_make_param(fns);

ctp.typ={'X', 'g', 'wu','u0','u1'};
ctp.tol=[1e-3 1e-4 1e-3 1e-3 1e-3];
ctp.T=0;

[m s p]=ded_stats_avrg(fns,ctp);
[Re f]=sort([p.Re]);
m=m(f);s=s(f);p=p(f);


Re = [p.Re];
Pe = round([p.Re].*[p.Scb]);
wu = [m.wu];
hu = [m.hu];
u0 = [m.u0];
u1 = [m.u1];
g  = [m.g];
db = [m.db];
bu = [m.bu];
subplot(5,1,1);plot(Re,wu,'s-');ylabel('wu');
subplot(5,1,2);plot(Re,u0,'s-');ylabel('u0');
subplot(5,1,3);plot(Re, g,'s-');ylabel( 'g');
subplot(5,1,4);plot(Re,db,'s-');ylabel('db');
subplot(5,1,5);plot(Re,log10(abs(bu)),'s-');ylabel('bu');




ded_make_avg(fns,'ay',ctp);
showb=true;
q.nbw=2;q.nbh=1;q.nuw=2;q.nuh=1;q.nu1=1;
dc;b=ded_gc_extrap_uw(fns,'ray',[5 30],[0.1 35],showb,q,'wh30');
guw = [b.guw];
gu1 = [b.gu1];
guh = [b.guh];


dc;
p2.mfc=[];

figure;clf;[h lh uc]=groupplot(Re,guw,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('uw');groupplot(Re,wu,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,gu1,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('U1');groupplot(Re,u0,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,guh,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('Uh');groupplot(Re,hu,Pe,p2);


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


























