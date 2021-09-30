function ded_process_gc_f7_xxxxh(bn)
% h simulations;

nn={'xxxxa','xxxxb','xxxxc','xxxxd','xxxxe','xxxxf','xxxxg','xxxxh','[0-9]*'};
on={'xxxxa','xxxxb','xxxxc','xxxxd','xxxxe','xxxxf','xxxxg','xxxxh','fixed'};
nn={    'xxxxh',     '[0-9]*',     'xxxxh',    '[0-9]*'};
nx={'2304', '2304',  '4608', '4608'};
on={'xxxxh','fixed','xxxxh','fixed'};


bn  = {     'i',      'i',      'i',     'i'};
bnn = {'*[ig]', '*[ig]', '*[ig]','*[ig]'};
nx  = {  '2304',   '2304',  '4608',   '4608'};
on  = { 'fixed', 'fixed'  ,'xxxxh',  'xxxxh'};
nn  = {'[0-9]*','[0-9]*',  'xxxxh',  'xxxxh'};

% $$$ if nargin<1;
% $$$   bn={'fnsf','g','h','i','j','mg'};
% $$$ end;
showb=true;
ctp.typ={'X','g'};
ctp.tol=[1e-3,1e-4];
q.nbw=2;q.nbh=1;q.nuw=2;q.nuh=1;q.nu1=1;
dd=[ded_dedalus_data_dir '/results/gc/f6'];
cd('~/');
for i=1:length(bn)
  fns=cellstr_ls(sprintf('gc/f6/%s/*/*/%s/%s/time',bnn{i},nx{i},nn{i}));
  if isempty(fns);
    continue;
  end;
  fnmat=[dd '/i-' nx{i} '-' on{i} '.mat'];
  fnc=[dd '/' bn{i} '-' nx{i} '-' on{i} '-classify.mat'];
  fnr=[dd '/' bn{i} '-' nx{i} '-' on{i} '-RJ.mat'];
  if all(file_nt(fnmat,fns)) & isfile(fnmat);
    continue;
  end;
  disp(fnmat);
  fns=cellstr_ls(sprintf('gc/f6/%s/*/*/%s/%s',bn{i},nx{i},nn{i}));
  p=ded_read_param(fns);
  ded_make_param(fns);
  ded_make_avg(fns,'ay',ctp);
  dc;b=ded_gc_extrap_uw(fns,'ray',[5 30],[0.1 35],showb,q,'wh30');
  if isempty(b);
    continue;
  end;
  p=ded_read_param({b.nm});
  save(fnmat','b','p');
  xmin=1;
  a=ded_gc_f7_g_classify(fns,xmin,true);   save(fnc,'a');
  a=ded_gc_f7_g_RJ(fns,xmin,false,'front');save(fnr,'a');
end;

if false;


  if false;

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
  end;

end;


load(fnmat);
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

;
dc;
p2.mfc=[];

figure;clf;[h lh uc]=groupplot(Re,guw,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('uw');groupplot(Re,uw,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,gu1,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('U1');groupplot(Re,U1,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,guh,Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('Uh');groupplot(Re,uh,Pe,p2);


ctp.typ={'X','g'};
ctp.tol=[1e-3,1e-4];

fns={'gc/f6/mg/1000/*/2*/[0-9]*',...
     'gc/f6/mg/2000/*/2*/[0-9]*',...
     'gc/f6/mg/2000/*/4*/[0-9]*',...
     'gc/f6/mg/4000/*/2*/[0-9]*',...
     'gc/f6/mg/4000/*/4*/[0-9]*',...
     'gc/f6/mg/8000/*/4*/[0-9]*',...
     'gc/f6/i/8000/*/4*/[0-9]*'};
for k=1:length(fns)
  disp('');;
  [m s]=ded_stats_avrg(fns{k},ctp);
  disp(sprintf('%30s %5s %5s %7s %7s','name','t2','dt','X','sX'));
  for j=1:length(m)
    disp(sprintf('%30s %6.1f %6.1f %7.4f %7.5f',m(j).nm,m(j).t2,m(j).dt,m(j).X-40,s(j).X'));
  end
end


ctp.typ={'X','g'};
ctp.tol=[1e-3,1e-4];

fns={'gc/f6/mg/4000/*/4*/*h'};
for k=1:length(fns)
  disp('');;
  [m s p]=ded_stats_avrg(fns{k},ctp);
  disp(sprintf('%30s (%6s %6s %6s) (%7s %6s) (%6s %6s)','name','t1','t2','dt','X','sX','mg','sg'));
  for j=1:length(m)
    dX=sqrt((m(j).X-40)^2+s(j).X^2);
    disp(sprintf('%30s (%6.1f %6.1f %6.1f) (%7.4f %6.4f) (%6.4f %6.4f) ',m(j).nm,m(j).t1,m(j).t2,m(j).dt,m(j).X-40,dX,m(j).g,s(j).g));
  end
end
ded_plot_X(fns,-2);

mg=[m.g];
sg=[s.g];
Re=[p.Re];
X=[m.X];
figure;
subplot(2,1,1);plot(Re,mg,'s-',Re,mg+sg,'s-',Re,mg-sg,'s-');
subplot(2,1,2);plot(Re,X,'s-');

