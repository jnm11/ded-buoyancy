% 2d Steady State DNS
if false
showb=false;
  
if false
  cd(ded_dedalus_data_dir);
  fns=cellstr_ls('gc/f6/g/*/*/*/[23]*/param.h5',[],'dir');
  %ded_make_param(fns);
  ctp.typ={'X'};
  ctp.tol=[1e-1];
  ded_make_avg(fns,'ay',ctp);
  ded_gc_fit_erf_profiles(fns,[],'ray');
  
  
  
  fns=cellstr_ls('gc/f6/i/*/*/*/*/param.h5',[],'dir');
  ded_make_param(fns);
  ctp.typ={'X','g'};
  ctp.tol=[1e-3,1e-4];
  ded_make_avg(fns,'ay',ctp);
  

  q.nbw=1;q.nbh=1;q.nuw=1;q.nuh=1;q.nu1=1;
end 
if false
  nm1=cellstr_ls('gc/f6/i/4000/0*/2304/xxxxc');
  nm1{end+1}='gc/f6/i/4000/1100/2304/xxxxc';
  b=ded_gc_extrap_uw(nm1,'ray',[1 20],[0.1 35],showb,q);
  p=ded_read_param(nm1);
  ded_gc_fit_Re_hw1([p.Re],[b.guw],[b.gu1],[p.g],nm1);
  
  nm2=cellstr_ls('gc/f6/i/4000/1*/2304/xxxxc');
  nm2={nm2{1} nm2{3:end}};
  b=ded_gc_extrap_uw(nm2,'ray',[10 30],[0.1 35],showb,q);
  p=ded_read_param(nm2);
  ded_gc_fit_Re_hw1([p.Re],[b.guw],[b.gu1],[p.g],nm2);
  
  
  fns=cellstr_ls('gc/f6/i/*/*/*/*');
  ctp.typ={'X','g'};
  ctp.tol=[1e-3,1e-4];
  ded_make_avg(fns,'ay',ctp);
  ded_gc_extrap_uw(fns,'ray',[ 1 15],[0.1 35],showb,q1);
  
  p=ded_read_param(fns);
end


fns=cellstr_ls('gc/f6/i/*/*/2304/xxxx*');
ded_gc_fit_erf_profiles(fns,ctp,'ray',true,true)

p=ded_read_param(fns);

ctp.tol=[1e-3,1e-4];

q1.nbw=1;q1.nbh=1;q1.nuw=2;q1.nuh=2;q1.nu1=2;
q2.nbw=1;q2.nbh=1;q2.nuw=1;q2.nuh=1;q2.nu1=1;

f1=find([p.U1]< 0.5);b1=ded_gc_extrap_uw({fns{f1}},'ray',[ 1 30],[0.1 35],showb,q1);p1=[b1.param];
f2=find([p.U1]>=0.5);b2=ded_gc_extrap_uw({fns{f2}},'ray',[10 30],[0.1 35],showb,q2);p2=[b2.param];
figure;ded_gc_fit_Re_hw1([p1.Re],[b1.guw],[b1.gu1],[p1.g],{p1.name});
figure;ded_gc_fit_Re_hw1([p2.Re],[b2.guw],[b2.gu1],[p2.g],{p2.name});
  
  

clear('ctp');ctp.typ={'X','g'};
ctp.T=100;  
fns=cellstr_ls('gc/f6/mg/*/*/*/xxxxg');
ded_make_avg(fns,'ay',ctp);
q1.nbw=2;q1.nbh=1;q1.nuw=2;q1.nuh=1;q1.nu1=1;
dc;
b=ded_gc_extrap_uw(fns,'ray',[5 30],[0.1 35],showb,q1,'wh30')
p=ded_read_param({b.nm});
Re=[p.Re];
Pe=round([p.Re].*[p.Scb]);
f=find(Pe==2000 & Re<= 2600);
a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name});


f=find(Pe==1000 & Re<=2000);
a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name});

a=ded_gc_fit_Re_hw1([p.Re],[b.guh],[b.guw],[b.gu1],[p.g],{p.name});



nm='gc/f6/i/4000/1800/2304/xxxxe'
ded_gc_fit_erf_profiles(fns,ctp,'ray',true,true)
ded_gc_extrap_uw(nm,'ray',[ 1 20],[0.1 22],showb,q1)


cd('~/');
fns=cellstr_ls('gc/f6/i/*/*/*/xxxxg');
%ded_plot_X(fns);
ctp.T=30;
ded_make_avg(fns,'ay',ctp);
dc;
b=ded_gc_extrap_uw(fns,'ray',[ 10 20],[0.1 35],showb,q1,'wh30')
p=ded_read_param({b.nm});
f=find([p.Re]<=800);
a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name});

groupplot([p.Re],[b.guw],round([p.Re].*[p.Scb]));

groupplot([p.Re],[b.guw],round([p.Re].*[p.Scb]));



ded_plot_X('gc/f6/i/*/*/*/[0-9]*',-2);
ded_plot_X('gc/f6/i/8/*/*/*',-2);



Calculate values extrapolated to uh = 0.8







f=find(Pe==1000 & Re<= 1100);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'all');
f=find(Pe==1000 & Re>= 1100);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'all');
f=find(Pe==2000)            ;a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'all');


f=find([b.guw]<0.66);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'mid');


f=find(Pe==2000 & Re<= 500 & Re> 100);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'low');
R=[100 150 200 250 400 500];disp(sprintf('--wu %6.4f --u1 %6.4f %05.0f\n',[a.fuw(R);a.fu1(R);R]));
%disp(sprintf('--uw %6.4',fa.fuw([100 150 200 250 400 500]));

dc;f=find(Pe==2000 & Re<= 1000 & Re>=  600);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'poly3');
dc;f=find(Pe==2000 & Re<= 2000 & Re>= 1100);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'high');



f=find(Pe==2000 & Re<= 2500 & Re>= 600);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'all');

f=find(Pe==2000 & Re>= 1000);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'all');




clear('ctp'); ctp.typ={'X','g'};  ctp.tol=[1e-3,1e-4];
cd('~/');fns=cellstr_ls('gc/f6/*/[124]*/*/*/xxxxg');
ded_make_avg(fns,'ay',ctp);
q1.nbw=2;q1.nbh=1;q1.nuw=2;q1.nuh=1;q1.nu1=1;
dc;b=ded_gc_extrap_uw(fns,'ray',[5 30],[0.1 35],showb,q1,'wh30');
p=ded_read_param({b.nm});
Re=[p.Re];
Pe=round([p.Re].*[p.Scb]);
uw=[p.wu];
U1=[p.U1];
uh=[p.hu];
dc;
p2.mfc=[];

f=find([b.guw]<0.66 & [b.gu1]<0.8);a=ded_gc_fit_Re_hw1([p(f).Re],[b(f).guh],[b(f).guw],[b(f).gu1],[p(f).g],{p(f).name},'mid',Pe(f));
Re=unique([p.Re]);
Pe=unique([p.Pe]);
[Re Pe]=ndgrid(Re,Pe);
Re=Re(:)';Pe=Pe(:)';
disp(sprintf('--pfn  -rfn --reset --wu %6.4f --U1 %6.4f %04.0f/%04.0f/2304/xxxxh\n',[a.fuw(Re,Pe);a.fu1(Re,Pe);Pe;Re]));






end


clear('ctp'); 
ctp.typ={'X'};  
ctp.tol=1e-3;
cd('~/');fns=cellstr_ls('gc/f6/mg/*/*/*/[0-9]*');
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

figure;clf;[h lh uc]=groupplot(Re,[b.guw],Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('uw');groupplot(Re,uw,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,[b.gu1],Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('U1');groupplot(Re,U1,Pe,p2);
figure;clf;[h lh uc]=groupplot(Re,[b.guh],Pe);legend(lh,cellsprintf('%.0f',uc));ylabel('Uh');groupplot(Re,uh,Pe,p2);

ded_plot_X('gc/f6/mg/*/*/*/[0-9]*');
