function ded_process_gc_f7_mgi

showb=true;
ctp.typ={'X','g'};
ctp.tol=[1e-3,1e-4];
ctp.mint=30;
q.nbw=2;q.nbh=1;q.nuw=2;q.nuh=1;q.nu1=1;
dd=[ded_dedalus_data_dir '/results/gc/f6'];
cd('~/');

fn1=cellstr_ls('gc/f6/mg/1000/*/2304/[0-9]*');
fn2=cellstr_ls('gc/f6/mg/2000/*/2304/[0-9]*');
fn3=cellstr_ls('gc/f6/mg/4000/*/2304/[0-9]*');
fn4=cellstr_ls('gc/f6/mg/2000/*/4608/[0-9]*');
fn5=cellstr_ls('gc/f6/mg/4000/*/4608/[0-9]*');
fn6=cellstr_ls('gc/f6/mg/8000/*/4608/[0-9]*');
fns=cat(1,fn1,fn2,fn3,fn4,fn5,fn6);
fnt=cellstrpostfix('/time',fns);
fnmat=[dd '/mgi.mat'];
fnc=[dd '/mgi-classify.mat'];
fnr=[dd '/mgi-RJ.mat'];
p=ded_read_param(fns);
ded_make_param(fns);
ded_make_avg(fns,'ay',ctp);
ded_gc_fit_erf_profiles(fns,[],'ray',false,false);
ded_gc_stability(fns);
b=ded_gc_extrap_uw(fns,'ray',[5 30],[0.1 35],showb,q,'wh30');
p=ded_read_param({b.nm});

display=false;
save(fnmat','b','p');
xmin=1;
a=ded_gc_f7_g_classify(fns,xmin,display);save(fnc,'a');
a=ded_gc_f7_g_RJ(fns,xmin,'front',display);save(fnr,'a');


% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ showb=true;
% $$$ ctp.typ={'X','g'};
% $$$ ctp.tol=[1e-3,1e-4];
% $$$ ctp.mint=30;
% $$$ q.nbw=2;q.nbh=1;q.nuw=2;q.nuh=1;q.nu1=1;
% $$$ dd=[ded_dedalus_data_dir '/results/gc/f6'];
% $$$ cd('~/');
% $$$ ded_make_avg('gc/f6/i/8000/0500/4608/27264','ay',ctp);
% $$$ 
% $$$ ded_make_avg('gc/f6/i/8000/0700/4608/23904','ay',ctp);
