function CPID=ded_gc_pid_tune(nm,X,p,mint)
if nargin<3
  p=[];
end
if nargin<4
  mint=[];
end
if isempty(p)
  p=0.5;
end
if isempty(mint)
  mint=0;
end

[sys,b]=ded_gc_calc_larx(nm,false,mint);

% $$$ Y=mean(b.Y((abs(b.X-X)<0.2)));
% $$$ if isempty(Y) | ~isfinite(Y)
% $$$   Y=b.Y(findmin(abs(b.X-X)));
% $$$ end

%nn=sys.na+sys.nb;

%lsys = linearize(sys,X,repmat(Y,nn,1));
opt = pidtuneOptions('PhaseMargin',70,'DesignFocus','reference-tracking');
[CPID,info]=pidtune(sys,'PIDF',0.5,opt);
% 1 selects a fast response time
%TPID = feedback(CPID*lsys, 1);
%step(TPID);
%disp(CPID);
%disp(info);
disp(sprintf('X=%5.2f, Kp=%8.5f, Ki=%8.5f, Kd=%8.5f',X,CPID.Kp,CPID.Ki,CPID.Kd));

return;
nm='gc/610';
ded_gc_calc_larx(nm); 
CPID=ded_gc_pid_tune(nm,20);



CPID=ded_gc_pid_tune(nm,6);
Kp: -0.8776
Ki: -0.1708
Kd: -0.0435
 

nms={'gc/600','gc/606','gc/622','gc/644','gc/650','gc/601','gc/607','gc/631','gc/645',...
     'gc/651','gc/602','gc/608','gc/640','gc/646','gc/652','gc/603','gc/641','gc/647',...
     'gc/604','gc/642','gc/648','gc/605','gc/643','gc/649'};
sys=ded_gc_calc_larx(nms);

nm='gc/4';
ded_gc_calc_larx(nm,[],.5); 
CPID=ded_gc_pid_tune(nm,6);
             Kp: -0.7809
              Ki: -0.1234
              Kd: -0.0387
 
CPID=ded_gc_pid_tune(nm,7);
             Kp: -0.7809
              Ki: -0.1234
              Kd: -0.0387
 


Kp: -0.7996
Ki: -0.0805
Kd: -0.2406

nm='gc/3';
ded_gc_calc_larx(nm,[],50); 
CPID=ded_gc_pid_tune(nm,6);
Kp: -0.8124
Ki: -0.0801
Kd: -0.2686

A="--PIDP 0.1 --PIDI 0.01 --PIDD 0.03 --PIDX 6 --PIDW 1"
mpiexec -n 8 ded_gc-1.31.py --sType gc -r -p  $A 3

A="--PIDP 0.1 --PIDI 0.02 --PIDD 0.04 --PIDW 1"
mpiexec -n 8 ded_gc-1.31.py --sType gc --rfn 3 --pfn 3 --PIDX 7 $A --reset 4
mpiexec -n 8 ded_gc-1.31.py --sType gc --rfn 4 --pfn 4 --PIDX 6 $A --reset 3

A="--PIDP 0.5 --PIDI 0.01 --PIDD 0.04 --PIDW 1"
mpiexec -n 8 ded_gc-1.31.py --sType gc --rfn 4 --pfn 4 --PIDX 7 $A --reset 5
mpiexec -n 8 ded_gc-1.31.py --sType gc --rfn 5 --pfn 5 --PIDX 5 $A --reset 6
mpiexec -n 8 ded_gc-1.31.py --sType gc --rfn 6 --pfn 6 --PIDX 7    --reset 7

nm='gc/1';sys=ded_gc_calc_larx(nm,[],1);CPID=ded_gc_pid_tune(nm,6);
nm='gc/2';sys=ded_gc_calc_larx(nm,[],1);CPID=ded_gc_pid_tune(nm,6);
nm='gc/3';sys=ded_gc_calc_larx(nm,[],1);CPID=ded_gc_pid_tune(nm,6);
nm='gc/4';sys=ded_gc_calc_larx(nm,[],1);CPID=ded_gc_pid_tune(nm,6);
nm='gc/5';sys=ded_gc_calc_larx(nm,[],1);CPID=ded_gc_pid_tune(nm,6);

ded_gc_calc_larx('gc/4');


[sys,lsys,B,fnm,s]=ded_gc_calc_larx({'gc/3','gc/4'},[],1); 
CPID=ded_gc_pid_tune(sys,6);

k=1e4;
clf;hold('on');
Y=linspace(0,max(B.Y),100);
W=0*Y;
warning('off');
for j=1:length(YY);
  WW = sim(sys,repmat(Y(j),k,1));
  W(j)=WW(end);
  if j>1
    if isfinite(W(j-1)) & ~isfinite(W(j))
      break;
    end
  end
end
W=W(1:j-1);
Y=Y(1:j-1);




Z=s(end,1,[]); 

lsys = linearize(sys,X,repmat(b.Y(end),4,1));

C6=ded_gc_pid_tune(sys,6);

sys=ded_gc_calc_larx('gc/618',[],1);
CPID=ded_gc_pid_tune('gc/618',20);
Kp: -0.8928
Ki: -0.2091
Kd: -0.0441
 

lsys = linearize(sys,20,repmat(b.Y(end),4,1));


dc;
nn=610:618;
for j=1:length(nn)
  nm=sprintf('gc/%03u',nn(j));
  ded_gc_calc_larx(nm); 
  CPID=ded_gc_pid_tune(nm,20);
  Kp(j)=-CPID.Kp;
  Ki(j)=-CPID.Ki;
  Kd(j)=-CPID.Kd;
  drawnow
end

disp([Kp(:) Ki(:) Kd(:)]);

nm='gc/10';
dc;ded_gc_calc_larx(nm,[],1,[],[3 3 1]);C1=ded_gc_pid_tune(nm,5);C2=ded_gc_pid_tune(nm,6);C3=ded_gc_pid_tune(nm,7);

pidTuner(sys);


%pidTuner(sys);
%ded_gc-1.31.py -r -p               --PIDP 0.80 --PIDI 0.056 --PIDD 0.04 --PIDST 25  gc/499
%ded_gc-1.31.py --rfn 499 --pfn 499 --PIDP 0.42 --PIDI 0.022 --PIDD 0.021 gc/498
%mpiexec -n 32 ded_gc-1.31.py -r -p --PIDP 0.40 --PIDI 0.020 --PIDD 0.020 --reset gc/498

nm='gc/498';
dc;ded_gc_calc_larx(nm,[],30);C1=ded_gc_pid_tune(nm,19.5);C2=ded_gc_pid_tune(nm,20);C3=ded_gc_pid_tune(nm,20.5);
