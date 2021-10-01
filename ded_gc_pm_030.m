function ded_gc_pm_030(nm,alpha,suffix)
%rsync hamilton:gc/f6/f/09 -vaP ~/gc/f6/f --exclude "b*"   --exclude "check*"   --exclude "final*"  
%p=ded_3d_param(nm);;
%ded_gc_pm_030('pm/029',0.05,'3');
%ded_gc_pm_030('pm/029',0.02,'4');
%ded_gc_pm_030('pm/029',0.01,'5');


if nargin<2
  alpha=0.05;
end
if nargin<2
  suffix='2';
end


if nargin==0
  nms={'pm/030','pm/031','pm/032','pm/033','pm/034','pm/035','pm/036','pm/037'};
  for j=1:length(nms)
    ded_gc_pm_030(nms{j},alpha,suffix);
  end
  return
end
pp=ded_read_param(nm);
if isempty(pp)
  return;
end

p.text=[];
p.flipy= 0;
p.png= 1;
p.gdx= 1;
p.gdy= 1;
p.box= 1;
p.bottom= 1;
p.isocap= 'all';
p.X= NaN;
p.Y= NaN;
p.Z= NaN;
p.fminx= Inf;
p.fmaxx= -Inf;
p.fminy= Inf;
p.fmaxy= -Inf;
p.fminz= Inf;
p.fmaxz= -Inf;
p.dq= 10;
p.top= 1;
p.y0= 0;
p.z0= 0;
p.bh= 0.5000;
p.pcx= [2 1 3];
p.incremental= 0;
p.col= {'red'  'blue'};
p.maxfn= Inf;
p.nrot= 1;
p.zmax= 0.7000;
p.trackfront= 0;
p.rd= 'opengl';
p.tsn= 1;
p.za= 0;
p.ctar= 2;
p.cont= 1;

p.rgx= [0 28];
p.x0= 2;
p.cva = 4.5;

p.sz= 2*[480 512];
p.reduce= 1;
p.aa = 0;

p.sz= 2*[480 512];
p.reduce= 0;
p.aa = 4;



p.rgy=pp.W/2*[-1 1];
p.rgz=pp.H/2*[-1 1];
p.xtick=p.rgx(1):5:p.rgx(end);
p.ytick=p.rgy(1):5:p.rgy(end);
p.ztick=p.rgz(1):5:p.rgz(end);

p.rg= 2*[pp.S 0];
p.isoval=alpha*abs(pp.S);
p.ftype= 's';
p.nmbnm=[nm '-' suffix 's'];
ded_3d(nm,p);

p.rg= 2*[0 pp.B];
p.isoval= alpha*abs(pp.B);
p.ftype= 'b';
p.nmbnm=[nm '-' suffix 'b'];
ded_3d(nm,p);




