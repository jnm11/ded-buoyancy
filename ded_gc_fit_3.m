function [c f pp]=ded_gc_fit_3(nms,fnm,typ,mm,p0)

if nargin<4
  typ=[]
end
if nargin<5
  mm=[];
end
if nargin<6
  p0=[];
end
%rsync -vaP hamilton:gc/gc2d7n ~/gc --exclude "final*" --exclude "check*" --exclude "b*" --exclude "force*" 
%rsync -vaP hamilton:gc/f6/g ~/gc/f6 --exclude "final*" --exclude "check*" --exclude "b*" --exclude "force*" 
nms='gc/gc2d7n/*';
fnm={'Re','g','X'};
mm1=[5 2];
f=[];
if isstruct(nms)
  c=nms;
  nms={c.nm};
else
  DD=[ded_dedalus_data_dir '/'];
  if ~iscell(nms)
    if any(nms=='*')
      nms=cellstr_ls([DD nms ],'ls -d');
    else
      if ~isfile(nms)
        nms=[DD nms];
      end
      nms={nms};
    end
  end
  
  j=0;
  for k=1:length(nms)
    nm=nms{k};
    s=ded_read_stats(nm);
    p=ded_read_param(nm);
    T=ded_convergence_T(nm);
    
    if isempty(s)
      disp(sprintf('ded_plot_PID: %s could not read stats',nm));
      continue
    end
    if isempty(s.t)
      disp(sprintf('ded_plot_PID: %s no time in stats',nm));
      continue
    end
    if isempty(p)
      disp(sprintf('ded_plot_PID: %s could not read param',nm));
      continue
    end
    
    j=j+1;
    f=find(s.t>=min(T,(min(s.t)+max(s.t))/2));
    nn={'g','X','U'};
    c.P(j)=p.PIDX;
    c.Re(j)=p.Re;
    c.U(j)=p.U;
    c.g(j)=p.g;
    c.X(j)=p.PIDX;
    for kk=1:length(nn)
      c.(['s' nn{kk}])(j)=0;
      if isfield(s,nn{kk})
        if max(f)<=length(s.(nn{kk}))
          c.(nn{kk})(j)=mean(s.(nn{kk})(f));
          c.(['s' nn{kk}])(j)=std(s.(nn{kk})(f));
        end
      end
    end
    if p.PIDX==0
         c.g(j)=p.g;
         c.U(j)=p.U;
%	          disp(sprintf('ded_plot_PID: PID not active for %s',nm));
    end
   if isfield(p,'hu')
      c.hu(j)=p.hu;
    else
      c.hu(j)=p.h;
    end      
    if isfield(p,'hb')
      c.hb(j)=p.hb;
    else
      c.hb(j)=p.h;
    end      
    c.Nx(j)=p.Nx;
    c.Ny(j)=p.Ny;
    c.Nz(j)=p.Nz;
    c.L(j)=p.L;
    c.H(j)=p.H;
    c.W(j)=p.W;
    c.N(j)=p.noise;
    c.sr{j}=p.series;
    c.nm{j}=nm;
  end
end
uNx  = unique(c.Nx);
ug   = unique(c.g);
uX   = unique(c.X);
uRe  = unique(c.Re);
uP   = unique(c.P);
uH   = unique(c.H);
uhu  = unique(c.hu);
uhb  = unique(c.hb);
uU   = unique(c.U);

for j=1:length(fnm)
  switch(fnm{j})
    case('Nx'); X{j}=c.Nx;
    case('g');  X{j}=c.g;
    case('X');  X{j}=c.X;
    case('R');  X{j}=c.R;
    case('Re'); X{j}=c.Re;
    case('sg'); X{j}=c.sg;
    case('sX'); X{j}=c.sX;
    case('sR'); X{j}=c.sR;
    case('P');  X{j}=c.P;
    case('H');  X{j}=c.H;
    case('hu'); X{j}=c.hu;
    case('hb'); X{j}=c.hb;
    case('U');  X{j}=c.U;
  end
end


f=find(X{3}>18 & X{3}<22);
X{1}=X{1}(f);
X{2}=X{2}(f);
X{3}=X{3}(f);

groupplot(X{1},X{3}); 

uX1=unique(X{1});
uX2=unique(X{2});
X1=uX1(round(linspace(1,length(uX1)),mm(1)));
X2=uX2(round(linspace(1,length(uX2)),mm(2)));
rg1=1:mm(1);
rg2=mm(1)+(1:mm(2));

f1==@(p,x) interp1(X1,p(rg1),x,'pchip');
f2==@(p,x) interp1(X2,p(rg2),x,'pchip');
f3==@(p,x) interp1(p(rg1),X1,x,'pchip');


X = f(Re,g);   %20=f1(Re,g0)   X = f1(Re,g0)*(1+g2(Re,g0)*(g-g0)); 

f = @(p,x,y) f1(p(rg1),x)*( f2(p(rg2),y)      )

p0(rg2)=0;
opt=optimoptions('lsqnonlin','display','none','tolx',1e-8,'tolfun',1e-8);
pp=lsqnonlin(f,pp,lp,[],opt);


xx=linspace(min(X),max(X),1e3);
  yy=ff(pp,xx);
  disp([typ sprintf(' %4.3f',pp)]);
  for j=1:length(X)
    disp(sprintf('%.0f %4.3f',X(j),ff(pp,X(j))));
  end
  f = ff;
else
  xx=[];
end

gpp.lt='-';
[h lh uc]=groupplot(X,Y,c.sr,gpp);
xlabel(Xs);
ylabel(Ys);
if ~isempty(xx)
  line(xx,yy);
end
legend(lh,uc,'location','ne');


return
nms='gc/qgc*/*';
figure(1);c=ded_gc_fit_Y_X(nms,'g','Re');
figure(2);ded_gc_fit_Y_X(c,'X','Re');

