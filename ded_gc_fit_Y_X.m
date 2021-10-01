function [c f pp err]=ded_gc_fit_Y_X(nms,Ys,Xs,typ,p0)

nmsin=nms;
if nargin<4
  typ=[];
end
if nargin<5
  p0=[];
end
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
  if isempty(nms)
    disp(sprintf('ded_gc_fit_Y_X: no files found for %s',nmsin));
    disp([DD nmsin]);
    c=[];
    f=[];
    pp=[];
    return
  end
  

  
  j=0;
  for k=1:length(nms)
    nm=nms{k};
    s=ded_read_stats(nm);
    p=ded_read_param(nm);
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
    T=ded_convergence_T(nm);
    T=min(T,0.8*s.t(end));
    
    j=j+1;
    name{j}=p.name;
    f=find(s.t>=min(T,(min(s.t)+max(s.t))/2));
    nn={'g','X','U'};
    W(j)=1+s.t(end)-T;
    
    c.PIDX(j) = get_field(p,'PIDX');
    c.Re(j)   = get_field(p,'Re');
    c.Scb(j)  = get_field(p,'Scb');
    c.Scs(j)  = get_field(p,'Scs');
    c.U(j)    = get_field(p,'U');
    c.g(j)    = get_field(p,'g');
    c.X(j)    = get_field(p,'PIDX');
    for kk=1:length(nn)
      c.(['s' nn{kk}])(j)=0;
      if isfield(s,nn{kk})
        if max(f)<=length(s.(nn{kk}))
          c.(nn{kk})(j)=mean(s.(nn{kk})(f));
          c.(['s' nn{kk}])(j)=std(s.(nn{kk})(f));
        end
      end
    end
    if isfinite(c.PIDX(j))
         c.g(j)=p.g;
         c.U(j)=p.U;
         %disp(sprintf('ded_plot_PID: PID not active for %s',nm));
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
    c.sr{j}=p.series;
    c.nm{j}=nm;
  end
end
uNx  = unique(c.Nx);
ug   = unique(c.g);
uX   = unique(c.X);
uRe  = unique(c.Re);
uP   = unique(c.PIDX);
uH   = unique(c.H);
uhu  = unique(c.hu);
uhb  = unique(c.hb);
uU   = unique(c.U);

[X,Xstd]=ded_gc_get_var(c,Xs);
[Y,Ystd]=ded_gc_get_var(c,Ys);

lp=[];
pp=[];
Xmin = min(X);
Xmax = max(X);
Xrg  = Xmax-Xmin;
Ymin = min(Y);
Ymax = max(Y);
Yrg  = Ymax-Ymin;

if ~isempty(typ)
  switch(typ)
    case 'rat23'
      ff = @(p,x) polyval(p(1:2),x/Xmax)./max(eps,polyval([1 p(3:4)],x/Xmax));
      pp0=[0 mean(Y) 0 1];
      lp=repmat(-inf,size(pp0));
    case 'rat34'
      ff = @(p,x) polyval(p(1:3),x/Xmax)./max(eps,polyval([1 p(4:6)],x/Xmax));
      pp0=[0 0 mean(Y) 0 0 1];
      lp=repmat(-inf,size(pp0));
    case 'iexp'
      ff = @(p,x) p(1)+p(2)*1e4./x-p(3)*exp(-p(4)*x/1e4);
      pp0=[2 10 0.1 5];
      lp=repmat(-inf,size(pp0));
    case 'rat45'
      ff = @(p,x) polyval(p(1:4),x/Xmax)./max(eps,polyval([1 p(5:8)],x/Xmax));
      pp0=[0 0 0 mean(Y) 0 0 0 1];
      lp=repmat(-inf,size(pp0));
    case 'linear'
      ff = @(p,x) polyval(p,x);
      pp=polyfit(X,Y,1);
    case 'quadratic'
      ff = @(p,x) polyval(p,x);
      pp=polyfit(X,Y,2);
   case 'cubic'
      ff = @(p,x) polyval(p,x,3);
      pp=polyfit(X,Y,3);
    case 'hpow'
      ff = @(p,x) p(1)+p(2)./sqrt(x);
      pp0=[Ymin (Yrg)*sqrt(Xmin)];
    case 'dpow'
      ff = @(p,x) p(1)+p(2)*(x/Xrg).^(-p(3));
      pp0(1)=Ymin;
      pp0(2)=Yrg;
      pp0(3)=1;
   case 'gpow'
      ff = @(p,x) p(1)+p(2)*max((x-Xmin)/Xrg+p(4),0).^(-p(3));
      pp0(1)=Ymin;
      pp0(2)=Yrg;
      pp0(3)=1;
      pp0(4)=0.1;
    case 'exp'
      ff = @(p,x) p(1)+p(2)*exp(-p(3)*(x-Xmin)/Xrg);
      pp0(1)=Ymin;
      pp0(2)=Yrg;
      pp0(3)=1;
    case 'd2pow' % Discontinuous double power fits
      f3 = @(x) min(1,max(0,x)).^2.*(3-2*min(1,max(0,x)));
      ff = @(p,x) (p(1)+p(2)*(x/1e4).^(-p(3))).*f3(p(8)*(p(4)-x/1e4)) + (p(5)+p(6)*(x/1e4).^(-p(7))).*f3(1-p(8)*(p(4)-x/1e4));
      pp0([1 5])=Ymin;
      pp0([2 6])=(Yrg)*.5;
      pp0([3 7])=1;
      pp0(4)=1;
      pp0(8)=10;
    case 'd3pow' % Continuous power fits
      ff = @(p,x) p(1)-p(2)*(1-min(1,x/p(4)/1e4).^(-p(3))) + p(5)*(max(1,x/p(4)/1e4).^(-p(6))-1);
      pp0(1)=Ymin;
      pp0([2 5])=(Yrg)*.5;
      pp0([3 6])=1;
      pp0(4)=1;
    case 'ipow'
      ff = @(p,x) p(1)+p(2)*x.^(p(3));
      pp0(1)=Ymin;
      pp0(2)=Yrg;
      pp0(3)=1;
    otherwise
      error(sprintf('Unknown fit type %s',typ))
  end
  
  
  if isempty(pp)
    opt=optimoptions('lsqnonlin','display','none','tolx',1e-8,'tolfun',1e-8);  
    if ~isempty(p0)
      pp0=p0;
    end
    if isempty(lp)
      lp=0*pp0;
    end
    f = @(p) W.*(ff(p,X)-Y);
    pp=lsqnonlin(f,pp0,lp,[],opt);
  end
  xx=linspace(Xmin,Xmax,1e3);
  yy=ff(pp,xx);
  disp([typ sprintf(' %4.3f',pp)]);
  disp(sprintf('%20s %9s %9s %9s','name',Xs,Ys,Ys));
  [XXX XXf]=sort(X);
  for k=1:length(X)
    j=XXf(k);
    disp(sprintf('%20s %9.4f %9.4f %9.4f',name{j},X(j),Y(j),ff(pp,X(j))));
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
  err=Y-f(pp,X);
end
legend(lh,uc,'location','ne');


return
nms='gc/qgc*/*';
figure(1);c=ded_gc_fit_Y_X(nms,'g','Re');
figure(2);ded_gc_fit_Y_X(c,'X','Re');

