function [c f pp]=ded_gc_fit_Z_Y_X(nms,Zs,Ys,Xs,typ,p0)

nmsin=nms;
if nargin<5
  typ=[];
end
if nargin<6
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
    disp(sprintf('ded_gc_fit_Z_Y_X: no files found for %s',nmsin));
    disp([DD nmsin]);
    c=[];
    f=[];
    pp=[];
    return
  end
  
  clf; 
  
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
    c.Re(j)=p.Re;
    c.V(j)=p.V;
    c.g(j)=p.g;
    if isfield(s,'X')
      c.X(j)=mean(s.X(f));
      c.sX(j)=std(s.X(f));
    end
    if isfield(s,'g')
      c.g(j)=mean(s.g(f));
      c.sg(j)=std(s.g(f));
    end
    if isfield(s,'U')
      c.V(j)=mean(s.V(f));
      c.Vg(j)=std(s.V(f));
    end
    tnm=cellstrtok(p.name,'/');
    
    fldnms={'hu','hb','Nx','Ny','Nz','L','H','W','Scb','Scs','Peb','Pes','fldnms'};
    for kk=1:length(fldnms)
      c.(fldnms{kk})(j)= get_field(p,fldnms{kk});
    end
    c.sr{j}=p.series;
    c.nm{j}=nm;
    for kk=1:length(tnm)
      c.(['nm' num2str(kk)]){j}=tnm{kk};
    end
    for kk=1:length(tnm)
      c.(['nm' 'a'-1+kk]){j}=tnm{length(tnm)+1-kk};
    end
  end
end
if isfield(c,'nm1');unm1=unique(c.nm1);end
if isfield(c,'nm2');unm2=unique(c.nm2);end
if isfield(c,'nm3');unm3=unique(c.nm3);end
if isfield(c,'nm4');unm4=unique(c.nm4);end
if isfield(c,'nm5');unm5=unique(c.nm5);end
if isfield(c,'nma');unma=unique(c.nma);end
if isfield(c,'nmb');unmb=unique(c.nmb);end
if isfield(c,'nmc');unmc=unique(c.nmc);end
if isfield(c,'nmd');unmd=unique(c.nmd);end
if isfield(c,'nme');unme=unique(c.nme);end
uNx  = unique(c.Nx);
ug   = unique(c.g);
uX   = unique(c.X);
uRe  = unique(c.Re);
uScb  = unique(c.Scb);
uScs  = unique(c.Scs);
uH   = unique(c.H);
uhu  = unique(c.hu);
uhb  = unique(c.hb);
uU   = unique(c.V);

X=ded_gc_get_var(c,Xs);
Y=ded_gc_get_var(c,Ys);
Z=ded_gc_get_var(c,Zs);

lp=[];
pp=[];
if ~isempty(typ)
  switch(typ)
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
      pp0=[min(Y) (max(Y)-min(Y))*sqrt(min(X))];
    case 'dpow'
      ff = @(p,x) p(1)+p(2)*(x/1e4).^(-p(3));
      pp0(1)=min(Y);
      pp0(2)=max(Y)-min(Y);
      pp0(3)=1;
    case 'd2pow' % Discontinuous double power fits
      f3 = @(x) min(1,max(0,x)).^2.*(3-2*min(1,max(0,x)));
      ff = @(p,x) (p(1)+p(2)*(x/1e4).^(-p(3))).*f3(p(8)*(p(4)-x/1e4)) + (p(5)+p(6)*(x/1e4).^(-p(7))).*f3(1-p(8)*(p(4)-x/1e4));
      pp0([1 5])=min(Y);
      pp0([2 6])=(max(Y)-min(Y))*.5;
      pp0([3 7])=1;
      pp0(4)=1;
      pp0(8)=10;
    case 'd3pow' % Continuous power fits
      ff = @(p,x) p(1)-p(2)*(1-min(1,x/p(4)/1e4).^(-p(3))) + p(5)*(max(1,x/p(4)/1e4).^(-p(6))-1);
      pp0(1)=min(Y);
      pp0([2 5])=(max(Y)-min(Y))*.5;
      pp0([3 6])=1;
      pp0(4)=1;
    case 'ipow'
      ff = @(p,x) p(1)+p(2)*x.^(p(3));
      pp0(1)=min(Y);
      pp0(2)=max(Y)-min(Y);
      pp0(3)=1;
    otherwise
      error(sprintf('Unknown fit type %s',typ))
  end
  

  if isempty(pp)
    opt=optimoptions('lsqnonlin','display','none','tolx',1e-8,'tolfun',1e-8);  if ~isempty(p0)
      pp0=p0;
    end
    if isempty(lp)
      lp=0*pp0;
    end
    f = @(p) ff(p,X)-Y;
    pp=lsqnonlin(f,pp0,lp,[],opt);
  end
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
[h lh uc]=groupplot(X,Y,Z,gpp);
xlabel(Xs);
ylabel(Ys);
if ~isempty(xx)
  line(xx,yy);
end
if ~iscell(uc)
  minuc=floor(log10(min(uc)));
  maxuc=ceil(log10(max(uc)));
  
  if maxuc-minuc<=2
    fmt=sprintf('%%%d.%df',6,max(0,5-maxuc));
  else
    fmt='%6.1e';
  end
  uc=cellsprintf(fmt,uc);
end
if length(lh<12)
  llh=legend(lh,uc,'location','best');
  title(llh,Zs);
end
return
