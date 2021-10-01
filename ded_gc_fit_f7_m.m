function [c f pp]=ded_gc_fit_f7_m(nms,XX)


Ys='g';
Xs='X';
nmsin=nms;
f=[];
if isstruct(nms)
  c=nms;
  nms={c.nm};
else
  nmsin=nms;
  nms=ded_parse_nms(nms);
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
    W(j)=100+s.t(end)-T;
    
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

for j=1:length(uRe)
  f=find([c.Re]==uRe(j));
  WW(j)=sum(W(f));
  XXX=c.X(f);
  ggg=c.g(f);
  if max(XXX)>XX-2 & min(XXX)<XX+2 & length(XXX)>1
    np=min(2,length(XXX)-1);
    while(1)
      %G=linspace(2.6,2.8,1e3);plot(G,polyval(p,G),ggg,XXX-XX,'s')
      p=polyfit(ggg,XXX-XX,np);
      dp=(np:-1:1).*p(1:end-1);
      rr=roots(p);
      rr=min(max(ggg),max(min(ggg),rr(imag(rr)==0)));
      ff=findmin(abs(rr-mean(ggg)));
      if isfinite(ff);break;end
      np=1;
    end
        
    Sg(j)=rr(ff);
    dXdg(j)=polyval(dp,Sg(j));
  else
    ff=findmin(abs(XX-XXX));
    Sg(j)=ggg(ff);
    dXdg(j)=NaN;
    disp([sprintf('ded_gc_fit_g: no sign change Re=%5.0f g=',uRe(j)) sprintf('%5.3f ',sort(ggg)) sprintf('%5.3f ',sort(XXX))]);
  end
end

X=uRe;
Y=Sg;

lp=[];
pp=[];
Xmin = min(X);
Xmax = max(X);
Xrg  = Xmax-Xmin;
Ymin = min(Y);
Ymax = max(Y);
Yrg  = Ymax-Ymin;

WW=1;
fs = @(x,y) (abs(x)<abs(y)).*x + (abs(x)>=abs(y)).*y;
%ff1 = @(p,x) (p(1)*Xmax/100+p(2)*x)./x;
ff1 = @(p,x) polyval(p(1:3),x/100)./polyval([1 0],x/100);
%ff2 = @(p,x) polyval(p(3:5),x/Xmax)./polyval([1 p(6:7)],x/Xmax);
%ff2 = @(p,x) polyval(p(4:end),log(x));
ff2 = @(p,x) polyval(p(4:6),log(x))./polyval([1 p(7:8)],log(x));
ff2 = @(p,x) polyval(p(4:7),log(x))./polyval([1 p(8:10)],log(x));
ff = @(p,x) ff1(p,x).*(x<175) + ff2(p,x).*(x>175);
p0=[2.103 0.611 2.217 -0.206 0.059 -0.078 0.026];
f1=find(X<=175);
f2=find(X>175);
p0(1:3)=polyfit(X(f1)/100,Y(f1).*X(f1)/100,2);
%p0(4:9)=polyfit(log(X(f2)),Y(f2),5);
p0(4:8)=[2.241 -30.797 109.206 -13.729 48.744];
p0(4:10)=[2.186 -41.643 261.772 -535.108 -19.094 120.375 -246.891];
lp=repmat(-inf,size(p0));
opt=optimoptions('lsqnonlin','display','none','tolx',1e-8,'tolfun',1e-8);  
f = @(p) WW.*fs(ff1(p,X)-Y,ff2(p,X)-Y); %W
pp=lsqnonlin(f,p0,lp,[],opt);

xx=linspace(0,1e4,1e3);
disp(sprintf(' %4.3f',pp));
if 0
  disp(sprintf('%20s %9s %9s','name',Xs,Ys));
  [XXX XXf]=sort(X);
  yy=[ff1(pp,xx(:)) ff2(pp,xx(:))];
  for k=1:length(X)
    j=XXf(k);
    disp(sprintf('%20s %9.4f %9.4f %9.4f %9.4f',name{j},X(j),Y(j),ff1(pp,X(j)),ff2(pp,X(j))));
  end
end

ming=min(c.g);
maxg=max(c.g);
grg=[ming maxg];
Rerg=[0 max(c.Re)*1.05];
gpp.lt='-';
clf;

grg1=[2.44    2.75];
grg2=[2.08    2.22];
Rerg1=[90 175];
Rerg2=[180/1.1 6000*1.1];
Xrg=[min(c.X) max(c.X)];
Xrg=[7.4 8.6];

clf;
jh=jsubplot([2 5],[0.05 0.05],[0.05 0.05],[0.02 0.02]);
set(jh,'box','on');

axes(jh(1,1));groupplot(c.g,c.X,c.Re,gpp);axis([grg1 Xrg]);
axes(jh(2,1));groupplot(c.g,c.X,c.Re,gpp);axis([grg2 Xrg]);

axes(jh(1,2));groupplot(c.Re,c.X,c.g,gpp);axis([Rerg1 Xrg]);
axes(jh(2,2));groupplot(c.Re,c.X,c.g,gpp);axis([Rerg2 Xrg]);


axes(jh(1,3));groupplot(c.Re,c.g,c.Re,gpp);hold('on');line(xx,ff1(pp,xx(:)),'linewidth',2);axis([Rerg1 grg1]);
axes(jh(2,3));groupplot(c.Re,c.g,c.Re,gpp);hold('on');line(xx,ff2(pp,xx(:)),'linewidth',2);axis([Rerg2 grg2]);
f1=find(uRe<=175);
f2=find(uRe>175);

axes(jh(1,4));plot(uRe(f1),Sg(f1)-ff1(pp,uRe(f1)),'s-',Rerg1,[0 0]);axis([Rerg1 -0.01 0.01]);
axes(jh(2,4));plot(uRe(f2),Sg(f2)-ff2(pp,uRe(f2)),'s-',Rerg2,[0 0]);axis([Rerg2 -0.01 0.01]);



axes(jh(1,5));plot(uRe,dXdg,'s-');axis([Rerg1 0  80]);
axes(jh(2,5));plot(uRe,dXdg,'s-');axis([Rerg2 0 200]);
set(jh(:,2:5),'xscale','log')

fRe1=ff1(pp,uRe);
fRe2=ff2(pp,uRe);
disp(sprintf('%7s %7s %7s %7s %7s %7s','Re','g',     'g',    'g',     'g',    'g'));
disp(sprintf('%7s %7s %7s %7s %7s %7s','Re','extrap','fit 1','round', 'fit2', 'round'));

nnm={};
g=[];
for j=1:length(uRe)
  if uRe(j)<175
    F=@(x) round(x*200)/200;
    gg=F(fRe1(j));
  else           
    F=@(x) round(x*400)/400;
    gg=F(fRe2(j)); 
  end
  if uRe(j)<1000 
    dd=sprintf('~/gc/f7/m/%05.0f/%05.0f',uRe(j),10000*gg);
  else            
    dd=sprintf('~/gc/f7/jb/%05.0f/%05.0f',uRe(j),10000*gg);
  end
  if ~isdir(dd)
    nnm{end+1}=dd(3:end);
    g(end+1)=gg;
  end
  disp(sprintf('%7.0f %7.4f %7.4f %7.4f %7.4f %7.4f',uRe(j),Sg(j),fRe1(j),F(fRe1(j)),fRe2(j),F(fRe2(j)) ));
end
for j=1:length(nnm)
  disp(sprintf('-g %7.4f %20s',g(j),nnm{j}));
end
