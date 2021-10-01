function [c ff pp]=ded_gc_fit_g_Re(nms,XX,typ,p0)


Ys='g';
Xs='X';
nmsin=nms;
if nargin<3
  typ=[];
end
if nargin<4
  p0=[];
end
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
    if ~isfield(s,'t')
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
    nn={'g','X','V'};
    W(j)=100+s.t(end)-T;
    
    c.PIDX(j) = get_field(p,'PIDX');
    c.Re(j)   = get_field(p,'Re');
    c.Scb(j)  = get_field(p,'Scb');
    c.Scs(j)  = get_field(p,'Scs');
    c.V(j)    = get_field(p,'V');
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
      c.V(j)=p.V;
      %disp(sprintf('ded_plot_PID: PID not active for %s',nm));
    end
    if isfield(p,'hu')
      c.hu(j)=p.hu;
    elseif isfield(p,'h')
      c.hu(j)=p.h;
    end      
    if isfield(p,'hb')
      c.hb(j)=p.hb;
    elseif isfield(p,'h')
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
if isfield(c,'hu'); uhu  = unique(c.hu); end
if isfield(c,'hb'); uhb  = unique(c.hb); end
uU   = unique(c.V);

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
      rr=rr(imag(rr)==0);
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
ff1=[];
ff2=[];
if ~isempty(typ)
  switch(typ)
    case 'poly2'
      pp=polyfit(X/Xmax,Y,2);
      ff = @(p,x) polyval(p,x/Xmax);
    case 'poly3'
      pp=polyfit(X/Xmax,Y.*Xmax./X,3);
      ff = @(p,x) x/Xmax.*polyval(p,x/Xmax);
     case 'poly4'
      pp=polyfit(X/Xmax,Y.*Xmax./X,4);
      ff = @(p,x) x/Xmax.*polyval(p,x/Xmax);
    case 'rat12'
      ff = @(p,x) polyval(p(1),x/Xmax)./max(eps,polyval([1 p(2)],x/Xmax));
      pp0=[mean(Y) 0];
      lp=repmat(-inf,size(pp0));
    case 'rat23'
      ff = @(p,x) polyval(p(1:2),x/Xmax)./max(eps,polyval([1 p(3:4)],x/Xmax));
      pp0=[0 mean(Y) 0 1];
      lp=repmat(-inf,size(pp0));
    case 'rat34'
      ff = @(p,x) polyval(p(1:3),x/Xmax)./max(eps,polyval([1 p(4:6)],x/Xmax));
      pp0=[40.5393 13.4220 0.3022 16.5469  7.2921 0.0717];
      lp=repmat(-inf,size(pp0));
    case 'split'
      %ff1 = @(p,x) (p(1)*Xmax+p(2)*x)./x;
      %ff2 = @(p,x) polyval(p(3:5),x/Xmax)./polyval([1 p(6:7)],x/Xmax);
      ff = @(p,x) (p(1)*Xmax/100+p(2)*x)./x.*(x<175) + polyval(p(3:5),x/Xmax)./polyval([1 p(6:7)],x/Xmax).*(x>175);
      pp0=[1.1 2.067 2.386046511, .4046511622, 0.5813953488e-1, ...
           .333333333, 0.2325581395e-1];
      pp0=[1.015 2.104 2.209 -0.331 0.064 -0.138 0.029];
      %pp0=[0.009 2.155 -28896.234 -61860.153 -18181.073 -37302.999 -7539.860];
      lp=repmat(-inf,size(pp0));
    case 'Ref'
            ff = @(p,x) polyval(p(1:3),x/Xmax)./max(eps,polyval([1 p(4:6)],x/Xmax))+...
                 polyval([p(8:9) 0 0],x/Xmax).*exp(-x/Xmax/p(7)) ;
      pp0=[40.5393 13.4220 0.3022 16.5469  7.2921 0.0717 0.0610   0 10.9936   -1.1153];
      lp=repmat(-inf,size(pp0));
      lp(7)=0.05;
    case 'rat45'
      ff = @(p,x) polyval(p(1:4),x/Xmax)./max(eps,polyval([1 p(5:8)],x/Xmax));
      pp0=[0 0 0 mean(Y) 0 0 1 0];
      lp=repmat(-inf,size(pp0));
    case 'iexp'
      ff = @(p,x) p(1)+p(2)*1e4./x-p(3)*exp(-p(4)*x/1e4);
      pp0=[2 10 0.1 5];
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
    fs = @(x,y) (abs(x)<abs(y)).*x + (abs(x)>=abs(y)).*y;
    if isempty(ff1)
      f = @(p) WW.*(ff(p,X)-Y); %W
    else
      f = @(p) WW.*fs(ff1(p,X)-Y,ff2(p,X)-Y); %W
    end
    pp=lsqnonlin(f,pp0,lp,[],opt);
  end
  xx=linspace(Xmin,Xmax,1e3);
  disp([typ sprintf(' %4.3f',pp)]);
  disp(sprintf('%20s %9s %9s','name',Xs,Ys));
  [XXX XXf]=sort(X);
  if isempty(ff1)
    yy=ff(pp,xx);
    for k=1:length(X)
      j=XXf(k);
      disp(sprintf('%20s %9.4f %9.4f %9.4f',name{j},X(j),Y(j),ff(pp,X(j))));
    end
  else
    yy=[ff1(pp,xx(:)) ff2(pp,xx(:))];
    for k=1:length(X)
      j=XXf(k);
      disp(sprintf('%20s %9.4f %9.4f %9.4f %9.4f',name{j},X(j),Y(j),ff1(pp,X(j)),ff2(pp,X(j))));
    end
    f = ff;
  end
else
  xx=[];
end
ming=min(c.g);
maxg=max(c.g);
grg=[ming maxg]+0.05*(maxg-ming)*[-1 1];
minRe=min(c.Re);
maxRe=max(c.Re);
Rerg=[minRe maxRe]+0.05*(maxRe-minRe)*[-1 1];
minX=min(c.X);
maxX=max(c.X);
Xrg=[minX maxX]+0.05*(maxX-minX)*[-1 1];

gpp.lt='-';
clf;

subplot(5,1,1);
[h lh uc]=groupplot(c.g,c.X,c.Re,gpp);
ylabel('X');
xlabel('g');
axis([grg Xrg]);

% $$$ subplot(5,1,2);
% $$$ [h lh uc]=groupplot(c.Re,c.X,c.Re,gpp);
% $$$ axis([Rerg Xrg]);
% $$$ ylabel('X');

subplot(5,1,2);
[h lh uc]=groupplot(c.g,c.Re,c.Re,gpp);
axis([grg Rerg]);
ylabel('Re');
h=line(yy,xx);
cc=hsv(length(h));
for j=1:length(h)
  set(h(j),'color',cc(j,:),'linewidth',2);
end

subplot(5,1,3);
[h lh uc]=groupplot(c.Re,c.g,c.Re,gpp);
hold('on');

h=line(xx,yy);
cc=hsv(length(h));
for j=1:length(h)
  set(h(j),'color',cc(j,:),'linewidth',2);
end
plot(uRe,Sg,'s-');%plot(uRe,Sg-ff(pp,uRe),'s-');
axis([Rerg grg]);
ylabel('g');xlabel('Re');

subplot(5,1,4);
if isempty(ff1)
  plot(uRe,Sg-ff(pp,uRe),'s-',uRe,0*uRe);
else
  plot(uRe,fs(Sg-ff1(pp,uRe),Sg-ff2(pp,uRe)),'s-',uRe,0*uRe);
end

axis([Rerg -inf inf]);
ylabel('g error');xlabel('Re');

subplot(5,1,5);
plot(uRe,dXdg,'s-');
axis([Rerg -inf inf]);
xlabel('Re');ylabel('dX');
if isempty(ff1)
  fRe=ff(pp,uRe);
  for j=1:length(uRe)
    disp(sprintf('Re=%5.0f g=%5.4f  g=%5.4f g=%5.3f',uRe(j),Sg(j),fRe(j),round(fRe(j)*200)/200));
  end
else
  fRe1=ff1(pp,uRe);
  fRe2=ff2(pp,uRe);
  for j=1:length(uRe)
    disp(sprintf('Re=%5.0f g=%5.4f  g=%5.4f g=%5.3f g=%5.4f g=%5.3f',uRe(j),Sg(j),fRe1(j),round(fRe1(j)*200)/200,fRe2(j),round(fRe2(j)*200)/200));
  end
end

if ~isempty(xx)
  line(xx,yy);
  %err=Y-f(pp,X);
end
%legend(lh,uc,'location','ne');


return
nms='gc/qgc*/*';
figure(1);c=ded_gc_fit_Y_X(nms,'g','Re');
figure(2);ded_gc_fit_Y_X(c,'X','Re');

X=uRe;
Y=Sg-ff(pp,uRe);

f=@(p,x) polyval([p(2:3) 0],x/Xmax).*exp(-x/Xmax/p(1));
ff=@(p) f(p,X)-Y;
p=[0.1 1 -2];
p=[0.0500   40.0000   -3.0000];
p=lsqnonlin(ff,p);
XX=linspace(0,X(end),1e3);
plot(X,Y,'s',XX,f(p,XX));
