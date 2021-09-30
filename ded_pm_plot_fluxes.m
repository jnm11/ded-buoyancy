function [a f]=ded_pm_plot_fluxes(nm,trg)
%ded_pm_plot_fluxes('pm/spm/004');
%ded_pm_plot_fluxes('pm/f7/d/03',[70 inf]);
%ded_pm_plot_fluxes('pm/f7/d/04',[100 inf]);
%rsync -vap hamilton:pm/f7/d/ ~/pm/f7/d --exclude b --exclude final --exclude a 
if nargin<2
  trg=[];
end
if isempty(trg)
  trg=[-inf inf];
end
p=ded_read_param(nm);
c=ded_coord(nm);
b=ded_read_javrg(nm,'ayz',trg);%,[],[p.t/2 inf]);
a=ded_read_g(nm,'yz',[],trg);%,[],[p.t/2 inf]);
if isempty(a)
  return;
end

if ~isempty(b)
  na=2*(isfield(b,'u')+isfield(b,'b')+isfield(b,'s'));
end

p.U=p.V;

A=pi*p.Radius^2;
U=p.U;
if strcmp(p.Ty,'SinCos');A=A/2;end
if strcmp(p.Tz,'SinCos');A=A/2;end


  
if 1
  bf=ded_read_g(nm,'force');%,[],[p.t/2 inf]);
  bf=[];
  if isempty(bf)
    switch(p.Forcing)
      case 4
        xu=[p.x1 p.x2];
        if p.x0==p.x1
          f.u  = interp1([-1 p.x1-p.Wx p.x1 p.x2 p.x3 p.L],[0 0 1 1 0 0],c.Jx,'pchip');
        else
          f.u  = interp1([-1 p.x0 p.x1 p.x2 p.x3 p.L],[0 0 1 1 0 0],c.Jx,'pchip');
        end
        f.u(c.Jx<p.x1 | c.Jx>p.x2)=NaN;
        f.b  = interp1([-1 p.x4 p.x5 p.x6 p.x7 p.x7+p.Wx p.L],[0 0 0 1 1 0 0],c.Jx,'pchip');
      case 6
        xu=[p.x1 p.x2];
        f.u  = interp1([-1 p.x0 p.x1 p.x2 p.x3 p.L],[0 0 1 1 0 0],c.Jx,'pchip');
        f.wb = interp1([-1 p.x4 p.x4+p.Wx p.x7-p.Wx p.x7 p.L],[0 0 1 1 0 0],c.Jx,'pchip');
        f.b  = interp1([-1 p.x5 p.x6      p.x7 p.x7+p.x6-p.x5  p.L],[0 0 1 1 0 0],c.Jx,'pchip');
        f.b(f.wb==0)=NaN;
      case 7
        fx=ded_hdf(nm,'force/fx');
        if isfield(fx,'wul') & isfield(fx,'wur')  
          f.u = zeros(c.NAx);
          f.u(1:length(fx.wul))         = p.U*fx.wul;
          f.u(end+1-length(fx.wur):end) = p.U*fx.wur;
          urg=[length(fx.wul)+1 c.NAx-length(fx.wur)];
          brg=urg;          
          xrg=urg;          
        end
        if isfield(fx,'Bwxx')
          f.b = p.B*fx.Bwxx;
        end
        if isfield(fx,'wux5') & isfield(fx,'wux2')
          ff=fx.wux5+fx.wux2;
          ff=c.Ax(ff<1e-3*max(ff));
          urg=[min(ff) max(ff)];
          ff=c.Ax(fx.wux5<1e-3*max(fx.wux5));
          brg=[min(ff) max(ff)];
          xrg=[max(brg(1),urg(1)) min(brg(2),urg(2))];
          f.u = midpoint(p.U*cumsum([0;fx.wux1])*c.dAx);
          f.b = p.B*(fx.wux5>1e-3*max(fx.wux5));
        end
        if isfield(p,'S')
          f.s = p.S*f.b/p.B;
        else
          f.s=NaN*f.b;
        end
    end
    f.bu = f.b.*f.u;
    f.su = f.s.*f.u;
    f.uu = f.u.*f.u;
    f.bb = f.b.*f.b;
    f.ss = f.s.*f.s;
  else
    dA=(bf.y(2)-bf.y(1))*(bf.z(2)-bf.z(1))/A;
    f.b  = dA*squeeze(sum(sum(bf.fb,1),2));
    f.b  = dA*squeeze(sum(sum(bf.fb,1),2));
    if isfield(b,'fu')
      f.u  = dA*squeeze(sum(sum(bf.fu,1),2));
      f.v  = dA*squeeze(sum(sum(bf.fv,1),2));
      f.w  = dA*squeeze(sum(sum(bf.fw,1),2));
      f.bu = dA*squeeze(sum(sum(bf.fu.*bf.fb,1),2));
      f.uu = dA*squeeze(sum(sum(bf.fu.*bf.fu,1),2));
     else
      f.u = p.U*interp1([0 p.x4 p.x4+1e-5 p.L],[0 1 0 0],c.Jx,'pchip');
      f.v = 0*c.Jx;
      f.w = 0*c.Jx;
      if isfield(f,'b')
        f.bu = f.b.*f.u;
        f.bb = f.b.*f.b;
      else
        f.bu=NaN*f.u;
        f.bb=NaN*f.u;
      end
      if isfield(f,'s')
        f.su = f.s.*f.u;
        f.ss = f.s.*f.s;
      else
        f.su=NaN*f.u;
        f.ss=NaN*f.u;
      end
      f.uu = f.u.*f.u;
    end
  end
  f.x=c.Jx;
  f.t=a.t;
end

f.s=f.b;
f.v=0*f.u;
f.w=0*f.u;
f.ss=f.s.*f.s;
f.su=f.s.*f.u;

%if 1
%  s=ded_read_state([ded_dedalus_data_dir '/' nm '/final/final_s1.hdf5']);
% s=ded_read_state([ded_dedalus_data_dir '/' nm '/checkpoint/checkpoint_s106.hdf5']);

%xrg=[2   p.L-3];
fx=find(c.Jx>=xrg(1) & c.Jx <=xrg(2));
fx=1:length(c.Jx);
x=c.Jx(fx);
f.u=f.u(fx);
f.v=f.v(fx);
f.w=f.w(fx);
f.b=f.b(fx);
f.s=f.s(fx);
f.uu=f.uu(fx);
f.bb=f.bb(fx);
f.ss=f.ss(fx);
f.bu=f.bu(fx);
f.su=f.su(fx);


if 0
  u=mean(a.u(fx,ft),2);
  %[pu nu]=fit_pow(x,u,'3',[],1);
  sz=[prod(size(a.u(fx,ft))) 1];
  ff= @(p) u./max(1,x-p(1)).^(5/3)/A-p(2);
  pu=[X0 p.U];
  opt=optimset('display','none');
  pu=lsqnonlin(ff,pu,[-inf 0],[inf inf],opt);
  plot(x,u/A,x,pu(2)*max(1,x-pu(1)).^(5/3));
  X=max(1,x-pu(1));
end




arg=[0 p.L 0 inf];
%trg=[min(a.t(ft)) max(a.t(ft))];

for k=1:3
  switch(k)
    case 1
      ft=find(a.t>=trg(1) & a.t<=trg(2));
      aa=a;
      arg=[0 p.L 0 2];
      yticks=[0:2];
    case 2
      ft=find(a.t>=trg(1) & a.t<=trg(2));
      aa=a;
      arg=[0 p.L 0 inf];
      yticks=[];
   case 3
      ft=1;
      aa=b;
      arg=[0 p.L 0 inf];
      yticks=[];
  end
  if ~isempty(aa)
    au  =  aa.u(fx,ft)./(A*U);
    auu = aa.uu(fx,ft)./(A*U*U);
    maxu=max(max(max(au)),max(f.u(:)));
    maxuu=max(max(max(auu)),max(f.uu(:)));
    if isfield(aa,'b');
      ab  =  aa.b(fx,ft)./(A);
      abb = aa.bb(fx,ft)./(A);
      abu = aa.bu(fx,ft)./(A*U);
      maxb=max(max(max(ab)),max(f.b(:)));
      maxbu=max(max(max(abu)),max(f.bu(:)));
    end
    if isfield(aa,'s')
      as  =  aa.s(fx,ft)./(A);
      ass = aa.ss(fx,ft)./(A);
      asu = aa.su(fx,ft)./(A*U);
      maxs=max(max(max(as)),max(f.s(:)));
      maxsu=max(max(max(asu)),max(f.su(:)));
    end
  end
  

  figure(k);clf;
  ah=jsubplot([1 na],[0.1 0.1],[0.04 0.04],[0.01 0.1]);
  j=0;
  if isfield(a,'u')
    j=j+1;
    axes(ah(j));
    plot(x,au,x,f.u);
    line([1;1]*urg,repmat([0;maxu],1,2),'color',0.7*[1 1 1] );
    line(arg(1:2),[1 1],'color',0.7*[1 1 1] );
    ylabel('u/A');
    title(sprintf('%s t=[%7.2f %7.2f]',nm,trg)); 
    if ~isempty(yticks);set(gca,'ytick',yticks);end;
    axis(arg);
  end
  
  if isfield(b,'b')
    j=j+1;
    axes(ah(j));
    plot(x,ab);
    line([1;1]*brg,repmat([0;maxb],1,2),'color',0.7*[1 1 1] );
    line(arg(1:2),[1 1],'color',0.7*[1 1 1] );
    ylabel('b/A');
    axis(arg);
    if ~isempty(yticks);set(gca,'ytick',yticks);end;
  end
  
  if isfield(b,'s')
    j=j+1;
    axes(ah(j));
    plot(x,as);
    line([1;1]*brg,repmat([0;maxs],1,2),'color',0.7*[1 1 1] );
    line(arg(1:2),[1 1],'color',0.7*[1 1 1] );
    ylabel('s/A');
    axis(arg);
    if ~isempty(yticks);set(gca,'ytick',yticks);end;
  end
  
  if isfield(b,'b') & isfield(b,'u')
    j=j+1;
    axes(ah(j));
    plot(x,abu)
    line([1;1]*brg,repmat([0;maxbu],1,2),'color',0.7*[1 1 1] );
    line([1;1]*urg,repmat([0;maxbu],1,2),'color',0.7*[1 1 1] );
    line(arg(1:2),[1 1],'color',0.7*[1 1 1] );
    ylabel('bu/A');
    axis(arg);
    if ~isempty(yticks);set(gca,'ytick',yticks);end;
  end
  
  if isfield(b,'s') & isfield(b,'u')
    j=j+1;
    axes(ah(j));
    plot(x,asu)
    line([1;1]*brg,repmat([0;maxsu],1,2),'color',0.7*[1 1 1] );
    line([1;1]*urg,repmat([0;maxsu],1,2),'color',0.7*[1 1 1] );
    line(arg(1:2),[1 1],'color',0.7*[1 1 1] );
    ylabel('bs/A');
    axis(arg);
    if ~isempty(yticks);set(gca,'ytick',yticks);end;
  end
  
  if isfield(b,'u')
    j=j+1;
    axes(ah(j));
    plot(x,auu);
    line([1;1]*urg,repmat([0;maxuu],1,2),'color',0.7*[1 1 1] );
    line(arg(1:2),[1 1],'color',0.7*[1 1 1] );
    ylabel('uu/A');
    axis(arg);
    if ~isempty(yticks);set(gca,'ytick',yticks);end;
  end
  set(ah(1:end-1),'xticklabels',[]);
end

if ~isempty(b)
  ftolb=max(find(ab/max(ab)>1e-2)):length(ab);
  wu=au./sqrt(auu);
  wu(ftolb)=NaN;
  
  if isfield(b,'b')
    wb=ab./sqrt(abb);
    ub=abu./ab;
    ub(ftolb)=NaN;
    wb(ftolb)=NaN;
  end
  if isfield(b,'s')
    ws=as./sqrt(ass);
    us=asu./as;
    ws(ftolb)=NaN;
    us(ftolb)=NaN;
  end
  figure;clf
  subplot(2,1,1);
  h=plot(x,wu);axis(arg);ylabel('w');
  l={'wu'};
  hold('on');
  if isfield(b,'b');h(end+1)=plot(x,wb);l{end+1}='wb';end
  if isfield(b,'s');h(end+1)=plot(x,ws);l{end+1}='ws';end
  legend(h,l);
  subplot(2,1,2);
  l={};h=[];
  cla;hold('on');
  if isfield(b,'b');h(end+1)=plot(x,ub);l{end+1}='ub';end
  if isfield(b,'s');h(end+1)=plot(x,us);l{end+1}='us';end
  axis(arg);ylabel('u');
  legend(h,l);
end


if isfield(a,'b');Tb=sum(a.b,1)*c.dx;end;
if isfield(a,'s');Ts=sum(a.s,1)*c.dx;end;
if isfield(a,'u');Tu=max(a.u,[],1);  end;
Tt=a.t;
s=ded_read_stats(nm);
if ~isempty(s)
  n=0;
  n=1+isfield(s,'divx')+isfield(s,'mun');
  j=0;
  figure;clf
  if isfield(s,'divx')
    j=j+1;
    subplot(n,1,j);
    plot(s.t,s.divx,Tt,Tu);
    xlabel('t');
    ylabel('divx');
  end
  if isfield(s,'mun')
    j=j+1;
    subplot(n,1,j);
    plot(s.t,s.mun);
    xlabel('t');
    ylabel('mun');
  end
  j=j+1;
  subplot(n,1,j);
  cla;hold('on');
  if isfield(a,'b');plot(a.t,Tb);end;
  if isfield(a,'s');plot(a.t,Ts);end;
  xlabel('t');
  ylabel('Tb & Ts');
end

return;


nm='pm/f7/10';
p=ded_read_param(nm);
a=ded_read_g(nm,'yz');
f=ded_read_hdf([ded_dedalus_data_dir '/' nm '/final/state-00007.hdf5']);
IA=p.Radius^2*pi;

dA=(f.y(2)-f.y(1))*(f.z(2)-f.z(1));
g.b=squeeze(sum(sum(f.b)))*dA;
g.u=squeeze(sum(sum(f.u)))*dA;
g.bb=squeeze(sum(sum(f.b.*f.b)))*dA;
g.uu=squeeze(sum(sum(f.u.*f.u)))*dA;
g.bu=squeeze(sum(sum(f.b.*f.u)))*dA;

plot(f.x,g.u/IA,f.x,a.u(:,end)/IA);
axis([0 5 0 2]);

plot(f.x,g.b/IA,f.x,a.b(:,end)/IA);
axis([0 5 0 2]);


plot(f.x,g.u/IA);
plot(x,a.b(:,end)/IA);
axis([0 5 0 2]);
