function ded_plot_stats(nms,typ,srt)
%ded_plot_X(nm) Plot results of X convergence
%ded_plot_X('pgc/t4/0102');i
%ded_plot_X({'gc/gc2d7n/30','gc/gc2d7n/70'});
%ded_plot_X({'gc/test/21','gc/test/22'});

if nargin<2
  typ=[];
end

if nargin<3
  srt=[];
end

if isempty(typ)
  typ={'X'};
end
if isempty(srt)
  srt={'Re'};
end

if ~iscell(typ)
  typ={typ};
end

if ~iscell(srt)
  srt={srt};
end

nms=ded_parse_nms(nms);

if isempty(nms)
  return;
end

[a nms]=ded_read_stats(nms);
p=ded_read_param(nms);

if isempty(p)
  disp(sprintf('ded_plot_X: No matching simulations'));
  disp(nms);
  return
end

ntyp=length(typ);

PIDX=zeros(size(p));
if isfield(p,'PIDX')
  for j=1:length(p)
    if ~isempty(p(j).PIDX)
      PIDX(j)=p(j).PIDX;
    end
  end
end
PIDX(~isfinite(PIDX))=0;
for j=1:length(p)
  p(j).Peb=round(p(j).Re.*p(j).Scb);
end

for k=1:length(srt)
  s(k,:)=[p.(srt{k})];
end
[S f]=sortrows(s');

p=p(f);
a=a(f);
nms={nms{f}};

uS=unique(S,'rows');


ntyp=length(typ);
for j=1:length(uS)
  f=find(all(S==uS(j,:),2));
  figure(j);clf;
  subplot(ntyp,1,1);
  for k=1:ntyp
    subplot(ntyp,1,k);
    title(nms{f(1)});
    for i=1:length(f)
      switch(typ{k})
        case('Xg')
          y=a(f(i)).X;
          x=a(f(i)).g;
          xl='g';
          yl='X';
        otherwise
          y=a(f(i)).(typ{k});
          x=a(f(i)).t;
          xl='t';
          yl=typ{k};
      end
      ny=round(length(y)/2);my=mean(y(ny:end));
      lhs=sprintf('%s %s %s:%9.6f',p(j).num,p(j).status(1:2),typ{k},my);
      plot(x,y,[min(x) max(x)],my([1 1]));
      th=text(x(end),y(end),lhs,'HorizontalAlignment','right','color',[0 0 0],'FontSize',8,'BackgroundColor',[0.7*[1 1 1] 0.2]);
      hold('on');
     end
     xlabel(xl);
     ylabel(yl);
     axis('tight');
  end
end

return;

keyboard
clf;
if ~iscell(nms)
  nms={nms};
end

[s nms]=ded_read_stats(nms);

se=[];
p=ded_read_param(nms);
if isempty(p)
  return;
end

nms={p.name};
m=length(p);

m=length(p);

sss={'Terminated','Aborted','Running','Suspended','Submitted','Unknown'};
fs=zeros(1,m);
for j=1:m
  ff=find(strcmp(p(j).status,sss));
  if isempty(ff)
    disp(sprintf('%s status is not known %s',p(j).name,p(j).status));
    ff=6;
  end
  fs(j)=ff;
end
cols=[[0.5 0.5 0.5];[1 0 0];[0 0 1];[0.7 0.7 0.7];[0 0 0];[0 1 0]];
col=cols(fs,:);

TT=repmat(NaN,2,m);
tend=repmat(NaN,1,m);
XX=repmat(NaN,2,m);
YY=repmat(NaN,2,m);

Re=[p.Re];
Peb=round([p.Re].*[p.Scb]);
ss=[];
GC=repmat(0==1,m,1);
UC=GC;
XC=UC;
mg=[p.g];
if isfield(p,'U');mU=[p.U];else;mU=repmat(NaN,size(p));end;
mX=NaN*mg;

T=ded_convergence_T(nms);
for j=1:m
  if ~isfield(s(j),'t')
    continue;
  end
  t=s(j).t;
  if ~isempty(t)
    tend(j)=t(end);
  end
  if isnan(T(j))
    T(j)=tend(j);
  end
  if isfield(s,'g')
    if length(s(j).g)==length(s(j).t)
      GC(j)=any(diff(s(j).g)~=0);
      mg(j)=mean(s(j).g(t>=T(j)));
    end
  end
  if isfield(s,'U')
    if length(s(j).U)==length(s(j).t)
      UC(j)=any(diff(s(j).U)~=0);
      mU(j)=mean(s(j).U(t>=T(j)));
    end
  end
  if isfield(s,'X')
    if length(s(j).X)==length(s(j).t)
      XC(j)=any(diff(s(j).X)~=0);
      mX(j)=mean(s(j).X(t>=T(j)));
    end
  end
end

aXC=any(XC);
aUC=any(UC);
aGC=any(GC);

nsf = aXC + aUC + aGC + (aXC&aGC) + (aXC&aUC);
if nsf==0
  return;
end
lha=zeros(5,1);
lha([aXC  aUC  aGC  (aXC&aGC)  (aXC&aUC)])=1:nsf;

ylb={'X','U','g','X','X'};
xlb={'t','t','t','g','U'};

ih=1:5;
ih(lha==0)=[];

ah=zeros(1,nsf);
for k=1:nsf
  ah(k)=subplot(nsf,1,k);
  cla;
  hold('on');
  xlabel(xlb{ih(k)},'interpreter','tex');
  ylabel(ylb{ih(k)},'interpreter','tex');
end
se=[];
for j=1:m
  
  lhs{j}=sprintf('%s %s g:%8.6f X:%9.6f',p(j).num,p(j).status(1:2),mg(j),mX(j)-PIDX(j));
  
  nsf = XC(j) + UC(j) + GC(j) + (XC(j)&GC(j)) + (XC(j)&UC(j));
  X=cell(5,1);Y=cell(5,1);
  if XC(j);       X{1}=s(j).t;Y{1}=s(j).X;end;
  if UC(j);       X{2}=s(j).t;Y{2}=s(j).U;end;
  if GC(j);       X{3}=s(j).t;Y{3}=s(j).g;end;  
  if XC(j)&GC(j); X{4}=s(j).g;Y{4}=s(j).X;end;  
  if XC(j)&UC(j); X{5}=s(j).U;Y{5}=s(j).X;end;  
 
  for k=1:5
    if isempty(X{k})
      continue;
    end
    axes(ah(lha(k)));
    nn=min(length(X{k}),length(Y{k}));
    line(X{k}(1:nn),Y{k}(1:nn),'color',col(j,:));
    th=text(X{k}(nn),Y{k}(nn),lhs{j},'HorizontalAlignment','right','color',[0 0 0],'FontSize',8,'BackgroundColor',[0.7*[1 1 1] 0.2]);
    if nrm==0 | nn<10
      gu_setalim(0,.1);
    else
      nn2=ceil(nn/2);
      axis([min(X{k}(nn2:nn)) max(X{k}(nn2:nn)) min(Y{k}(nn2:nn)) max(Y{k}(nn2:nn))]);
    end
    
    if isfield(p,'PIDX')
      if ylb{k}=='X' & p(j).PIDX>0
        line(X{k}([1 end]),p(j).PIDX*[1 1],'color',0.7*[1 1 1]);
      end
    end
    %line(t([1 end]),Zmean(j,k)*[1 1],'color',0.5*[1 1 1]);
    %line(T([1 1]),[Zmin(j,k) Zmax(j,k)],'color',0.7*[1 1 1]);
  end
  ss=sprintf('%20s, t=%7.2f, Peb=%5.0f, Re=%5.0f, T=%7.2f, g=%8.6f, U=%5.3f, X=%9.6f %s',p(j).name,tend(j),Peb(j),Re(j),T(j),mg(j),mU(j),mX(j)-PIDX(j),p(j).status);  
  if m>1
    if j==1
      se=sprintf('Peb=%5.0f, Re=%5.0f, ',Peb(j),Re(j));
    end
    %  se=[se sprintf('%s g=%8.6f X=%9.6f, ',p(j).num,mg(j),mX(j)-PIDX(j))];
  end
  disp(ss);
end


if 0
  for j=1:k
    legend(h1(:,k),lhs,'location','best')
  end
end

axes(ah(1));

hhh=[];
for k=1:length(sss)
  hhh(k)=line(NaN,NaN,'color',cols(k,:));
end
if m>1
  legend(hhh,sss,'location','best');
end

if ~isempty(se)
  title(se,'fontsize',8);
else
  title(ss,'fontsize',8);
end

% $$$ tmin=min(tmin,1);tmax=max(tmax,1);
% $$$ Zmin=min(Zmin,1);tmax=max(Zmax,1);
% $$$ trg=[tmin tmax];
% $$$ Zrg=[Zmin' Zmax'];
  
return;
Xrg=Xrg + max(1e-5,diff(Xrg,1,2)/20)*[-1 1];
Yrg=Yrg + max(1e-5,diff(Yrg,1,2)/20)*[-1 1];

DD=[ded_dedalus_data_dir '/'];
if nf
  for j=1:length(nms)
    figure(j);
    subplot(nsf,1,1);
    title(sprintf('%s g=%8.6f X=%9.6f Re=%5.0f',nms{j}(length(DD)+1:end),mg(j),mX(j)-PIDX(j),Re(j)));
    xlabel('t','interpreter','tex');
    ylabel('X','interpreter','tex');
    axis([trg(j,:),Xrg(j,:)]);
    line(TT(:,j),Xrg(j,:),'color',0.7*[1 1 1],'color',[0 0 0]);
    line(trg(j,:),XX(:,j),'color',0.7*[1 1 1],'color',[0 0 0]);
    subplot(nsf,1,2);
    xlabel('t','interpreter','tex');
    ylabel(ylbl,'interpreter','tex');
    axis([trg(j,:),Yrg(j,:)]);
    line(TT(:,j),Yrg(j,:),'color',0.7*[1 1 1],'color',[0 0 0]);
    line(trg(j,:),YY(:,j),'color',0.7*[1 1 1],'color',[0 0 0]);
    subplot(nsf,1,3);
    ylabel('X','interpreter','tex');
    xlabel(ylbl,'interpreter','tex');
    axis([Yrg(j,:),Xrg(j,:)]);
    line(YY(:,j),Xrg(j,:),'color',0.7*[1 1 1],'color',[0 0 0]);
    line(Yrg(j,:),XX(:,j),'color',0.7*[1 1 1],'color',[0 0 0]);
  end
else
  subplot(nsf,1,1);
  xlabel('t','interpreter','tex');
  ylabel('X','interpreter','tex');
  line(TT,Xrg,'color',0.7*[1 1 1],'linestyle','--');
  line(trg,XX,'color',0.7*[1 1 1],'linestyle','--');
  axis([trg Xrg]);
  subplot(nsf,1,2);
  xlabel('t','interpreter','tex');
  ylabel(ylbl,'interpreter','tex');
  line(TT,Yrg,'color',0.7*[1 1 1],'linestyle','--');
  line(trg,YY,'color',0.7*[1 1 1],'linestyle','--');
  axis([trg Yrg]);
  subplot(nsf,1,3);
  ylabel('X','interpreter','tex');
  xlabel(ylbl,'interpreter','tex');
  line(Yrg,XX,'color',0.7*[1 1 1],'linestyle','--');
  line(YY,Xrg,'color',0.7*[1 1 1],'linestyle','--');
  axis([Yrg Xrg]);
  %legend(h1,nms,'location','northeastoutside');
  %legend(h2,nms,'location','northeastoutside');
  %legend(hnsf,nms,'location','northeastoutside');
end

