function [c f pp err]=ded_gc_fit_X_Re_g(nms,XX,typ,p0)


nmsin=nms;
if nargin<3
  typ=[];
end
if nargin<4
  p0=[];
end
f=[];

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
  c.W(j)  = 1+s.t(end)-T;
  c.Re(j) = p.Re;
  c.g(j)  = p.g;
  c.U(j)  = p.U;
  c.X(j)  = mean(s.X(f));
  c.nm{j} = nm;
end



maxRe=max(c.Re);
minRe=min(c.Re);
maxg=max(c.g);
ming=min(c.g);
Re=c.Re;
g=c.g;
X=c.X; % Take of target point;

N=4;
fg   = @(p,x)ming+polyval(p(1:N),x/maxRe)./max(eps,polyval([1 p(N+1:2*N-1) 0],x/maxRe));
pfg  = zeros(1,2*N-1);pfg(1)=1;
fdg  = @(p,x) max(0,polyval(p(1:2),x/maxRe)./max(eps,polyval([p(3) 1],x/maxRe)));
pdg  = [0 1 0];
fX   = @(p,g,Re) XX+(g-fg(p(1:5),Re)).*fdg(p(6:10),Re);

ReRe=linspace(minRe,maxRe,1e3);
ff = @(p) fg(p,Re)-g;
pfg=lsqnonlin(ff,pfg,[],[],opt);
figure(1);clf;
subplot(2,1,1);
[h lh uc]=groupplot(c.Re,c.g,c.Re,gpp);
hold('on');
line(ReRe,fg(pfg,ReRe),'color',[0.7 0.7 0.7],'linewidth',2);
subplot(2,1,2);
[h lh uc]=groupplot(c.Re,c.g-fg(pfg,c.Re),c.Re,gpp);


ff = @(p) fX(p,g,Re)-X;
opt=optimoptions('lsqnonlin','display','none','tolx',1e-8,'tolfun',1e-8);  
p0=[pfg pdg];
pp=lsqnonlin(ff,p0,[],[],opt);
pfg=pp(1:5);
pdg=pp(6:8);

line(ReRe,fg(pfg,ReRe),'color',[1 0 0],'linewidth',2);


figure;clf;
plot(ReRe,fdg(pdg,ReRe));


subplot(3,1,1);
[h lh uc]=groupplot(c.Re,c.X,c.g,gpp);
axis('tight')
line(ReRe,fX(pp,g,ReRe),'color',[0.7 0.7 0.7],'linewidth',4);
ylabel('X');

subplot(3,1,2);
plot(ReRe,100*fdg(pdg,ReRe),'color',[0.7 0.7 0.7],'linewidth',4);
ylabel('dX/dg');

subplot(3,1,3);
[h lh uc]=groupplot(c.Re,c.g,c.Reg,gpp);
axis('tight')
subplot(3,1,3);
line(ReRe,ming+(maxg-ming)*fg(pfg,ReRe),'color',[0.7 0.7 0.7],'linewidth',2);
ylabel('g');
xlabel('Re');
