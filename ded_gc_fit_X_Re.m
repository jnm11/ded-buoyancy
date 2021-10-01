function g=ded_gc_fit_X_Re(nm,tol,fr)
%ded_gc_fit_U_Re('9[0-2]*');
if nargin<2
  tol=[];
end
if nargin<3
  fr=[];
end
if isempty(tol)
  tol=0.01;
end
if isempty(fr)
  fr=0;
end


a=ded_all_stats('gc',nm);

A=[[a.h];[a.H];[a.U];[a.forcing]]';

AA=unique(A,'rows');
if size(AA,1)>1
  disp('ded_gc_fit_X_U: non matching simulations');
  f=find(all(A(1,:)==AA,2));
  a=a(f);
end

x=[a.Re];
[x f]=sort(x);
a=a(f);

y=[a.fXm]; 
f=isfinite(y);
a=a(f);

x=[a.Re];
y=[a.fXm]; 
maxx=max(x);
minx=min(x);
maxy=max(y);
miny=min(y);


f = @(p,x) p(1)-p(2)*(x/maxx).^(-p(3));
ff= @(p) f(p,x)-y;
p=[0 0 1];
%p(1)-p(2)=maxy
%p(1)-p(2)*minx/maxx=miny
p(2)=(maxy-miny)/(1+maxx/minx);
p(1)=maxy+p(2);

p=lsqnonlin(ff,p,[0 0 0]);

xx=linspace(minx,maxx,1e3);
clf;
plot(x,y,'-s',xx,f(p,xx));
xlabel('U');
ylabel('X');
title(nm);
e=abs(ff(p)./y);

disp(sprintf('%8s %7s %7s %7s %8s %7s %12s','name','Re','U','X','t','err','status'));
for j=1:length(e)
  disp(sprintf('%8s %7.0f %7.4f %7.3f %8.2f %7.4f %12s',a(j).name,a(j).Re,a(j).U,a(j).fXm,a(j).t,e(j),a(j).status));
  ft(j)=  strcmp(a(j).status,'Terminated') & (e(j)>tol   |  ~isfinite(a(j).T));  % run longer
  fr(j)= ~strcmp(a(j).status,'Terminated') & (e(j)<tol/2 &   isfinite(a(j).T));  % Terminate
end

f=find(ft);
for k=1:length(f)
  j=f(k);
  t=round(a(j).t+max(a(j).t*0.5,50));
  cmd=sprintf('ded_setparam.py ~/gc/%s/param.h5 T %.1f',a(j).name,t);
  disp(cmd); if fr>0; unix(cmd);end
  cmd=sprintf('echo Aborted > ~/gc/%s/status',a(j).name);
  disp(cmd); if fr>0; unix(cmd);end
end

f=find(fr);
for k=1:length(f)
  j=f(k);
  cmd=sprintf('ded_setparam.py ~/gc/%s/param.h5 T %.1f',a(j).name,floor(a(j).t));
  disp(cmd); if fr>0; unix(cmd);end
  cmd=sprintf('echo Terminated > ~/gc/%s/status',a(j).name);
  disp(cmd); if fr>0; unix(cmd);end
  cmd=sprintf('echo Terminated > ~/gc/%s/abort',a(j).name);
  disp(cmd); if fr>0; unix(cmd);end
end

g= @(x) f(p,x);
