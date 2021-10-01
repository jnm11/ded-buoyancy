function a=ded_gc_find_head(x,h,xmin,display)

if nargin <4
  display=[];
end

if isempty(display)
  display=nargout==0;
end

if ischar(x)
  nm=x;
   pr=load([nm '/profile.mat']);
   x=pr.a.x;
   h=pr.a.b;
end


rg = find(x>xmin);  % Ignore forcing region
x=x(rg);x=x(:);
h=h(rg);h=h(:);


nh=length(h);
while sum(abs(diff(sign(diff(h))))) > 50
  h=(h([1 1:end-1])+2*h(1:end)+h([2:end end]))/4;
end

  

hmax=max(h);

H=1.05*hmax;

f1=max(find(h>hmax/5));
f2=min(find(h<hmax/1e4 & x>x(f1)));
f2=2*f2-f1;

xx=x(f1:f2);
mx=mean(xx);
xx=xx-mx;

ff = @(p,x) p(1).*( exp(-p(2)^2.*(p(3)-x).^2)+ p(2).*(p(3)-x).*sqrt(pi).*(1+erf(p(2).*(p(3)-x))))/(2.*p(2)^2);

fff = @(p) ff(p,xx)-h(f1:f2);
opt=optimset('display','none');
p=[1 500 0 ];
p(1)=h(f1(1))/ff(p,xx(1));
p=lsqnonlin(fff,p,[],[],opt);
%plot(xx,h(f1:f2),xx,ff(p,xx));

Xf=p(3)+mx;
hf=ff(p,Xf);

if false
  plot(xx,h(f1:f2),xx,ff(p,xx));
  
  h1=0;
  h2=hmax/10;
  hf=hmax;
  while 1
    oh=hf;
    hf=(h1+h2)/2;
    if abs(oh-hf)<1e-4 & n==1 & h1>0
      break;
    end
    dh=h-hf;
    n = sum(dh(1:end-1).*dh(2:end)<=0);
    if n<=1
      h2=hf;
    else
    h1=hf;
    end
  end
  hf=hf*2;
  keyboard
  f1=max(find(h>hf));
  f2=f1+1;
  %(x(f2)-x)*h(f1)+(x-x(f1))*h(f2)=hf*(x(f2)-x(f1));
  %x*(h(f2)-h(f1)) = hf*(x(f2)-x(f1)) + x(f1)*h(f2)-x(f2)*h(f1);
  Xf = (hf*(x(f2)-x(f1)) + x(f1)*h(f2)-x(f2)*h(f1))/(h(f2)-h(f1));
  
end

f0=(0:2)+max(find(h(2:end-1)>h(1:end-2) & h(2:end-1)>h(3:end)&x(2:end-1)<Xf&h(2:end-1)>0.3*hmax));
[Xh hh]=qmax(x(f0),h(f0));
f=(0:2)+max(find(h(2:end-1)<h(1:end-2) & h(2:end-1)<h(3:end)&x(2:end-1)<Xh-hh & h(2:end-1)<hh*0.9));
if isempty(f)
  a.Xf=NaN;
  a.hf=NaN;
  a.Xh=NaN;
  a.hh=NaN;
  a.Xt=NaN;
  a.ht=NaN;
  a.alpha=NaN;
  a.X0=NaN;
  a.h0=NaN;
  a.alpha1=NaN;
  a.alpha2=NaN;
  return;
end
[Xt ht]=qmax(x(f),h(f));

xmin=x(1);
f=find(x>=xmin & x<=Xt);
pp=polyfit(x(f),h(f),1);



  
a.Xf=Xf;
a.hf=hf;
a.Xh=Xh;
a.hh=hh;
a.Xt=Xt;
a.ht=ht;
a.alpha=-pp(1); % entrainment coefficient


p(1)=0.01;
p(2)=1;
p(3)=h(1);
p(4)=a.alpha;
p(5)=0;
p(6)=Xf;
p(7)=Xt;
p(8)=1;
X1=Xf-Xh;
X2=Xf-Xt;
p(9:12)=[0 0 -(X1*ht-X2*hh) -(-X1^2*ht+X2^2*hh)]/(X2*X1*(X1-X2));


x1=x(1);
f =@(p) gc_h_fun(p,x,x1)-h;

dX=max(0.1,Xf-Xh);

opt=optimset('display','none');
lb=[-inf   1   0   0 -inf Xf-dX/2 Xh-10*dX  1/dX -inf -inf -inf -inf];
ub=[ inf inf inf inf  inf Xf+dX/2 Xh-dX/2 10/dX  inf  inf  inf  inf];


fb=false;
while(true)
  p=lsqnonlin(f,p,lb,ub,opt);
  pt=p(1:5);
  Xf=p(6);
  ps=p(7:8); % switch parameters
  ph=p(9:end);
  Xs=ps(1);
  ws=ps(2);
  r=roots(poly_diff([ph 0]));
  r=r(imag(r)==0);
  switch(length(r))
    case(0)
      p(end+1)  = 0;
      lb(end+1) = -inf;
      ub(end+1) =  inf;
      fb=true;
    case(1)
      break;
    otherwise
      if fb; break; end;
      p(end)=[];
      lb(end)=[];
      ub(end)=[];
  end
end

% $$$ r=roots(ph);
% $$$ while any(r>=0 & r>=Xf-Xs)
% $$$   ph(1)=[];
% $$$   r=roots(ph);
% $$$ end


[f f1 f2 f3 Xt ht]=gc_h_fun(p,x,x1);




ht=min(ht,max(h(x>=xmin&x<=Xt)));

Xh=fminbnd(@(x) -gc_h_fun(p,x,x1), Xs,Xf);
%Xt=fminbnd(@(x) +gc_h_fun(p,x,x1), 2*Xs-Xh,Xh);
%ht=gc_h_fun(p,a.Xt,x1);

a.X0=x(1);
a.Xf=Xf;
a.Xs=Xs;
a.Xt=Xt;
a.Xh=Xh;
a.ht=ht;
a.hh=gc_h_fun(p,a.Xh,x1);
a.hf=gc_h_fun(p,a.Xf,x1);
a.h0=gc_h_fun(p,a.X0,x1);
a.hs=gc_h_fun(p,a.Xs,x1);

a.alpha=(a.h0-a.ht)/(a.Xt-a.X0); % Average entrainment gradient

if display
  figure;clf;hold('on');
  plot(x,h,'-')
  plot(a.Xf([1 1]),[0 H],'k',a.Xf,a.hf,'ks');
  plot(a.Xh([1 1]),[0 H],'r',a.Xh,a.hh,'rs');
  plot(a.Xt([1 1]),[0 H],'b',a.Xt,a.ht,'bs');
  plot(a.Xs([1 1]),[0 H],'y',a.Xs,a.hs,'ys');
  plot(x([1 end]),a.h0([1 1]));
  plot(x,f);%,x,(1-f3).*f1,x,f3.*f2);
  axis([x(1) x(end) 0 max(h)*1.05]);

  ft=find(x>=a.X0&x<=a.Xt);
  fh=find(x>=a.Xt&x<=a.Xf);
  plot(x(ft),f1(ft),x(fh),f2(fh))
end

% $$$ if a.alpha<5e-3;
% $$$   keyboard
% $$$ end


% The switch isnt the best approach
% Better to have 3 polynomial sections
% 1 tail
% 2 behind the head with a minimum
% 3 the head with a maximum and a zero
function [f f1 f2 f3 Xt ht]=gc_h_fun(p,x,x1)


pt=p(1:5);
Xf=p(6);
ps=p(7:8); % switch parameters
ph=p(9:end);

Xs=ps(1);
ws=ps(2);



q1=[ph 0];




f1 = pt(1)*exp(-pt(2)*(x-x1))  + pt(3) - pt(4)*(x-x1)-pt(5)*(x-x1).^2;
f2 = polyval(q1,max(Xf-x,0)); 
f3 = (1+erf((x-Xs)*ws))/2;
f  = (1-f3).*f1 +f3.*f2;

if nargout>4
  x2=x1-Xf;
  q3=[-pt(5) -2*pt(5)*x2+pt(4) -pt(5)*x2^2+pt(4)*x2+pt(3)];
  q2=q1;
  q2(end-2:end)=q2(end-2:end)-q3;
  r  =  roots(q2);
  r  = Xf-r(imag(r)==0);
  if isempty(r)
    r=NaN;
    Xt=Xs;
  else
    Xt = r(findmin(abs(r-Xs)));
  end
  Xt=min(Xt,Xs);
  ht=polyval(q1,Xf-Xt);
  ht=pt(1)*exp(-pt(2)*(Xt-x1))  + pt(3) - pt(4)*(Xt-x1)-pt(5)*(Xt-x1).^2;
  %  plot(x,f1,x,max(0,f2),x,max(0,polyval(q1,Xf-x)),'--',x,max(0,polyval(q3,Xf-x)),'--',Xt,ht,'s')
  %keyboard
end

return;
dc;fns=cellstr_ls('gc/f6/*/*/*/2304/xxxxh');
dc;fns=cellstr_ls('gc/f6/mg/*/*/2304/[0-9]*');
for j=1:length(fns);
  nm=fns{j};
  fn=[ded_dedalus_data_dir '/results/' nm '/profile-ray.mat'];
  load(fn);
  xmin=1;
  fr=ded_gc_find_head(a.x,a.bm0,xmin,true);
end
