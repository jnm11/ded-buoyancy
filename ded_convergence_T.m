function [T b a]=ded_convergence_T(nm,p,display)

if nargin<2
  p=[];
end
if nargin<3
  display=[];
end
if isempty(display)
  display=nargout==0;
end
if iscell(nm)
  for j=1:length(nm)
    T(j)=ded_convergence_T(nm{j},p,display);
  end
  return;
end
T=NaN;
b=[];
a=[];
if isfield(p,'T')
  T=p.T;
  return;
end
if isfield(p,'trg')
  T=p.trg(1);
  return;
end
ip.fd='/tmp';
ip.fnt=10;
ip.sz=[4 2];
ip.nc=5;
ip.typ={'b','bx','bz','bxx','bxz','bzz'};
ip.typ={'X','b','g','divz','U','divx'};
ip.t=NaN;
ip.mint=0;
ip.tol=1e-4;
p=combine_struct(ip,p);
if ~iscell(p.typ)
  p.typ={p.typ};
end
%b=ded_read_mom(nm);
if 0
  tnm={'xyz','axyz','yz','ayz','y','ay'};
  for j=1:length(tnm)
    b=ded_read_g(nm,tnm{j});
    if ~isempty(b) ; break; end
  end
end
b=ded_read_stats(nm);
if isempty(b)
  disp(sprintf('ded_convergence_T: No stats data for %s',nm));
  dt=NaN;
  a=[];
  T=NaN;
  return;
end

if 0
  if isempty(b)
    T=NaN;
    dt=0;
    bb.mean=NaN;
    bb.std=NaN;
    bb.min=NaN;
    bb.max=NaN;
    for j=1:length(p.typ)
    a.(p.typ{j})=bb;
    end
    a.MB=NaN;
    a.SB=NaN;
    a=orderfields(a);
    return;
  end
end

if isempty(b)
  return;
end
if ~isfield(b,'t');
  return;
end

switch length(b.t)
  case 0
    T=NaN;
    dt=0;
    return;
  case 1
    T=b.t(1);
    dt=0;
    return;
end
p.typ=intersect(p.typ,fieldnames(b));

ntyp=length(p.typ);
if length(p.tol)==1
  p.tol=repmat(p.tol,ntyp,1);
end
for j=1:ntyp
  f1(j)=find_stationary_crossings(b.(p.typ{j}),p.nc,0);
  f2(j)=find_tol(b.(p.typ{j}),p.tol(j),false,ntyp,j);

end
ff=max(min(f1,f2));
dt=max(p.mint,b.t(end)-b.t(ff));
T=b.t(end)-dt;

if display
  figure;
  clf;
  hold('on');
  for j=1:length(p.typ)
    x=b.(p.typ{j});
    x=x-min(x);
    x=x/max(x);
    plot(b.t,x,b.t(f1(j))*[1 1],[0 1],b.t(f2(j))*[1 1],[0 1],'color',0.5*[1 1 1]);
  end
  plot([T T],[0 1],'color',[0 0 0],'linewidth',2);
  axis('tight');
  xlabel('t');
  title(sprintf('Convergence time for %s is %6.1f/%6.1f',nm,T,b.t(end)));
end

if nargout >3
  f=find(b.t>=T | b.t>=b.t(end)*.99);
  for j=1:length(p.typ)
    a.(p.typ{j}).mean=mean(b.(p.typ{j})(f));
    a.(p.typ{j}).std=std(b.(p.typ{j})(f));
    a.(p.typ{j}).min=min(b.(p.typ{j})(f));
    a.(p.typ{j}).max=max(b.(p.typ{j})(f));
  end
  %b    = [a.b];   b    =  [b.mean];
  %bx   = [a.bx];  bx   = [bx.mean];
  %bxx  = [a.bxx]; bxx = [bxx.mean];
  %a.mb =  bx/b;
  %a.sb = sqrt(max(0,bxx/b-a.mb^2));
  a=orderfields(a);
end

return


% Find region where std deviation is less than tol
function f=find_tol(x,tol,display,ntyp,j);
n=length(x);
x=flip(x(:)');
mc=cumsum(x)./(1:n);       % Cumulative mean
ms=cumsum(x.^2)./(1:n);    % Cumulative sum of squares
s=sqrt(max(0,ms-mc.^2));   % Cumulative std
f=n+2-max(find(s<=tol));   % s(1) is always  0 return n+1 for this

if display
  maxs=1.05*max(max(s),tol);
  subplot(ntyp,1,j);cla;
  plot(1:n,s,[f f],[0 maxs],[1 n],[tol tol]);
  axis([0 n+1 0 maxs]);
end

