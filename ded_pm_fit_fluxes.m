function ded_pm_fit_fluxes(nm,trg,xrg,fitg) 
if nargin==0;
  plot_all_e;
  return;
end

if nargin<4
  fitg=[];
end
if nargin<3
  xrg=[];
end
if nargin<2
  trg=[];
end
if isempty(trg)
  trg=[-inf inf];
end
if isempty(fitg)
  fitg=1;
end

%ded_pm_plot_fluxes('pm/f7/e/01',[200 inf],[3 28]);
%nm='pm/f7/e/01';
%trg=[300 inf];
p=ded_read_param(nm);
c=ded_coord(nm);
b=ded_read_javrg(nm,'ayz',trg,'combine');
if isempty(b)
  return;
end

%  IC/pi = C*x^(1/3), 
%  IU/pi = U*x^(5/3), 
% ICC/pi = C^2/x^(4/3), 
% IUU/pi = U^2*x^(4/3), 
% ICU/pi = C*U
%

if isempty(xrg)
  xrg=[3 p.L-3];
end
x1=xrg(1);
x2=xrg(2);
f =find(c.Jx>=x1 & c.Jx<=x2);
x=c.Jx;
xf=c.Jx(f);

B=0;
BB=0;
hb=isfield(b,'b');
hs=isfield(b,'s');
if hb;    B  = B  +   p.B*b.b      /pi; end; % Total buoyancy
if hs;    B  = B  +   p.S*b.s      /pi; end; % Total buoyancy
if hb;    BB = BB +   p.B*p.B*b.bb /pi; end; % Total buoyancy
if hs;    BB = BB +   p.S*p.S*b.ss /pi; end; % Total buoyancy
if isfield(b,'bs'); BB = BB + 2*p.S*p.B*b.bs /pi; end; % Total buoyancy
%if hb&hs; BB = BB + 2*p.S*p.B*b.bs /pi; end; % Total buoyancy
  
Q  = b.u /pi; % Volume flux
M  = b.uu/pi; % Momentum flux
F  = b.bu/pi; % Buoyancy flux

pn=3+fitg;
pp=[x1-1 ones(1,pn)];

if fitg
  fp = @(pp) ded_MTT(xf-pp(1),pp(2),pp(3),pp(4),pp(5))-[B(f);Q(f);M(f);F(f)];
else
  fp = @(pp) ded_MTT(xf-pp(1),pp(2),p.g,pp(3),pp(4))-[B(f);Q(f);M(f);F(f)];
end

pp=lsqnonlin(fp,pp,[-10 zeros(1,pn)],[x1-1 100*ones(1,pn)]);
if fitg
  [fB fQ fM fF fC]=ded_MTT(xf-pp(1),pp(2),pp(3),pp(4),pp(5));
else
  [fB fQ fM fF fC]=ded_MTT(xf-pp(1),pp(2),p.g,pp(3),pp(4));
  pp=[pp(1:2) p.g pp(3:end)];
end

clf;
ah=jsubplot([1 5],[0.1 0.05],[0.0 0.02],[0.02 0.1]);
xrg=[0 p.L];xr=[x1*[1;1] x2*[1;1]];

axes(ah(1));
plot(x,B,xf,fB);ylabel('B');
maxB=max(B(f));
axis([xrg 0 1.1*maxB]);
line(xr,[[0 0];max(B)*[1 1]],'color',0.7*[1 1 1]);
title(sprintf('%s: x0=%6.3f, F=%5.3f, g=%5.3f, a=%5.3f, b=%5.3f',nm,pp(1),pp(2),pp(3),pp(4),pp(5)));

axes(ah(2));
plot(x,Q,xf,fQ);
maxQ=max(Q(f));
axis([xrg 0 1.1*maxQ]);
ylabel('Q');
line(xr,[[0 0];max(Q)*[1 1]],'color',0.7*[1 1 1]);

axes(ah(3));
maxM=max(M(f));
plot(x,M,xf,fM);ylabel('M');
axis([xrg 0 1.1*maxM]);
line(xr,[[0 0];max(M)*[1 1]],'color',0.7*[1 1 1]);

axes(ah(4));
maxF=max(F(f));
plot(x,F,xf,fF);ylabel('F');
axis([xrg 0 1.1*maxF]);
line(xr,[[0 0];max(F)*[1 1]],'color',0.7*[1 1 1]);

axes(ah(5));
maxC=max(C(f));
plot(x,C,xf,fC);ylabel('BB');
axis([xrg 0 1.1*maxC]);
line(xr,[[0 0];max(C)*[1 1]],'color',0.7*[1 1 1]);

set(ah(1:end-1),'xticklabels',[]);


return;


function plot_all_e
dc;
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/01',[200     inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/01',[341-40  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/02',[351-40  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/03',[148-40  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/04',[ 86-40  inf],[4 26.5]);
figure;ded_pm_fit_fluxes('pm/f7/e/05',[ 20  inf],[4 26.5]);
figure;ded_pm_fit_fluxes('pm/f7/e/06',[ 1  inf],[4 26.5]);
figure;ded_pm_fit_fluxes('pm/f7/e/07',[ 1  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/08',[ 30-20  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/09',[  1     inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/10',[ 76-20  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/11',[233-40  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/12',[ 68-30  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/13',[ 33-10  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/14',[ 15-10  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/15',[ 57-20  inf],[4 26.5]);
% $$$ figure;ded_pm_fit_fluxes('pm/f7/e/16',[  3     inf],[4 26.5]);
