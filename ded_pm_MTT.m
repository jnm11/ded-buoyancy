function ded_pm_MTT(nm,trg,xrg) 
if nargin<3
  xrg=[];
end
if nargin<2
  trg=[];
end
if isempty(trg)
  trg=[-inf inf];
end

%ded_pm_plot_fluxes('pm/f7/e/01',[200 inf],[3 28]);
%nm='pm/f7/e/01';
%trg=[300 inf];
p=ded_read_param(nm);
c=ded_coord(nm);
b=ded_read_javrg(nm,'ayz',trg,'combine');
if isempty(b)
  fns=ded_get_fn(nm,'ayz');
  b=ded_read_hdf(fns{end});
  trg=[b.t1 inf];
  disp(sprintf('Changing time region to %6.2f %6.2f',trg(1),trg(2)));
  b=ded_read_javrg(nm,'ayz',trg,'combine');
end  
%  IC/pi = C*x^(1/3), 
%  IU/pi = U*x^(5/3), 
% ICC/pi = C^2/x^(4/3), 
% IUU/pi = U^2*x^(4/3), 
% ICU/pi = C*U
%


x1=2*min(c.Jx(b.bu>median(b.bu)/2));
x2=2*max(c.Jx(b.bu>median(b.bu)/2))-p.L;
x1=3;
x2=p.L-2;
if isempty(xrg)
  xrg=[3 p.L-2];
end
x1=xrg(1);
x2=xrg(2);

f =find(c.Jx>=x1 & c.Jx<=x2);
x=c.Jx;
xf=c.Jx(f);

% Follows nomenclature of
%Turbulent transport and entrainment in jets and plumes: A DNS study
% Phys Rev Fluids 1  074301 (2016)

Q = b.u/pi;  % Volume flux
M = b.uu/pi; % Momentum flux

rm=Q./sqrt(M);
wm=M./Q;

pn=3;



if 1
  Qfun=@(x,u,a,fp)     36/25*u^1*a^2*max(0,x).^fp;
  Mfun=@(x,u,a,fp)     36/25*u^2*a^2*max(0,x).^fp;
  dQfun=@(x,u,a,fp) fp*36/25*u^1*a^2*max(0,x).^(fp-1);
  dMfun=@(x,u,a,fp) fp*36/25*u^2*a^2*max(0,x).^(fp-1);
  p0=[1 1 1];
  ff= @(p) [Qfun(xf-p(1),p(2),p(3),5/3)-Q(f);Mfun(xf-p(1),p(2),p(3),4/3)-M(f)];
  pp=lsqnonlin(ff,p0,[-10 0 0],[xf(1) inf inf]);
  pp(4)=5/3;
  pp(5)=4/3;
  fQ  =  Qfun(x-pp(1),pp(2),pp(3),pp(4));
  fM  =  Mfun(x-pp(1),pp(2),pp(3),pp(5));
  fdQ = dQfun(x-pp(1),pp(2),pp(3),pp(4));
  fdM = dMfun(x-pp(1),pp(2),pp(3),pp(5));
  frm = fQ./sqrt(fM);
  alpha  = frm./(2*fQ).*fdQ;
  alpha2 = repmat(pp(3),size(x));
else
  pQ   = polyfit(xf,Q(f),pn);
  pM   = polyfit(xf,M(f),pn);
  prm  = polyfit(xf,rm(f),1);
  pdQ  = poly_diff(pQ,1);
  pdrm = poly_diff(prm,1);
  fQ   = polyval(pQ,x);
  fM   = polyval(pM,x);
  frm  = polyval(prm,x);
  fdQ  = polyval(pdQ,x);
  fdrm = polyval(pdrm,x);
  alpha  = frm./(2*fQ).*fdQ;
  alpha2 = 5/6*fdrm;
end

clf;
ah=jsubplot([1 4],[0.1 0.1],[0.0 0.02],[0.02 0.02]);

xrg=[0 p.L];
xr=[x1*[1;1] x2*[1;1]];

axes(ah(1));
plot(x,Q,x,fQ);
Qmax=1.05*max(Q(:));
axis([xrg 0 Qmax]);
line(xr,[[0 0];Qmax*[1 1]],'color',0.7*[1 1 1]);
ylabel('Q');

axes(ah(2));
plot(x,M,x,fM);
Mmax=1.05*max(M(:));
axis([xrg 0 Mmax]);
line(xr,[[0 0];Mmax*[1 1]],'color',0.7*[1 1 1]);
ylabel('M');

axes(ah(3));
plot(x,rm,x,frm,x,Q./sqrt(max(eps,M)),x,fQ./sqrt(max(eps,fM)));
rmmax=1.05*max(rm(:));
axis([xrg 0 rmmax]);
line(xr,[[0 0];rmmax*[1 1]],'color',0.7*[1 1 1]);
ylabel('$r_m$','interpreter','latex');

axes(ah(4));
h=plot(x,alpha,x,alpha2);
alphamax=1.05*max(alpha(:));
alphamax=0.2;
axis([xrg 0 alphamax]);
line(xr,[[0 0];alphamax*[1 1]],'color',0.7*[1 1 1]);
ylabel('$\alpha$','interpreter','latex');
legend(h,{'dQ','dr'});
xlabel(sprintf('$\\alpha=%5.3f$',mean(alpha2)),'interpreter','latex');


return;


ded_pm_MTT('pm/f7/e/01',[200     inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/01',[341-40  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/02',[351-40  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/03',[148-40  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/04',[ 86-40  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/05',[ 65-40  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/06',[ 51-20  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/07',[  9     inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/08',[ 30-20  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/09',[  1     inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/10',[ 76-20  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/11',[233-40  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/12',[ 68-30  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/13',[ 33-10  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/14',[ 15-10  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/15',[ 57-20  inf],[4 26.5]);
ded_pm_MTT('pm/f7/e/16',[  3     inf],[4 26.5]);





CC = b.bb/A;

CU = b.bu/A;

fC  = @(p,x) p(2)*max(0,x-p(1)).^( 1/3);
fU  = @(p,x) p(3)*max(0,x-p(1)).^( 5/3);
fCC = @(p,x) p(4)*max(0,x-p(1)).^(-4/3);
fUU = @(p,x) p(5)*max(0,x-p(1)).^( 4/3);
fCU = @(p,x) p(6)*max(0,x-p(1)).^(   0);

fp = @(p) [fC(p,xf);fU(p,xf);fCC(p,xf);fUU(p,xf);fCU(p,xf);]-[C(f);U(f);CC(f);UU(f);CU(f)];
pp=lsqnonlin(fp,[0 1 1 1 1 1],[-5 zeros(1,5)],[5 100*ones(1,5)]);

wU=(U/pi)./sqrt(UU/pi);
cU=polyfit(x,wU,1);
plot(x,wU,x,polyval(cU,x));disp(cU(1)*5/6);

wC=(C/pi)./sqrt(CC/pi);
cC=polyfit(x,wC,1);
plot(x,wC,x,polyval(cC,x),x,wU,x,polyval(cU,x));disp(wC(1)*5/6);disp(cU(1)*5/6);



clf;
ah=jsubplot([1 5],[0.1 0.1],[0.0 0.02],[0.02 0.02]);
xrg=[0 p.L];xr=[x1*[1;1] x2*[1;1]];
axes(ah(1));
plot(x,C,x,fC(pp,x));ylabel('C');axis([xrg 0 inf]);
line(xr,[[0 0];max(C)*[1 1]],'color',0.7*[1 1 1]);
axes(ah(2));
plot(x,U,x,fU(pp,x));ylabel('U');axis([xrg 0 inf]);
line(xr,[[0 0];max(U)*[1 1]],'color',0.7*[1 1 1]);
axes(ah(3));
plot(x,CU,x,fCU(pp,x));ylabel('CU');axis([xrg 0 inf]);
line(xr,[[0 0];max(CU)*[1 1]],'color',0.7*[1 1 1]);
axes(ah(4));
plot(x,CC,x,fCC(pp,x));ylabel('CC');axis([xrg 0 inf]);
line(xr,[[0 0];max(CC)*[1 1]],'color',0.7*[1 1 1]);
axes(ah(5));
plot(x,UU,x,fUU(pp,x));ylabel('UU');axis([xrg 0 inf]);
line(xr,[[0 0];max(UU)*[1 1]],'color',0.7*[1 1 1]);
set(ah(1:end-1),'xticklabels',[]);

return;


fC  = @(p,x) p(2)*max(0,x-p(1)).^( 1/3);
fU  = @(p,x) p(3)*max(0,x-p(1)).^( 5/3);
fCC = @(p,x) p(2)^2*max(0,x-p(1)).^(-4/3);
fUU = @(p,x) p(3)^2*max(0,x-p(1)).^( 4/3);
fCU = @(p,x) p(2)*p(3)*max(0,x-p(1)).^(   0);

fp = @(p) [fC(p,x);fU(p,x);fCC(p,x);fUU(p,x);fCU(p,x);]-[C;U;CC;UU;CU];
pp=lsqnonlin(fp,[2 1 1],[0 0 0],[10 10 10]);

return;

a=ded_read_javrg('pm/f7/e/01','a',[320 inf]);
c=ded_coord('pm/f7/e/01');
for j=1:c.NJx
  mesh(c.Jy,c.Jz,a.b(:,:,j));
  title(sprintf('%u',j));
  drawnow;
end




ded_pm_fit_fluxes('pm/f7/e/01',[200 inf],[3 28]);
ded_pm_fit_fluxes('pm/f7/e/11',[200 inf],[3 27]);