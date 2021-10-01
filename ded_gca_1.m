function [a b p]=ded_gca_1(nm,trg)
% ded_gca_1 analysis of gravity current
% analyses averaged data files  
%a=ded_gc_read(NM,'finalt');
%nm='gc/f6/g/22';
%nm='gc/test/50';
if nargin<2
  trg=[];
end
if nargin<1
  disp('ded_gca_1: requires at least one input argument');
  a=[];b=[];p=[];
  return;
end

f=ded_read_g(nm,'force');
p=ded_read_param(nm);
c=ded_coord(nm);
[a fld]=ded_read_javrg(nm,'ay',trg);
if isempty(a)
  [a fld]=ded_read_javrg(nm,'a',trg);
end
if isempty(a)
  disp(sprintf('ded_gca_1: no average data for %s [%6.3f %6.3f]',nm,trg(1),trg(2)));
  b=[];
  return;
end

  
if p.Ny>1
  for j=1:length(fld)
    a.(fld{j})=a.(fld{j})/p.W;
  end
end

ue=mean(mean(a.uu-a.u.^2));

if isempty(a)
  disp(sprintf('ded_gca_1: No averaged data for %s',nm));
  return;
end

[a b]=ded_gc_diff_int(a,p,c.dJx);
ue(2)=mean(mean(a.uu-a.u.^2));


PS=p.B*p.g*p.hb;
nu= 1/p.Re;
%xdmax=max([p.x0 p.x1 p.x2 p.x3 p.x4 p.x5 p.x6  p.x7]);
xrg=[0 p.L];
prg=[0 p.L -0.5 0.5];

if p.forcing==7
  div=abs(sum(b.cm,1));
  divtol=10*median(div);
  fx=find(div<divtol);
  fx=min(fx):max(fx);
else  
  xdmax=3.5;
  fx=find(c.Jx>xdmax);
end

x=c.Jx(fx);
z=a.z;
xrg=[x(1) x(end)];
zrg=[0 p.H];
fz=[findmin(abs(c.Jx-3)) fx(end)];

a.bIz   = a.bIz   - a.bIz(  end,:);
b.cpxIz = b.cpxIz - b.cpxIz(end,:);
b.cpzIz = b.cpzIz - b.cpzIz(end,:);

a.p     = a.p     - a.p(    :,fx(end));
a.P     = a.P     - a.P(    :,fx(end));
b.cpxIx = b.cpxIx - b.cpxIx(:,fx(end));
b.cpzIx = b.cpzIx - b.cpzIx(:,fx(end));

hsp=-p.g*a.bIz; % Hydrostatic pressure inegrated down
b.hs=a.p-hsp; 
%b.hs=b.hs(end,:)-b.hs;
b.Hs=a.P-hsp; % Hydrostatic pressure inegrated down
              %b.Hs=b.Hs(end,:)-b.Hs;

e.cpx.z.max=max(max(abs(b.cpx(:,fz))));
e.cpz.z.max=max(max(abs(b.cpz(:,fz))));
e.npx.z.max=max(max(abs(b.npx(:,fz))));
e.npz.z.max=max(max(abs(b.npz(:,fz))));
e.cpx.z.rms=sqrt(mean(mean(b.cpx(:,fz).^2)));
e.cpz.z.rms=sqrt(mean(mean(b.cpz(:,fz).^2)));
e.npx.z.rms=sqrt(mean(mean(b.npx(:,fz).^2)));
e.npz.z.rms=sqrt(mean(mean(b.npz(:,fz).^2)));

e.cpx.x.max=max(max(abs(b.cpx([1 end],fx))));
e.cpz.x.max=max(max(abs(b.cpz([1 end],fx))));
e.npx.x.max=max(max(abs(b.npx([1 end],fx))));
e.npz.x.max=max(max(abs(b.npz([1 end],fx))));
e.cpx.x.rms=sqrt(mean(mean(b.cpx([1 end],fx).^2)));
e.cpz.x.rms=sqrt(mean(mean(b.cpz([1 end],fx).^2)));
e.npx.x.rms=sqrt(mean(mean(b.npx([1 end],fx).^2)));
e.npz.x.rms=sqrt(mean(mean(b.npz([1 end],fx).^2)));


figure;clf;
gu_title('Conservative Momentum equation')
ah=jsubplot([1 2],[0.1 0.05],[0.01 .01],[0.01 0.05]);
axes(ah(1));h=plot(x,  b.cpx([1 end],fx));line(x,0*x,'color',0.7*[1 1 1]);
gu_setalim(0,0.1);ylabel(sprintf('rms   x   max%s%8.1e %8.1e',10,e.cpx.x.rms,e.cpx.x.max));     
axes(ah(2));  plot(x,  b.cpz([1 end],fx));line(x,0*x,'color',0.7*[1 1 1]);
gu_setalim(0,0.1);ylabel(sprintf('rms   z   max%s%8.1e %8.1e',10,e.cpz.x.rms,e.cpz.x.max));      
% $$$ axes(ah(3));  plot(x,b.cpxIx([1 end],fx));gu_setalim(0,0.1); ylabel('int x');
% $$$ axes(ah(4));  plot(x,b.cpzIx([1 end],fx));gu_setalim(0,0.1); ylabel('int z');
axes(ah(1));
legend(h,'bottom','top','location','ne');
set(ah(1:end-1),'xticklabels',[]);
set(ah,'yticklabels',[]);

figure;clf;
gu_title('Non-Conservative Momentum equation');
ah=jsubplot([1 2],[0.1 0.05],[0.02 .02],[0.01 0.05]);
axes(ah(1));h=plot(x,  b.npx([1 end],fx));line(x,0*x,'color',0.7*[1 1 1]);
gu_setalim(0,0.1);ylabel(sprintf('rms   x   max%s%8.1e %8.1e',10,e.npx.x.rms,e.npx.x.max));     
axes(ah(2));  plot(x,  b.npz([1 end],fx));line(x,0*x,'color',0.7*[1 1 1]);
gu_setalim(0,0.1);ylabel(sprintf('rms   z   max%s%8.1e %8.1e',10,e.npz.x.rms,e.npz.x.max));      
axes(ah(1));
legend(h,'bottom','top','location','ne');
set(ah(1:end-1),'xticklabels',[]);
set(ah,'yticklabels',[]);




figure;clf;
gu_title('Conservative Momentum equation')
%ah=jsubplot([4 1],[0.1 0.1],[0.01 .01],[0.01 0.05]);
ah=jsubplot([2 1],[0.1 0.1],[0.05 .05],[0.01 0.05]);
axes(ah(1));h=plot(  b.cpx(:,fz),z);gu_setalim(0.1,0);
line(0*z,z,'color',0.7*[1 1 1]);
xlabel(sprintf('rms   x   max\n%8.1e %8.1e',e.cpx.z.rms,e.cpx.z.max));
axes(ah(2));  plot(  b.cpz(:,fz),z);gu_setalim(0.1,0);
line(0*z,z,'color',0.7*[1 1 1]);
xlabel(sprintf('rms   z   max\n%8.1e %8.1e',e.cpz.z.rms,e.cpz.z.max));
%axes(ah(3));  plot(b.cpxIz(:,fz),z);gu_setalim(0.1,0);xlabel('int x');
%axes(ah(4));  plot(b.cpzIz(:,fz),z);gu_setalim(0.1,0);xlabel('int z');
set(ah(2:end),'yticklabels',[]);
axes(ah(2));
legend(h,'left','right','location','ne');


figure;clf;
gu_title('Non-Conservative Momentum equation');
ah=jsubplot([2 1],[0.1 0.1],[0.05 .05],[0.01 0.05]);
axes(ah(1));h1=plot(  b.npx(:,fz),z);
line(0*z,z,'color',0.7*[1 1 1]);
gu_setalim(0.1,0);xlabel(sprintf('rms   x   max%s%8.1e %8.1e',10,e.npx.z.rms,e.npx.z.max));
axes(ah(2));h2=plot(  b.npz(:,fz),z);
line(0*z,z,'color',0.7*[1 1 1]);
gu_setalim(0.1,0);xlabel(sprintf('rms   z   max%s%8.1e %8.1e',10,e.npz.z.rms,e.npz.z.max));
set(ah(2:end),'yticklabels');
set(ah,'xticklabels',[]);
axes(ah(2));
legend(h1,{sprintf('left %7.1e',max(abs(b.npx(:,fz(1))))),sprintf('right %7.1e',max(abs(b.npx(:,fz(2)))))},'location','ne');
legend(h2,{sprintf('left %7.1e',max(abs(b.npz(:,fz(1))))),sprintf('right %7.1e',max(abs(b.npz(:,fz(2)))))},'location','ne');


figure;clf
ah=jsubplot([1 4],[0.1 0.05],[0.01 .01],[0.01 0.05]);
b.vuu=(a.duudx+a.duudz)/p.Re;
b.vww=(a.dwwdx+a.dwwdz)/p.Re;
b.duw=a.udwdx-a.wdudz;
b.vorz=a.dudz-a.dwdx;
minP=min(min(a.P([1 end],fx)));
maxP=max(max(a.P([1 end],fx)));
rgP= (minP+maxP)/2+(maxP-minP)*1.05/2*[-1 1];
h=jplot(ah(1),x,  a.p([1 end],fx),'p');       
jplot(ah(2),x,  a.P([1 end],fx),'P');       
jplot(ah(3),x,  a.u([1 end],fx),'u');
jplot(ah(4),x,  b.vuu([1 end],fx),'visc');
legend(h,'bottom','top','location','ne');




figure;clf;
ah=jsubplot([5 1],[0.1 0.1],[0.01 .01],[0.01 0.05]);
axes(ah(1));h=plot(   b.hs(:,fz)-b.hs(end,fz),z);gu_setalim(0.1,0);title('b.hs');
axes(ah(2));plot(a.ww(:,fz)/2,z);gu_setalim(0.1,0);title('ww/2');
axes(ah(3));plot(a.udwdxIz(:,fz),z);gu_setalim(0.1,0);title('Iudwdx');
axes(ah(4));plot(  b.vww(:,fz),z);gu_setalim(0.1,0);title('vww');
XX=b.hs(:,fz)+a.ww(:,fz)/2+a.udwdxIz(:,fz)+b.vww(:,fz);
axes(ah(5));plot(XX,z);gu_setalim(0.1,0);title('residual');
set(ah(2:5),'yticklabels',[]);
legend(h,'left','right','location','best');

figure;clf;
ah=jsubplot([5 1],[0.1 0.1],[0.01 .01],[0.01 0.05]);
axes(ah(1));h=plot(a.u(:,fz)+p.U                  ,z);title('u-U');   gu_setalim(0.1,0);
axes(ah(2));  plot(a.w(:,fz)                      ,z);title('w');     gu_setalim(0.1,0);
axes(ah(3));  plot(a.uu(:,fz)-a.u(:,fz).*a.u(:,fz),z);title('uu-u*u');gu_setalim(0.1,0);
axes(ah(4));  plot(a.uw(:,fz)-a.u(:,fz).*a.w(:,fz),z);title('uw-u*w');gu_setalim(0.1,0);
axes(ah(5));  plot(a.ww(:,fz)-a.w(:,fz).*a.w(:,fz),z);title('ww-w*w');gu_setalim(0.1,0);
set(ah(2:5),'yticklabels',[]);
legend(h,'left','right','location','best');


e=[];

figure;clf
ah=jsubplot([1 4],[0.1 0.05],[0.01 .01],[0.01 0.05]);
%axes(ah(1));plot(x,a.dudx(:,fx)+a.dwdz(:,fx),x,b.cm(:,fx));gu_setalim(0,0.1); ylabel('div');
jplot(ah(1),c.Jx,a.dudx+a.dwdz-a.d,'div');
jplot(ah(2),c.Jx,a.uu-a.u.*a.u,'uu');
jplot(ah(3),c.Jx,a.ww-a.w.*a.w,'ww');
jplot(ah(4),c.Jx,a.uw-a.u.*a.w,'uw');
set(ah(1:end-1),'xticklabels',[]);set(ah,'ytick',[]);
axes(ah(1));title('Consistency laminar stationary');

figure;clf;
ah=jsubplot([1 4],[0.1 0.05],[0.01 .01],[0.01 0.05]);
jplot(ah(1),c.Jx, a.ududx-a.duudx/2,'ududx');       
jplot(ah(2),c.Jx, a.ududz-a.duudz/2,'ududz');       
jplot(ah(3),c.Jx, a.wdwdx-a.dwwdx/2,'wdwdx');       
jplot(ah(4),c.Jx, a.wdwdz-a.dwwdz/2,'wdwdz');       
set(ah(1:end-1),'xticklabels',[]);set(ah,'ytick',[]);
axes(ah(1));title('Quad derivatives 1')

figure;clf; % This is bad especially 2& 4
ah=jsubplot([1 4],[0.1 0.05],[0.01 .01],[0.01 0.05]);
jplot(ah(1),c.Jx, a.udwdz+a.ududx-a.u.*a.d,'u del u');       
jplot(ah(2),c.Jx, a.wdudx+a.wdwdz-a.w.*a.d,'w del u');       
jplot(ah(3),c.Jx, a.udwdx+a.wdudx-a.duwdx,'duwdx');       
jplot(ah(4),c.Jx, a.udwdz+a.wdudz-a.duwdz,'duwdz');       
set(ah(1:end-1),'xticklabels',[]);set(ah,'ytick',[]);
axes(ah(1));title('Quad derivatives 2')
% $$$  keyboard
% $$$ f1=a.udwdx+a.wdudx;
% $$$ f2=a.duwdx;
% $$$ plot(c.Jx,a.udwdx,c.Jx,a.u.*a.dwdx)
% $$$ plot(c.Jx,a.udwdx-a.u.*a.dwdx)
% $$$ plot(c.Jx,a.wdudx-a.w.*a.dudx)
% $$$ 
% $$$ plot(c.Jx,a.wdwdx-a.w.*a.dwdx)
% $$$ 
% $$$  C=sum(sum(f1.*f2))/sum(sum(f1.^2));
 
 
figure;clf; % This is good
ah=jsubplot([1 4],[0.1 0.05],[0.01 .01],[0.01 0.05]);
jplot(ah(1),c.Jx, a.ududx-a.u.*a.dudx,'ududx');       
jplot(ah(2),c.Jx, a.ududz-a.u.*a.dudz,'ududz');       
jplot(ah(3),c.Jx, a.udwdx-a.u.*a.dwdx,'udwdx');       
jplot(ah(4),c.Jx, a.udwdz-a.u.*a.dwdz,'udwdz');       
set(ah(1:end-1),'xticklabels',[]);set(ah,'ytick',[]);
axes(ah(1));title('Quad derivatives stationary 1')


figure;clf; % This is good
ah=jsubplot([1 4],[0.1 0.05],[0.01 .01],[0.01 0.05]);
jplot(ah(1),c.Jx, a.wdudx-a.w.*a.dudx,'wdudx');       
jplot(ah(2),c.Jx, a.wdudz-a.w.*a.dudz,'wdudz');       
jplot(ah(3),c.Jx, a.wdwdx-a.w.*a.dwdx,'wdwdx');       
jplot(ah(4),c.Jx, a.wdwdz-a.w.*a.dwdz,'wdwdz');       
set(ah(1:end-1),'xticklabels',[]);set(ah,'ytick',[]);
axes(ah(1));title('Quad derivatives stationary 2')

figure;clf; % This is good
ah=jsubplot([1 4],[0.1 0.05],[0.01 .01],[0.01 0.05]);


jplot(ah(1),c.Jx, a.dudx-a.cdudx,'dudx');
jplot(ah(2),c.Jx, a.dudz-a.cdudz,'dudz');              
jplot(ah(3),c.Jx, a.dwdx-a.cdwdx,'dwdx');              
jplot(ah(4),c.Jx, a.dwdz-a.cdwdz,'dwdz');              
set(ah(1:end-1),'xticklabels',[]);
set(ah,'ytick',[]);
axes(ah(1));title('Single derivatives')

return;


function h=jplot(a,x,f,s)
axes(a);
h=plot(x,f);
axis('tight');
aa=axis;
e=max(abs(aa([3:4])));
ylabel(sprintf('%s %5.0e',s,e));       
return;



[a b p]=ded_gca_1('gc/f7/g/1300/1000');


fx =find(c.x>2 & c.x<9.5);
x=c.x(fx);
P=a.p+a.uu/2+a.ww/2;
plot(x,a.P(end,fx));

a.hp=a.bIz*p.g;

p0=a.p(end,fx(end));
p1=a.p(end,fx(1))-p0;
p2=a.p(1,fx(1))-p0;
p3=a.p(1,fx(end))-p0;
[ps fp]=max(a.p(1,fx))-p0;
p0=p0-p0;

P0=a.P(end,fx(end));
P1=a.P(end,fx(1))-P0;
P2=a.P(1,fx(1))-P0;
P3=a.P(1,fx(end))-P0;
Ps=a.P(1,fp)-P0;
P0=P0-P0;
u0=a.u(end,fx(end));
u1=a.u(end,fx(1));
u2=a.u(1,fx(1));
u3=a.u(1,fx(end));
us=a.u(1,fx(fp));

hp12=a.hp(end,fz(1))-a.hp(1,fz(1));
hp03=a.hp(end,fz(end))-a.hp(1,fz(end));

E01=p0+u0^2/2-p1-u1^2/2;
E32=p3+u3^2/2-p2-u2^2/2;


fz=fx([1 end]);
plot(a.hp(:,fz)+a.p(:,fz)-[p1 p0],c.z);

plot(x,P(1,fx)-P(1,fx(end)));

E1=P(1,fx(1))-P(1,fx(end));
E2=p2-p1-hp12;

lp=p1+hp12;
rp=p3+p.U^2/2;

p2+E1-rp

if 1
  disp(sprintf('p0: %6.3f, p1: %6.3f, p2: %6.3f, ps: %6.3f, p3: %6.3f, ',p0, p1, p2, ps, p3));
  disp(sprintf('u0: %6.3f, u1: %6.3f, u2: %6.3f, us: %6.3f, u3: %6.3f, ',u0, u1, u2, us, u3));
  disp(sprintf('ps-p3 %6.3f Ps-P3     %6.3f',ps-p3,Ps-P3));
  disp(sprintf('p2-ps %6.3f P2-Ps     %6.3f',p2-ps,P2-Ps));
  disp(sprintf('p1-p2 %6.3f p1-p2+bIz %6.3f',p1-p2,p1-p2+hp12));
  disp(sprintf('p0-p1 %6.3f P0-P1     %6.3f',p0-p1,P0-P1));
  disp(sprintf('p3-p0 %6.3f p3-p0     %6.3f',p3-p0,p3-p0));
end


viscxx=diff(a.dudx(1,[fx(1) fx(fp) fx(end)])/p.Re);
viscxz=diff(a.dudzIx(1,[fx(1) fx(fp) fx(end)])/p.Re);


plot(x,Mxa(:,fx))


Bx=a.Pdx(1,:)-nu*(a.udxx(1,:)+a.udzz(1,:));

Mzb=a.uwdx+a.pdz-nu*(a.wdxx+a.wdzz)+p.g*a.b;

plot(c.x,a.uudx-duudx);
plot(c.x,a.uwdz-duwdz);
plot(c.x,a.pdx-dpdx);



ndiff=1;
dwdz=ichebdifff(a.w,dimz,p.H,ndiff);hist(dudx(:)+dwdz(:),1000);div=dudx+dwdz;imagesc(abs(div)>0.01);


f1=a.udwdx;
f2=pr_diff(a.uw,dx,2)-ichebdifff(a.ww/2,1,p.H,1);
plot(f1(:),f2(:),'.');

a.dudx=pr_diff(a.u,dx,2);
a.dwdx=pr_diff(a.w,dx,2);

if isempty(a)
  return;
end


dx=c.x(2)-c.x(1);

 


fnm={'b','u','p','pdz','uu','uudz','vv','vvdz','ww'};
for j=1:length(fnm)
  L.(fnm{j})=a.(fnm{j})(1,:);
  R.(fnm{j})=a.(fnm{j})(end,:);
end
fnm={'uI','vI','bI','pI','uuI','vvI','wwI'};
for j=1:length(fnm)
  I.(fnm{j})=(a.(fnm{j})(end,:)-a.(fnm{j})(1,:));
end
L.udx=pr_diff(L.u,dx,1);
R.udx=pr_diff(R.u,dx,1);
L.bdx=pr_diff(L.b,dx,1);
R.bdx=pr_diff(R.b,dx,1);
L.udxx=pr_diff(L.udx,dx,1);
R.udxx=pr_diff(R.udx,dx,1);

p0=a.p(end,end);
uu0=a.uu(end,end);
p0=R.p(end);
uu0=R.uu(end,end);
P0=p0+uu0^2/2;
L.P=L.p+L.uu/2-P0;
R.P=R.p+R.uu/2-P0;


P=a.p+a.uu/2+a.ww/2;
P0=P(end,end);
P=P-P0;



figure(1);
h=plot(x,P(1,fx),x,P(end,fx));
legend(h,'lower','upper');
xlabel('x');
ylabel('P');
axis('tight');


figure(2);
h=plot(x,P(1,fx)-P(end,fx),x,a.p(1,fx)-R.p(end,fx),x,p.g*I.bI(fx));
legend(h,{'dP','dp','bI'});
xlabel('x');
ylabel('P');
axis('tight');

figure(3);
h=plot((a.p(1,fx)-R.p(end,fx)-p.g*I.bI(fx))/PS);
legend(h,{'dP','dp','bI'});
xlabel('x');
ylabel('P');
axis('tight');


figure(3);
S=a.S(:,fx);
imagesc(x,z,S);
set(gca,'ydir','normal');
title('Dissipation');



%This should be identically zero
%plot(c.x,a.L.uudx/2+a.L.uwdz)

% $$$ fnm=setdiff(fieldnames(L),{'x','b','t','nt'});
% $$$ for j=1:length(fnm);
% $$$   L.(fnm{j})=L.(fnm{j})/PS;
% $$$   R.(fnm{j})=R.(fnm{j})/PS;
% $$$ end
% $$$   


x=c.x;
sfigure(1);clf
LP=L.p+L.uu/2;LP=LP(fx)-LP(end);
RP=R.p+R.uu/2;RP=RP(fx)-RP(end);
Ldu=pr_diff(L.u,dx,1);Ldu=Ldu(fx)/Re;
Rdu=pr_diff(R.u,dx,1);Rdu=Rdu(fx)/Re;
Rv=(pr_diff(R.u,dx,1)+pr_int(R.udzz,dx,1))/Re;Rv=Rv(fx)-Rv(end);
Lv=(pr_diff(L.u,dx,1)+pr_int(L.udzz,dx,1))/Re;Lv=Lv(fx)-Lv(end);
LNS=LP-Lv;LNS=LNS-LNS(end);
RNS=RP-Rv;RNS=RNS-RNS(end);

subplot(3,1,1);
%h=plot(x,L.p(fx),x,L.uu(fx)/2,x,L.vv(fx)/2,x,LP,x,Ldu);
%legend(h,{'p','u^2/2','v^2/2','P','u_x'},'location','best','interpreter','tex');
if isfield(L,'vv')
  h=plot(x,L.vv(fx)/2,'r',x,LP,'b',x,Ldu,'g--');
  legend(h,{'v^2/2','P','u_x'},'location','best','interpreter','tex');
else
  h=plot(x,LP,'b',x,Ldu,'g--');
  legend(h,{'P','u_x'},'location','best','interpreter','tex');
end

ylabel('bottom');title('Bernouilli');
gu_setalim
subplot(3,1,2);
if isfield(R,'vv')
  plot(x,R.vv(fx)/2,'r',x,RP,'b',x,Rdu,'g--');
else
  plot(x,RP,'b',x,Rdu,'g--');
end

%plot(x,R.p(fx),x,R.uu(fx)/2,x,R.vv(fx)/2,x,RP,x,Rdu);
ylabel('top');
gu_setalim;
subplot(3,1,3);cla
line(x,x*0,'linestyle','--','color',.7*[1 1 1]);
hold('on');
% $$$ h=plot(x,RP,'b',x,LP,'r-.',x,Rdu,x,Ldu);
% $$$ legend(h,{'T','B','dT','DB'},'location','best');
h=plot(x,RP,'b',x,LP,'r-.');
legend(h,{'T','B'},'location','best');
gu_setalim;

if isfield(a.y,'E')
  E=a.y.E(:,fx);
elseif  isfield(a.y,'Ey')
  E=a.y.Ey(:,fx);
else
  E=[];
end
if ~isempty(E)
  sfigure(2);clf
  imagesc(a.y.x,a.y.z,E);
  set(gca,'ydir','normal','position',[0.01 0.01 0.98 0.9],'box','on');
  set(gca,'clim',percentile(E(:),[1 99]));
  title('Dissipation rate')
  gu_setalim;
end


sfigure(3);clf;
subplot(3,1,1);
h=plot(x,LP,x,Lv);
gu_setalim;
legend(h,{'P','viscous'});
ylabel('bottom');

subplot(3,1,2);
h=plot(x,RP,x,Rv);
gu_setalim;
ylabel('top');

subplot(3,1,3);cla
line(x,x*0,'linestyle','--','color',.7*[1 1 1]);
hold('on')
h=plot(x,LNS,'b',x,RNS,'r-.');
legend(h(1:2),{'bottom','top'},'location','best');
title('x momentum along top and bottom');
gu_setalim;
