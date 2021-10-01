function b=ded_gc_extrap_uw(nm,typ,fitrg,plotrg,showb,pp,sfn)
% extrapolate uh,uw,uq,bh,bw to x is zero using serf fits on final state
% f fit range
% ff plot range

if nargin<2
  plotrg=[];
end
if nargin<3
  fitrg=[];
end
if nargin<4
  plotrg=[];
end
if nargin<5
  showb=[];
end
if nargin<6
  sfn=[];
end
if iscell(nm)
  b=[];
  for j=1:length(nm)
    bb=ded_gc_extrap_uw(nm{j},typ,fitrg,plotrg,showb,pp,sfn);
    b=struct_array_append(b,bb,nm{j});
  end
  return;
end

[b fne]=ded_gc_fit_erf_profiles(nm,[],typ,false,false);

if ~isempty(sfn)
  fnmat = [ded_dedalus_data_dir '/results/' nm '/' sfn '.mat'];
  if file_nt(fnmat,fne)
    load(fnmat);
    return;
  end
end

if isempty(b)
  return;
end
if isempty(showb)
  showb=false;
end
p=ded_read_param(nm);
H=p.H;
u=-p.V;
if isfield(p,'hu')
  hu=p.hu;
else
  hu=[];
end



if isempty(plotrg)
  f1=find((b.x<p.L/2)& abs (b.uI-median(b.uI))<1e-3);
  f1=min(f1):max(f1);
else
  f1= find(b.x>=plotrg(1)&b.x<=plotrg(2));
end
x1=b.x(f1);

if isempty(fitrg)
  f2 = f1(1:round(length(f1)/2));
else
  f2= find(b.x>fitrg(1)&b.x<=fitrg(2));
end

x2=b.x(f2);

 
  %subplot(3,1,1);plot(x,b.uw(f2),'.', x,b.fuw(x));ylabel('uw');
%subplot(3,1,2);plot(x,b.uh(f2),'.', x,b.fuh(x));ylabel('uh');
%subplot(3,1,3);plot(x,b.u1(f2),'.', x,b.fu1(x),x,f(pp,x,0));ylabel('u1');
%fs=100;plot(z,c.u(:,fs),z,f(pp,x(fs)));
  

ip.nbw=2;
ip.nbh=1;
ip.nuw=2;
ip.nuh=1;
ip.nu1=1;
pp=combine_struct(ip,pp);

puw=polyfit(x2,b.uw(f2),pp.nuw);
puh=polyfit(x2,b.uh(f2),pp.nuh);
pbw=polyfit(x2,b.bw(f2),pp.nbw);
pbh=polyfit(x2,b.bh(f2),pp.nbh);
pu1=polyfit(x2,b.u1(f2),pp.nu1);

b.x1=b.x(f2(1));
b.x2=b.x(f2(end));
x3=linspace(0,x1(end),100);
b.param = p;

fu1=@(p,x) polyval(p,     x);
if isempty(hu)
  fuh=@(p,x) polyval(p,x);
else
  fuh=@(p,x) polyval([p hu],x);
end
fuh=@(p,x) polyval(p,x);
fuw=@(p,x) polyval(p,     x);
fbh=@(p,x) polyval(p,     x);
fbw=@(p,x) polyval(p,     x);

b.euw=fuw(puw,0);
b.euh=fuh(puh,0);
b.ebw=fbw(pbw,0);
b.ebh=fbh(pbh,0);
b.eu1=fu1(pu1,0);

rguh=          1:(pp.nuh+1);
rguw=rguh(end)+(1:pp.nuw+1);
rgu1=rguw(end)+(1:pp.nu1+1);

if ~isfield(p,'wu');p.wu=NaN;end
if ~isfield(p,'hu');p.hu=NaN;end
if ~isfield(p,'U1');p.U1=NaN;end
if ~isfield(p,'wb');p.wb=NaN;end
if ~isfield(p,'hb');p.hb=NaN;end

b.nm=nm;

if true
  c=load(sprintf('%s/results/%s/ay.mat',ded_dedalus_data_dir,nm));
  c=c.a;
  w=ichebintw(size(c.u,1));w=w(:);

  cc=ded_coord(nm);
  z=cc.z;
  
  f = @(pp,x,z) -p.V+(fu1(pp(rgu1),x)+p.V).*serf1(z,fuh(pp(rguh),x),fuw(pp(rguw),x),H);
  fff = @(pp)  reshape(w.*(f(pp,x2,z)-c.u(:,f2)),length(z)*length(f2),1);
  
  pp0 = [puh puw pu1];
  disp(pp0);
  disp(pp0(rguh)-puh);
  disp(pp0(rguw)-puw);
  disp(pp0(rgu1)-pu1);
  
  lb=repmat(-inf,1,length(pp0));
  ub=repmat( inf,1,length(pp0));
  ub(rgu1(end-1))=0;
  %pp=lsqnonlin(fff,pp0,lb,ub);
  pp=pp0;
  e1=sqrt(mean(fff(pp0).^2));
  e2=sqrt(mean(fff(pp).^2));
  disp([e1 e2]);
  b.fuw = @(x)fuw(pp(rguw),x);
  b.fuh = @(x)fuh(pp(rguh),x);
  b.fu1 = @(x)fu1(pp(rgu1),x);
  b.guw=b.fuw(0);
  b.guh=b.fuh(0);
  b.gu1=b.fu1(0);
  if false;
    for F=f2(1):10:f2(end)
      fF= -p.V+b.uB(F)*serf1(z,b.uh(F),b.uw(F),p.H);
      plot(z,f(pp0,b.x(F),z),z,f(pp,b.x(F),z),z,c.u(:,F),z([1 end]),fu1(pp(rgu1),b.x(F)*[1 1]));
      drawnow;
    end;
    %disp([fu1(pp0(rgu1),x2(1)) fu1(pp(rgu1),x2(1))]);
    clf;
    subplot(3,1,1);plot(x2,polyval(puw,x2),x2,fuw(pp(rguw),x2));
    subplot(3,1,2);plot(x2,polyval(puh,x2),x2,fuh(pp(rguh),x2));
    subplot(3,1,3);plot(x2,polyval(pu1,x2),x2,fu1(pp(rgu1),x2));
    keyboard
    plot(z,f(pp,x2,z),z,c.u(:,f2));
  end

else
  b.fuw = @(x)NaN;
  b.fuh = @(x)NaN;
  b.fu1 = @(x)NaN;
end
b.gbw=NaN;
b.gbh=NaN;

clf;
if showb
  ah=jsubplot([1 5],[0.10 0.05],[0.01 0.01],[0.01 0.05]);
  axes(ah(4));
  plot(x1,b.bw(f1),x3,fbw(pbw,x3),'--');
  ylabel('bw');
  axis('tight');aa=axis;axis([aa(1:3) 1.1*aa(4)-0.1*aa(3)]);
  axes(ah(5));
  plot(x1,b.bh(f1),x3,fbh(pbh,x3),'--');
  ylabel('bh');
  axis('tight');aa=axis;axis([aa(1:3) 1.1*aa(4)-0.1*aa(3)]);
else
  ah=jsubplot([1 3],[0.10 0.05],[0.01 0.01],[0.01 0.05]);
end

axes(ah(1));
plot(x1,b.uw(f1),x3,fuw(puw,x3),'--',x3,b.fuw(x3));
ylabel('uw');
axis('tight');aa=axis;axis([aa(1:3) 1.1*aa(4)-0.1*aa(3)]);
title(nm);


axes(ah(2));
plot(x1,b.uh(f1),x3,fuh(puh,x3),'--',x3,b.fuh(x3));
ylabel('uh');
axis('tight');aa=axis;axis([aa(1:3) 1.1*aa(4)-0.1*aa(3)]);


[w1 w2]=ichebendsw(length(z));

axes(ah(3));
plot(x1,b.u1(f1),x3,fuw(pu1,x3),'--',x3,b.fu1(x3));%,x1,w1*c.u(:,f1));
ylabel('u1');
axis('tight');aa=axis;axis([aa(1:3) 1.1*aa(4)-0.1*aa(3)]);
set(ah(1:end-1),'xticklabels',[]);
drawnow;
disp('param local global');
disp(sprintf('uw %6.4f %6.4f %6.4f',p.wu,b.euw,b.guw));
disp(sprintf('uh %6.4f %6.4f %6.4f',p.hu,b.euh,b.guh))
disp(sprintf('u1 %6.4f %6.4f %6.4f',p.U1,b.eu1,b.gu1));
disp(sprintf('bw %6.4f %6.4f %6.4f',p.wb,b.ebw,b.gbw));
disp(sprintf('bh %6.4f %6.4f %6.4f',p.hb,b.ebh,b.gbh));


if ~isempty(sfn)
 save(fnmat,'b');
end
drawnow;

return;

b=ded_gc_extrap_uw('gc/f7/mg/4000/0250/2304/36593',[0.1 35],[5 30]);

hn=(b.uh-b.bh)./sqrt(b.uw.^2+b.bw.^2);
hp=(b.uh+b.bh)./sqrt(b.uw.^2+b.bw.^2);
plot(b.x,hn.*erf(hn),b.x,hp.*erf(hp));

q=(sqrt(pi).*(erf(hp).*hp-hn.*erf(hn))+exp(-hp.^2)-exp(-hn.^2))./(sqrt(b.bw.^2+b.uw.^2).*sqrt(pi).*(hn.^2-hp.^2));
plot(b.x,b.q);


gc/f6/i/4000/0250/4608/xxxxb

fns=cellstr_ls('gc/f6/i/*/*/2*/xxxxb')
%ded_make_avg(fns,'ay');
ded_gc_fit_erf_profiles(fns,[],'ray');

cd('~/');
a =ded_gc_extrap_uw(cellstr_ls('gc/f6/i/*/*/2*/xxxxb'),'ray',[5 20],[]);

ded_gc_fit_erf_profiles(nm,[],'ray',true,true);
a=ded_gc_extrap_uw(nm,'ray',[1 10],[0.1 30]);

[a b]=ded_gc_fit_erf_profiles(nm,[],'ray');
plot(a.x,a.u1);


 p=[a.param];
  Pe=round([p.Scb].*[p.Re]);
  Re=[p.Re];
  clf;
  subplot(3,1,1);
  [h lh s]=groupplot(Re,[a.euw],Pe);
  legend(lh,cellsprintf('%.0f',s));
  xlabel('Re');
  ylabel('uw');
  
  subplot(3,1,2);
  [h lh s]=groupplot(Re,[a.euh],Pe);
  legend(lh,cellsprintf('%.0f',s));
  xlabel('Re');
  ylabel('uh');

  subplot(3,1,3);
  [h lh s]=groupplot(Re,[a.eu1],Pe);
  legend(lh,cellsprintf('%.0f',s));
  xlabel('Re');
  ylabel('U1');


