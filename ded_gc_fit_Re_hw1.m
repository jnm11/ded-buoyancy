function a=ded_gc_fit_Re_hw1(Re,uh,uw,u1,g,nms,typ,Pe)

a.RS=max(Re);
a.PS=max(Pe);
P=Pe/a.PS;
R=Re/a.RS;
nPe=length(unique(Pe));

opt=optimset('display','none');

%fuh = @(p,R) polyval(p(1:3),R)./(1+p(4)*(R+p(5)).^2);
%fi =@(p) fuh(p,R)-uh;
[fuh a.puh] = fit_generic(R,uh,typ{1});
 
fs = @(p,R) (1+erf(p(1)*(R-p(2))))/2;
fwa = @(p,R) max(p(3),p(1)./max(eps,R.^p(2)));
fwb = @(p,R) p(1)-p(2)*exp(-p(3)*R);

if nPe==1
  [fuw a.puw]=fit_generic(R,uw,typ{2});
else
  switch(typ{2})
    case 'poly3'
      fuw = @(p,R,P) polyval(p,R);
      p0=[0 0 mean(R)];
    case 'low'
      fuw = @(p,R,P) p(1)./max(eps,R.^p(2));
      p0=[0.3174    0.2377];
    case 'mid'
      fuw = @(p,R,P) max(p(3),p(4)+(p(1)+p(5)*P.^p(6))./max(eps,R.^p(2)));
      p0=[0.3174    0.2377 0.4 0 0 1];
    case 'high'
      fuw = @(p,R,P)fwb(p,R);
      p0=[1    0.5   1];
    case 'all'
      fuw = @(p,R,P) fwa(p(1:3),R)+fs(p(7:8),R).*(fwb(p(4:6),R)-fwa(p(1:3),R));
      p0=[0.3174    0.2377    0.3658    1    0.5   1  10.0745    0.6256];
  end
  fi =@(p) fuw(p,R,P)-uw;
  %puw = lsqnonlin(fi,[1 1 1 0.4 0],[],[],opt);
  lb=repmat(-inf,size(p0));
  ub=repmat( inf,size(p0));
  a.puw = lsqnonlin(fi,p0,lb,ub,opt);
end


if nPe==1
  [fg a.pg]=fit_generic(R,g,typ{4});
else
  fg = @(p,R) polyval(p(1:3),R)./polyval([p(4:5) 1],R);
  fi =@(p) fg(p,R)-g;
  a.pg = lsqnonlin(fi,[1.4197    1.2041    0.4363    0.6308    0.6806]*100,[],[],opt);
end

fu1a = @(p,R) polyval(p(1:3),R);
fu1b = @(p,R) p(1)-p(2)*exp(-p(3)*R);
rg1=1:3;rg2=4:6;rgs=7:8;
if nPe==1
  [fu1 a.pu1]=fit_generic(R,u1,typ{3});
else
  switch(typ{3})
    case 'poly3'
      fu1 = @(p,R,P) polyval(p,R)
      p0=[0 0 mean(R)];
    case({'low'})
      p0=[0 0 mean(u1)];
      fu1  = @(p,R,P) fu1a(p,R);
    case({'mid'})
      p0=[0 0 mean(u1) 0 0 0];
      fu1  = @(p,R,P) polyval(p(1:3),R)+P.*polyval(p(4:6),R);
    case 'high'
      p0=[1 0.5 0];
      fu1  = @(p,R,P) fu1b(p,R);
    case 'all'
      p0=[0 0 mean(u1) 1 0.5 0 10 0.5];
    fu1  = @(p,R,P) fu1a(p(rg1),R) + fs(p(rgs),R).*(fu1b(p(rg2),R)-fu1a(p(rg1),R));
  end
  fi =@(p) fu1(p,R,P)-u1;
  a.pu1 = lsqnonlin(fi,[0 0 mean(u1) 1 0.5 0 10 0.5],[],[],opt);
end

  
  

RR =linspace(min(R),max(R),1e3);
uP=unique(P(:));



figure(1);clf;
subplot(4,1,1);groupplot(Re,uh,Pe);hold('on');plot(a.RS*RR,fuh(RR));ylabel('uh');
subplot(4,1,2);groupplot(Re,uw,Pe);hold('on');plot(a.RS*RR,fuw(RR));ylabel('uw');
subplot(4,1,3);groupplot(Re,u1,Pe);hold('on');plot(a.RS*RR,fu1(RR));ylabel('u1');
subplot(4,1,4);groupplot(Re, g,Pe);hold('on');plot(a.RS*RR, fg(RR));ylabel( 'g');

figure(2);clf;
subplot(4,1,1);groupplot(Re,uh-fuh(Re/a.RS),Pe);ylabel('duh');
subplot(4,1,2);groupplot(Re,uw-fuw(Re/a.RS),Pe);ylabel('duw');
subplot(4,1,3);groupplot(Re,u1-fu1(Re/a.RS),Pe);ylabel('du1');
subplot(4,1,4);groupplot(Re, g- fg(Re/a.RS),Pe);ylabel('dg');


cuw = fuw(R);
cu1 = fu1(R);
cuh = fuh(R);
cg  = fg(R);
for j=1:length(R)
  disp(sprintf('--wu %6.4f --U1 %6.4f --hu %6.4f --g %6.4f %s/%5.0f',cuw(j),cu1(j),cuh(j),cg(j),nms{j},cg(j)*1e4));
end



a.nm  = nms; 
a.uw  = cuw;
a.u1  = cu1;
a.g   = cg;
a.fuw = @(Pe) fuw(Re/a.RS);
a.fu1 = @(Pe) fu1(Re/a.RS);
a.fg  = @(Pe) fg(Re/a.RS);

