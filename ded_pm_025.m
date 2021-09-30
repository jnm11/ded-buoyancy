function  [b c p]=ded_pm_025(nm,T,display)
%nm='pm/025';

if nargin<3
  display=nargout==0;
end

a=ded_read_g(nm,'ayz',[],[T inf]);
p=ded_read_param(nm);

if ~isfield(a,'s')
  a.s =NaN*a.b;
  a.su=NaN*a.b;
  a.ss=NaN*a.b;
  p.S=NaN;
  p.SV=NaN;
end

tol=1e-3;

A=pi*p.radius^2;
b.A=A;
dt=sum(a.dt);
b.u  = sum(a.u,2) /(dt*A*p.U);
b.s  = sum(a.s,2) /(dt*A*p.S);
b.b  = sum(a.b,2) /(dt*A*p.B);
b.su = sum(a.su+p.SV*a.s,2)/(dt*A*p.S*p.U);
b.bu = sum(a.bu,2)/(dt*A*p.B*p.U);
b.ss = sum(a.ss,2)/(dt*A*p.S^2);
b.bb = sum(a.bb,2)/(dt*A*p.B^2);
b.uu = sum(a.uu,2)/(dt*A*p.U^2);     
b.sw=sqrt(b.s.^2./max(tol,b.ss));
b.bw=sqrt(b.b.^2./max(tol,b.bb));
b.uw=sqrt(b.u.^2./max(tol,b.uu));  
b.x=a.x;
b.dt=dt;
b.t1=a(1).t1;
b.t2=a(end).t2;
b.A  = A;
b.B  = p.B;
b.S  = p.S;
b.U  = p.U;
b.nm = nm;

% This is not effective since u can be negative

dt=repmat(a.dt(:)',[length(a.x) 1]);
c.u  = a.u  ./(dt*A*p.U);        
c.s  = a.s  ./(dt*A*p.S);        
c.b  = a.b  ./(dt*A*p.B);        
c.su = (a.su+p.SV*a.s) ./(dt*A*p.S*p.U);    
c.bu = a.bu ./(dt*A*p.B*p.U);    
c.ss = a.ss ./(dt*A*p.S^2);      
c.bb = a.bb ./(dt*A*p.B^2);      
c.uu = a.uu ./(dt*A*p.U^2);       
c.sw = sqrt(c.s.^2./max(tol,c.ss));
c.bw = sqrt(c.b.^2./max(tol,c.bb));
c.uw = sqrt(c.u.^2./max(tol,c.uu));
c.x  = a.x;
c.dt = a.dt;
c.t1 = a.t1;
c.t2 = a.t2;
c.t  = a.t;
c.A  = A;
c.B  = p.B;
c.S  = p.S;
c.U  = p.U;
c.nm = nm;

nt=length(a.t);

f=find(b.x>p.x6);
maxb=max(b.b(f));
maxs=max(b.s(f));
maxbu=max(b.bu(f));
maxsu=max(b.su(f));
maxbw=max(b.bw(f));
maxsw=max(b.sw(f));
maxuu=max(max(c.uu(f,:)));
maxu=max(max(c.u(f,:)));
maxbs=1.2*max(maxb,maxs);
maxbsu=1.2*max(maxbu,maxsu);
maxbsw=1.2*max(maxbw,maxsw);

if display >0
  clf;
  for j=1:nt;
    
    subplot(5,1,1);
    plot(a.x,c.s(:,j),a.x,c.b(:,j),c.x,b.s,c.x,b.b);
    line([p.x6 p.x6],[0 maxbs]);line([p.x7 p.x7],[0 maxbs]);
    ylabel('b s');
    title(sprintf('%s %6.3f',nm,c.t(j)));
    axis([0 p.L 0 maxbs]);
    
    subplot(5,1,2);
    plot(a.x,c.su(:,j),c.x,c.bu(:,j),c.x,b.su,c.x,b.bu);
    line([p.x2 p.x2],[0 maxbsu]);line([p.x3 p.x3],[0 maxbsu]);
    ylabel('bu su');
    axis([0 p.L 0 maxbsu]);
    
    subplot(5,1,3);
    plot(a.x,c.sw(:,j),c.x,c.bw(:,j),c.x,b.sw,c.x,b.bw);
    line([p.x6 p.x6],[0 maxbsw]);line([p.x7 p.x7],[0 maxbsw]);
    ylabel('r, sb');
    axis([0 p.L 0 maxbsw]);
    
% $$$   subplot(5,1,4);
% $$$   plot(a.x,c.su./max(tol,c.s),c.x,c.bu./max(tol,c.b),c.x,b.su./max(tol,b.s),c.x,b.bu./max(tol,b.b),c.x,b.uu./max(tol,b.u));
% $$$   ylabel('u, sbu');
% $$$   axis([0 30 0 5]);
% $$$   line([p.x2 p.x2],[0 5]);line([p.x3 p.x3],[0 5]);
    
    subplot(5,1,4);
    plot(c.x,c.u(:,j),b.x,b.u);
    ylabel('u');
    axis([0 p.L 0 maxu]);
    line([p.x2 p.x2],[0 maxu]);line([p.x3 p.x3],[0 maxu]);
 
    subplot(5,1,5);
    plot(c.x,sqrt(c.uu(:,j)),b.x,sqrt(b.uu));
    ylabel('uu');
    axis([0 p.L 0 maxuu]);
    line([p.x2 p.x2],[0 maxuu]);line([p.x3 p.x3],[0 maxuu]);
    
    
    drawnow;
  end
end


return;
[b1 c1 p1]=ded_pm_025('pm/030',[10 inf], 1);
[b1 c1 p1]=ded_pm_025('pm/031',[0 inf], 1);


[b1 c1 p1]=ded_pm_025('pm/033',[0 inf], 1);
[b2 c2 p2]=ded_pm_025('pm/025',[10 inf], 1);

[b1 c1 p1]=ded_pm_025('pm/024',[0 inf], 1);
[b2 c2 p2]=ded_pm_025('pm/025',[10 inf], 1);


rg=find(b1.x>=p1.x7 & b1.x<=p1.L-1);
ff=@(q,x)  q(2)*max(0,x-q(1)-p1.x7) + interp1([0 p1.x5+p1.Wx*[-0.5 0.5] p1.x6 p1.x7+(0:5)/4*q(end) p1.L],[0 0 1 1 1 q(3:6) 0 0],x,'pchip');
q=[2.2483    0.3168    0.8321    0.5819    0.3913    0.1820    9.5764];
fff=@(q) ff(q,b1.x(rg))-b1.bw(rg);
q1b=lsqnonlin(fff,q,[0 0 0 0 0 0 0],[inf inf inf inf inf inf 7]);
fff=@(q) ff(q,b2.x(rg))-b2.bw(rg);
q2b=lsqnonlin(fff,q,[0 0 0 0 0 0 0],[inf inf inf inf inf inf 7]);
fff=@(q) ff(q,b2.x(rg))-b2.sw(rg);
q2s=lsqnonlin(fff,q,[0 0 0 0 0 0 0],[inf inf inf inf inf inf 7]);
preprint([6 3],12);clf;
h(1:2)=plot(b1.x,b1.bw,'b-',b1.x,ff(q1b,b1.x),'b--');hold('on');
h(3:4)=plot(b2.x,b2.bw,'r-',b2.x,ff(q2b,b2.x),'r--');
h(5:6)=plot(b2.x,b2.sw,'k-',b2.x,ff(q2s,b2.x),'k--');
title(sprintf('Entrainment coefficients [%4.2f %4.2f %4.2f]',q1b(2),q2b(2),q2s(2)));
xlabel('z');
ylabel('plume width');
legend(h,{'b1','b1 fit','b2','b2 fit','s2','s2 fit'},'location','NorthWest');
axis([2 29 0 8]);
print('-depsc2','~/b1b2s2.eps');
unix('scp ~/b1b2s2.eps asahi:Dropbox/Jim-Claudia-Plumes/b1b2s2.eps');


rg=find(x>7 & x<27);
% $$$ p1=polyfit(x(rg),b.1(rg),1);
% $$$ p2=polyfit(x(rg),b.2(rg),1);
% $$$ plot(x,b.1,'r-',x,b.2,'b-',x,polyval(p1,x),'r--',x,polyval(p2,x),'b--');


