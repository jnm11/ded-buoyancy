%rsync -vap --progress hamilton:gc/f7/a/2000 ~/gc/f7/a  --exclude "*old*" --exclude b --delete
%rsync -vap --progress hamilton:gc/f6/mg/2000/0400/4608/29210 ~/gc/f6/mg/2000/0400/4608  --exclude "final" --exclude b --exclude ay --delete


nm='gc/f7/a/2000/2000/2304/016/xxxxa';
nm='gc/f7/a/2000/2000/2304/016/22637';
nm='gc/f7/a/2000/2000/2304/016/xxxxx';
nm='gc/f7/a/2000/2000/2304/016/xxxxb';
nm='gc/f7/a/2000/2000/2304/016/xxxxb';
nm='gc/f7/a/2000/2000/4608/001/xxxxb';
nm='gc/f7/a/2000/2000/4608/032/22200';
%nm='gc/f6/mg/2000/0400/4608/29210';
T=30;

p=ded_read_param(nm);
a=ded_read_javrg(nm,'ay',[T inf],'false',[],1/p.W);
% $$$ fns=ded_get_fn(nm,'final',[],'state');a=ded_read_hdf(fns{end});
% $$$ a=ded_meany(a);
fx=ded_read_hdf([nm '/force/fx.hdf5']);

%pp.W=2048;pp.ar=2;ded_gc_2db(nm,{'b','psi'},'ay',pp);
%pp.W=2048;pp.ar=2;ded_gc_2db(nm,{'b'},'b',pp);
c=ded_coord(nm);
if ~isfield(p,'wu');p.wu=0.4;end
if ~isfield(p,'hu');p.hu=0.8;end
if ~isfield(p,'U1');p.U1=0.5;end
h=p.hu;
w=p.wu;

minx=(c.Ax(fx.wdls(end)+2));
maxx=max(c.x(max(a.b>0.1,[],1)))/2;
f=find(c.x>=minx & c.x<=maxx);
x=c.x(f);
% $$$ u=eerf(c.z,p.hu,p.wu,p.U1,p.V,p.H);
% $$$ mu=mean(a.u(:,10:20),2);
% $$$ plot(c.z,u,c.z,mu);

au=a.u;
%au=squeeze(mean(a.u,2));
n=length(f);
hu=zeros(n,1);wu=zeros(n,1);U1=zeros(n,1);
for j=1:length(f)
  ffx= @(q) eerf(c.z,q(1),q(2),q(3),p.V,p.H)-au(:,f(j));
  q=lsqnonlin(ffx,[p.hu,p.wu,p.U1],[0.2 0.1 0.1],[1 1 1]);
  hu(j)=q(1);
  wu(j)=q(2);
  U1(j)=q(3);
  u=eerf(c.z,q(1),q(2),q(3),p.V,p.H);
  plot(c.z,u,c.z,au(:,j));
  drawnow
end
maxhu=max(hu(c.x>0.2&c.x<10));
[minU1 fU1]=min(U1(c.x>0.2&c.x<10));
[minwu fwu]=min(wu(c.x>0.2&c.x<10));
maxhu=hu(round(mean(fU1+fwu)));
subplot(3,1,1);plot(x,hu,[minx maxx],maxhu([1 1]),[minx maxx],p.hu([1 1]));ylabel('hu');axis('tight');
subplot(3,1,2);plot(x,wu,[minx maxx],minwu([1 1]),[minx maxx],p.wu([1 1]));ylabel('wu');axis('tight');
subplot(3,1,3);plot(x,U1,[minx maxx],minU1([1 1]),[minx maxx],p.U1([1 1]));ylabel('U1');axis('tight');
disp([maxhu minU1 minwu]);disp([p.hu p.U1 p.wu])



  
% $$$ hu=0.8526;
% $$$ wu=0.3705;    
% $$$ U1=0.5044;





%ded_gc_test_f7(nm,[20 100]);
%fs=ded_gc_test_f7_2(nm);


% $$$ 2021-01-11
% $$$ 
% $$$ Looking at this suggests U1=0.36 would be much better
% $$$ plot(fs.u(1,:))
% $$$ started gc/f7/a/2000/2000/2304/016/xxxxa
% $$$  with this value
% $$$  

 
% $$$ h=0.64;
% $$$ wu=0.5;
% $$$ U1=-0.8718;

minx=(c.Ax(fx.wdls(end)+2));
maxx=(c.Ax(fx.wdrs(1)-2));
f=find(c.x>=minx & c.x<=maxx);
x=c.x(f);

figure;
imagesc(a.SB(:,f));
title('SB');
figure;
imagesc(a.SR(:,f));
title('SR');

Iz=ichebintw(c.Nz,p.H);
Izz=ichebintw(c.NAz,p.H)/p.H;
[wb wt]=ichebendsw(c.Nz);
I0=(wt-wb);
I1=Iz;
I2=(p.H-c.z').*I1;


parity=ded_read_parity(nm);
c.dim=2;
c.c='xz';
d1=ded_diff(a,p,parity,c,{'u','w','b','uw','uu','ww','p'});
d2=ded_diff(d1,p,parity,c,{'dudx','dwdx','dudz','dwdz'});
Ia=ded_int(a,p,parity,c,{'uw','u'});
parity.duwdz=[-1 1 0];
parity.ddudzdz=[-1 1 0];
Ib=ded_int(d1,p,parity,c,{'duwdz'});
Ic=ded_int(d2,p,parity,c,{'ddudzdz'});

s=ded_read_g(nm,'sev');
av=ded_read_g(nm,'avar');
ww=sqrt(Iz*a.ww);
dwdz=sqrt(Iz*d1.dwdz.^2);
figure
if isfield(av,'ww')
  subplot(4,1,1);plot(av.t,sqrt(Izz*av.ww),[a.t1 a.t2],ww([1 1]));ylabel('ww');xlabel('t');axis([-inf inf 0 inf]);
  subplot(4,1,2);plot(av.t,sqrt(Izz*av.dwdz.^2),[a.t1 a.t2],dwdz([1 1]));ylabel('dwdz');xlabel('t');axis([-inf inf 0 inf]);
end
subplot(4,1,3);plot(c.x(1:50),ww(1:50));xlabel('x');ylabel('ww');axis([-inf inf 0 inf]);
subplot(4,1,4);plot(c.x(1:50),dwdz(1:50));xlabel('x');ylabel('dwdz');axis([-inf inf 0 inf]);


%a.dudz=ichebdifff(a.u,1);
x=c.x(f);
D1=Iz*a.SR(:,f);
D2=Iz*((d1.dudz(:,f)+d1.dwdx(:,f)).^2/2+d1.dudx(:,f).^2+d1.dwdz(:,f).^2);
D3=(a.u(end,f)-a.u(1,f)).^2/p.H;
B1=Iz*a.SB(:,f);
B2=Iz*(d1.dbdx(:,f).^2+d1.dbdz(:,f).^2);
B3=(a.b(end,f)-a.b(1,f)).^2/p.H;
figure;
subplot(2,1,1);plot(x,D1,x,D2,x,D3);
ylabel('SR');
subplot(2,1,2);plot(x,B1,x,B2,x,B3);axis([-inf inf 0 inf])
ylabel('SB');
xlabel('x');

figure
plot(x,Iz*a.bu(:,f));
xlabel('x');
ylabel('Ibu');

%plot(x,Iz*a.p(:,f)+Iz*a.uu(:,f));

% (ww)_z + (uw)_x + p_z + g*b -(w_xx+w_zz)/Re 
% Iww0   + Iuwx1  + Ip0 + gIb1 -(Iwxx1-Iwz0)/Re
% Iww1   + Iuwx2  + Ip1 + gIb2 -(Iwxx2)/Re

% (uu)_x + (uw)_z + p_x  - (u_xx+u_zz)/Re 
% Iuux1  +        + Ipx1 - (Iuxx1)/Re 
% Iuux2  + Iuw1   + Ipx2 - (Iuxx2+Iu0)/Re 

% uu    + (uwIx)_z + p   - (u_x+u_zzIx)/Re 
% Iuu0  +          + Ip0 - (Iux0+u_zzIx)
% Iuu1  +          + Ip1 - (Iux1)/Re 
% Iuu2  + Iuw1Ix   + Ip2 - (Iux2+Iu0Ix)/Re 
ddudzdzIx=Ic.ddudzdzIx(:,f)-Ic.ddudzdzIx(:,f(end));

ExI  = a.uu   + Ib.duwdzIx + a.p - (d1.dudx+Ic.ddudzdzIx)/p.Re; ExI =ExI(:,f)-ExI(:,f(end));
ExP1 = a.uu/2              + a.p;                               ExP1=ExP1(:,f)-ExP1(:,f(end));
ExP2 = a.uu   + Ib.duwdzIx + a.p;                               ExP2=ExP2(:,f)-ExP2(:,f(end));
ExP3 = a.uu/2              + a.p - (d1.dudx+Ic.ddudzdzIx)/p.Re; ExP3=ExP3(:,f)-ExP3(:,f(end));
figure;
subplot(2,1,1);plot(x,wt*ExI,x,wt*ExP1,x,wt*ExP2,x,wt*ExP3);ylabel('P top');
subplot(2,1,2);plot(x,wb*ExI,x,wb*ExP1,x,wb*ExP2,x,wb*ExP3);ylabel('P bot');
xlabel('x');







c1=(wb*ddudzdzIx)'\(wb*ExP2)'

% $$$ plot(x,wb*ExP1,x,wb*ExP2,x,wb*ExP2-c1*wb*ddudzdzIx)

a.p=a.p-mean(mean(a.p(:,f)));
Iu0   = I0*a.u(:,f);
Iu1   = I1*a.u(:,f);
Iu2   = I2*a.u(:,f);
Ib0   = I0*a.b(:,f);
Ib1   = I1*a.b(:,f);
Ib2   = I2*a.b(:,f);
Ip0   = I0*a.p(:,f);
Ip1   = I1*a.p(:,f);
Ip2   = I2*a.p(:,f);
Iww0  = I0*a.ww(:,f);
Iww1  = I1*a.ww(:,f);
Iww2  = I2*a.ww(:,f);
Iuu0  = I0*a.uu(:,f);
Iuu1  = I1*a.uu(:,f);
Iuu2  = I2*a.uu(:,f);
%Iuw0  = I0*a.uw(:,f);  identically zero
Iuw1  = I1*a.uw(:,f);
Iuw2  = I2*a.uw(:,f);
Iuwx0 = I0*d1.duwdx(:,f);
Iuwx1 = I1*d1.duwdx(:,f);
Iuwx2 = I2*d1.duwdx(:,f);
Iuux0 = I0*d1.duudx(:,f);
Iuux1 = I1*d1.duudx(:,f);
Iuux2 = I2*d1.duudx(:,f);
Iwxx0 = I0*d2.ddwdxdx(:,f);
Iwxx1 = I1*d2.ddwdxdx(:,f);
Iwxx2 = I2*d2.ddwdxdx(:,f);
Iuxx0 = I0*d2.ddudxdx(:,f);
Iuxx1 = I1*d2.ddudxdx(:,f);
Iuxx2 = I2*d2.ddudxdx(:,f);
Iux0  = I0*d1.dudx(:,f);
%Iux1  = I1*d1.dudx(:,f);  % Identically zero
Iux2  = I2*d1.dudx(:,f);
Iwz0  = I0*d1.dwdz(:,f);
Ipz0  = I0*d1.dpdz(:,f);
Iwzz0 = I0*d2.ddwdzdz(:,f);
Iwwz0 = I0*d1.dwwdz(:,f);
Iw0   = I0*a.w(:,f);
Iw1   = I1*a.w(:,f);
Iw2   = I2*a.w(:,f);
Pt    = wt*a.p(:,f);
Pb    = wb*a.p(:,f);

% Iuw0Ix = I0*Ia.uwIx(:,f); identically 0
Iuw1Ix = I1*Ia.uwIx(:,f);    
Iuw2Ix = I2*Ia.uwIx(:,f);
Iu0Ix  = I0*Ia.uIx(:,f);    
Iu1Ix  = I1*Ia.uIx(:,f);    
Iu2Ix  = I2*Ia.uIx(:,f);
Iduwdz0Ix=I0*Ib.duwdzIx(:,f);
Iddudzdz0Ix=I0*Ic.ddudzdzIx(:,f);
Ez   = +d1.duwdx +d1.dwwdz +d1.dpdz +p.g*a.b -(d2.ddwdxdx +d2.ddwdzdz)/p.Re;Ez=Ez(:,f);
rz0 = +Iuwx0    +Iwwz0    +Ipz0    +p.g*Ib0 -(Iwxx0      +Iwzz0     )/p.Re;
rz1 = +Iuwx1    +Iww0     +Ip0     +p.g*Ib1 -(Iwxx1      +Iwz0      )/p.Re;
rz2 = +Iuwx2    +Iww1     +Ip1     +p.g*Ib2 -(Iwxx2      +Iw0       )/p.Re     -p.H*Pb;

rz0 =                     +Ipz0    +p.g*Ib0 -(Iwzz0      )/p.Re;
rz1 = +Iuwx1              +Ip0     +p.g*Ib1 -(Iwxx1 +Iwz0)/p.Re;
rz2 = +Iuwx2    +Iww1     +Ip1     +p.g*Ib2 -(Iwxx2      )/p.Re    -p.H*Pb;

rx0 = Iuu0  + Iduwdz0Ix + Ip0 - (Iux0+Iddudzdz0Ix)/p.Re;
rx1 = Iuu1  +           + Ip1;
rx2 = Iuu2  + Iuw1Ix    + Ip2 - (Iux2+Iu0Ix)/p.Re; 

rz0 = +Ipz0    +p.g*Ib0;
rz1 = +Iuwx1              +Ip0     +p.g*Ib1;
rz2 = +Iuwx2    +Iww1     +Ip1     +p.g*Ib2 -p.H*Pb;

% $$$ e1=[rp1/Ib1 rp1/Ip0 rp1/Iuwx1 rp1/Iwxx1 rp1/Iwz0];
% $$$ e2=[rp2/Ib2 rp2/Ip1 rp2/Iuwx2 rp2/Iwxx2];
figure;
subplot(5,1,1);plot(x,rz1);ylabel('rz1');
subplot(5,1,2);plot(x,rz2);ylabel('rz2');
subplot(5,1,3);plot(x,rx0-rx0(end));ylabel('rx0');
subplot(5,1,4);plot(x,rx1-rx1(end));ylabel('rx1');
subplot(5,1,5);plot(x,rx2-rx2(end));ylabel('rx2');
xlabel('x');

% $$$ plot(I0*ExI-rx0)
% $$$ plot(I1*ExI-rx1)
% $$$ plot(I2*ExI-rx2)




if isfield(s,'divz')
  clear hu wu U1
  for j=1:size(s.divz,2)
    ffx= @(q) eerf(c.Az,q(1),q(2),q(3),p.V,p.H)-s.divz(:,j);
    q=lsqnonlin(ffx,[p.hu,p.wu,p.U1],[0.2 0.1 0.1],[2 2 2]);
    hu(j)=q(1);
    wu(j)=q(2);
    U1(j)=q(3);
    u=eerf(c.Az,q(1),q(2),q(3),p.V,p.H);
    plot(c.Az,u,c.Az,s.divz(:,j));
    drawnow
  end
  figure;clf;
  subplot(3,1,1);plot(s.t,hu,s.t,s.hu);ylabel('hu');axis('tight');
  subplot(3,1,2);plot(s.t,wu,s.t,s.wu);ylabel('wu');axis('tight');
  subplot(3,1,3);plot(s.t,U1,s.t,s.U1);ylabel('U1');axis('tight');
  xlabel('x');
% $$$   f=find(s.t>90 & s.t<inf);
% $$$   ffx= @(q) eerf(c.Az,q(1),q(2),q(3),p.V,p.H)-mean(s.divz(:,f),2);
% $$$   q=lsqnonlin(ffx,[p.hu,p.wu,p.U1],[0.2 0.1 0.1],[2 2 2]);
end

figure;
clf;
ded_plot_X(nm);

