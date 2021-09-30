function [b1 b2]=ded_test_jdiff(nm,typ,trg)
%ded_test_jdiff('pm/test/aino/01','state');
%ded_test_jdiff('gc/diff/01','state');
%ded_test_jdiff('pm/tttest','a');
if nargin<3
  trg=[];
end
if isempty(trg)
  trg=[0 inf];
end

p=ded_read_param(nm);

c=ded_coord(nm);
if strcmp(typ,'state')
  fns=sort(ded_get_fn(nm,'final',[],'state'));
else
  a=ded_read_javrg(nm,typ,trg,'combine');
end
if isempty(a)
  disp(sprintf('ded_test_jdiff: no files of type %s in range',typ));
  return;
end
if ~isfield(a,'d')
  fdfile=[ded_dedalus_data_dir '/' nm '/force/fd.hdf5'];
  if isfile(fdfile)
    fd=ded_read_hdf(fdfile);
    a.d=fd.fd;
    clear('fd');
  else
    a.d=0;
  end
end

parity=ded_read_hdf([ded_dedalus_data_dir '/' nm '/' 'parity.hdf5']);

f=ded_hdf(nm,'fx');
nm1={'d','s','b','p','u','v','w'};
nms={};
for j=1:length(nm1)
  for k=1:length(nm1)
    nms{j,k}=[nm1{j} nm1{k}];
  end
end
nms={nms{:},nm1{:}};

nms=intersect(fieldnames(a),nms);
da=ded_diff(a,p,parity,c,nms);
Ia=ded_int(a,p,parity,c,nms);

if c.dim==3
  da.d=da.dudx+da.dvdy+da.dwdz;
  if isfield(a,'ududx') & isfield(a,'vdudy') & isfield(a,'wdudz')
    a.cdu  =  a.ududx  +a.vdudy  +a.wdudz +a.ud;
    a.cdv  =  a.udvdx  +a.vdvdy  +a.wdvdz +a.vd;
    a.cdw  =  a.udwdx  +a.vdwdy  +a.wdwdz +a.wd;
  end
else
  da.d=da.dudx+da.dwdz;
  if isfield(a,'ududx') & isfield(a,'wdudz')
    a.cdu  =  a.ududx   +a.wdudz +a.ud;
    a.cdw  =  a.udwdx   +a.wdwdz +a.wd;
  end
end


cc='xyz';
for j=1:length(nm1)
  for k=1:3
    w1=[nm1{j} 'd' nm1{j} 'd' cc(k)];
    w2=['d' nm1{j} nm1{j} 'd' cc(k)];
    if isfield(a,w1) & isfield(da,w2)
      da.(w1)=da.(w2)/2;
    end
  end
end

nms=intersect(fieldnames(a),fieldnames(da));

n=length(nms);
m=4;
n1=ceil(n/m);
j=0;
for j1=1:n1
  figure; ahx=jsubplot([1 m],[0.15 0.05],[0.02 0.02],[0.01 0.01]);
  figure; ahz=jsubplot([m 1],[0.05 0.15],[0.02 0.02],[0.01 0.01]);
  for j2=1:m;
    j=j+1;
    if j>n
      break;
    end
    e1= a.(nms{j});
    e2=da.(nms{j});  
    f=(e1-e2)./(1+abs(e1)+abs(e2));
    maxf=max(abs(f(:)));
    rmsf=sqrt(mean(f(:).^2));
    s1=sprintf('%8s %5.0e %5.0e',nms{j},rmsf,maxf);
    s=sprintf('%8s\n %5.0e %5.0e',nms{j},rmsf,maxf);
    disp(s1)
    
    if ndims(f)==3
      f=squeeze(sqrt(mean(f.^2,2)));
    end
    axes(ahx(j2));
    plot(c.Jx,f);axis('tight');
    ylabel(s);
    if j2==m
      xlabel('x');
    end
    
    axes(ahz(j2));
    plot(f,c.Jz);axis('tight');
    xlabel(s);
    if j2==1
      zlabel('z');
    end
  end
  set(ahx(1:end-1),'xticklabel',[]);
  set(ahz(2:end  ),'yticklabel',[]);
  set(ahx,'yticklabel',[]);
  set(ahz,'xticklabel',[]);
end

return;

%dwdz  = ichebdifff(a.w,dimz,p.H,1,'Dirichlet','Dirichlet',c.Jz);

a.duudx = pr_diff(a.uu,dx,dimx,1,[], 1);
a.duwdx = pr_diff(a.uw,dx,dimx,1,[],-1);
a.dwwdx = pr_diff(a.ww,dx,dimx,1,[], 1);
a.duudz = ichebdifff(a.uu,dimz,p.H);
a.duwdz = ichebdifff(a.uw,dimz,p.H);
a.dwwdz = ichebdifff(a.ww,dimz,p.H);


W=f.wux4/max(f.wux4)+f.wux5/max(f.wux5);
Wtol=1e-3;
Wx1=min(c.Ax(W<Wtol));
Wx2=max(c.Ax(W<Wtol));
JW=interp1([0 Wx1 2*Wx1 2*Wx2-p.L Wx2 p.L],[0 0 1 1 0 0],x','pchip');

%a.ududz = a.u.*a.dudz;
%a.udwdz = a.u.*a.dwdz;
%a.wdwdx = a.w.*a.dwdx;
%a.wdudx = a.w.*a.dudx;

figure;
subplot(3,2,1);plot(x,a.duudx-      2*a.ududx); ylabel('duudx');
subplot(3,2,3);plot(x,a.duwdx-a.wdudx-a.udwdx); ylabel('duwdx');
subplot(3,2,5);plot(x,a.dwwdx-      2*a.wdwdx); ylabel('dwwdx');
subplot(3,2,2);plot(x,JW.*(a.duudx-      2*a.ududx)); ylabel('duudx');
subplot(3,2,4);plot(x,JW.*(a.duwdx-a.wdudx-a.udwdx)); ylabel('duwdx');
subplot(3,2,6);plot(x,JW.*(a.dwwdx-      2*a.wdwdx)); ylabel('dwwdx');

figure;
subplot(3,2,1);plot(x,a.duudz-      2*a.ududz); ylabel('duudz');
subplot(3,2,3);plot(x,a.duwdz-a.wdudz-a.udwdz); ylabel('duwdz');
subplot(3,2,5);plot(x,a.dwwdz-      2*a.wdwdz); ylabel('dwwdz');
subplot(3,2,2);plot(x,JW.*(a.duudz-      2*a.ududz)); ylabel('duudz');
subplot(3,2,4);plot(x,JW.*(a.duwdz-a.wdudz-a.udwdz)); ylabel('duwdz');
subplot(3,2,6);plot(x,JW.*(a.dwwdz-      2*a.wdwdz)); ylabel('dwwdz');

figure;
subplot(5,1,1);plot(x,a.uu    - a.u.^2);      ylabel('uu');
subplot(5,1,2);plot(x,a.ududx - a.u.*a.dudx); ylabel('ududx');
subplot(5,1,3);plot(x,a.ududz - a.u.*a.dudz); ylabel('ududz');   
subplot(5,1,4);plot(x,a.udwdx - a.u.*a.dwdx); ylabel('udwdx');
subplot(5,1,5);plot(x,a.udwdz - a.u.*a.dwdz); ylabel('udwdz');   

figure;
subplot(5,1,1);plot(x,a.ww    - a.w.^2);      ylabel('ww');   
subplot(5,1,2);plot(x,a.wdudx - a.w.*a.dudx); ylabel('wdudx');
subplot(5,1,3);plot(x,a.wdudz - a.w.*a.dudz); ylabel('wdudz');   
subplot(5,1,4);plot(x,a.wdwdx - a.w.*a.dwdx); ylabel('wdwdx');
subplot(5,1,5);plot(x,a.wdwdz - a.w.*a.dwdz); ylabel('wdwdz');   

b=ded_zgrid(a,c.NAz,{'b','SS','divu','WX','WY','WZ','Eb'},{},{},{},p.H,dimz);
figure;
subplot(5,1,1);imagesc(c.x,b.z,b.b   );set(gca,'ydir','normal');ylabel('b');
subplot(5,1,2);imagesc(c.x,b.z,b.SS  );set(gca,'ydir','normal');ylabel('SS');
subplot(5,1,3);imagesc(c.x,b.z,b.divu);set(gca,'ydir','normal');ylabel('divu');
subplot(5,1,4);imagesc(c.x,b.z,b.WY  );set(gca,'ydir','normal');ylabel('WY');
subplot(5,1,5);imagesc(c.x,b.z,b.Eb  );set(gca,'ydir','normal');ylabel('Eb');


return;
plot([a.dwwdx(70,:) -flip(a.dwwdx(70,:)) a.dwwdx(70,:)],'s-')
plot([a.ww(70,:) flip(a.ww(70,:)) a.ww(70,:)],'s-')

plot([a.ww(70,:) flip(a.ww(70,:)) a.ww(70,:)],'s-')
axis([700 900 0 0.25])
plot([a.w(70,:) flip(a.w(70,:)) a.w(70,:)],'s-')   
axis([700 900 -inf inf])                        

a=ded_read_hdf('pm/test/aino/01/a/a-00001.hdf5');
p=ded_read_param('pm/test/aino/01');
pr=ded_read_hdf('pm/test/aino/01/parity.hdf5');
c=ded_coord('pm/test/aino/01');
dx=[c.dJx c.dJy c.dJz];
nm={'dudx','dudy','dudz','dvdx','dvdy','dvdz','dwdx','dwdy','dwdz'};

b=ded_diff(a,p,pr,c,{'u','v','w'});

nm='gc/diff/01';
p=ded_read_param(nm);
c=ded_coord(nm);
pr=ded_read_hdf([ded_dedalus_data_dir '/' nm '/parity.hdf5']);
fns=cellstr_ls([ded_dedalus_data_dir '/' nm '/diff/*']);
A=struct();J=struct();dn={},
cc=flip(c.c)';
a=repmat(struct(),1,2);
for j=1:length(fns)
  [dd nn]=fileparts(fns{j});
  b=ded_read_hdf(fns{j});
  nm=fieldnames(b);
  if nn(1)=='A'  
    a(1).(nn(2:end))=b.(nm{1});
  else          
    a(2).(nn(2:end))=b.(nm{1});
  end
  rg=min(find(nn=='d'))+1:max(find(nn=='d'))-1;
  if ~isempty(rg)
    dn{end+1}=nn(rg);
  end
end
a(1).x=c.Ax;
a(1).y=c.Ay;
a(1).z=c.Az;
a(1).dd=flip(c.ddA);
a(1).dd=flip(c.ddA);
a(2).dd=flip(c.ddJ);
a(1).X=flip(c.ddA);
a(2).X=flip(c.ddJ);
a(2).x=c.Jx;
a(2).y=c.Jy;
a(2).z=c.Jz;
dn=unique(dn);
dim=length(cn);
for i=1:2
  for j=1:length(dn)
    for k=1:dim
      switch(p.T{k})
        case 'Cheb'
          f1=ichebdifff(a(i).(dn{j}),k,a(i).X(k)); % Assume spacing is Chebychev points
        case {'SinCos','Fourier'}
          f1=pr_diff(a(i).(dn{j}),a(i).dd(k),k,1,[],pr.(dn{j})(k));
      end
      f2=a(i).(['d' dn{j} 'd' cc(k)]); 
     
      switch(k)
        case(1)
          subplot(4,1,3); plot(c.Az,f1(:,20),c.Az,f2(:,20));
          title(['d' dn{j} 'd' cc(k)]);
          subplot(4,1,4);plot(c.Az,f1(:,20)-f2(:,20));
          xlabel(cc(k));
        case(2)
          subplot(4,1,1);plot(a(i).x,f1(20,:),c.Ax,f2(20,:));
          title(['d' dn{j} 'd' cc(k)]);
          subplot(4,1,2);plot(a(i).x,f1(20,:)-f2(20,:));
          xlabel(cc(k));
      end
    end
  end
end


for j=1:length(nm)
  f=(a.(nm{j})-b.(nm{j}))./(1+abs(a.(nm{j}))+abs(b.(nm{j})));
  disp(sprintf('%s %8.1e %8.1e',nm{j},max(abs(f(:))),sqrt(mean(f(:).^2))));
end


keyboard

dudx=pr_diff(a.u,c.dJx,3,1,[],-1);
f=(a.dudx-dudx)./(realmin+abs(a.dudx)+abs(dudx));
disp(sprintf('dudx %8.1e',max(abs(f(:)))));
subplot(2,1,1);
plot(c.Jx,squeeze(a.u(20,20,:)));
subplot(2,1,2);
plot(c.Jx,squeeze(a.dudx(20,20,:))-squeeze(dudx(20,20,:)));

dudy=pr_diff(a.u,c.dJy,2,1,[],1);
f=(a.dudy-dudy)./(1+abs(a.dudy)+abs(dudy));
disp(sprintf('dudy %8.1e %8.1e',max(abs(f(:))),sqrt(mean(f(:).^2))));
subplot(2,1,1);
plot(c.Jy,squeeze(a.u(20,:,20)));
subplot(2,1,2);
plot(c.Jy,squeeze(a.dudy(20,:,20))-squeeze(dudy(20,:,20)));



a=ded_read_hdf('pm/test/aino/01/Au.hdf5');
b=ded_read_hdf('pm/test/aino/01/Adudx.hdf5');
u=squeeze(a.u(20,20,:));
dudx=squeeze(b.dudx(20,20,:));

dudx2=pr_diff(u,c.dAx,1,1,[],-1);
subplot(2,1,1);
plot(c.Ax,u);
subplot(2,1,2);
plot(c.Ax,dudx-dudx2);
size(dudx2)





a=ded_read_hdf('pm/test/aino/01/Ju.hdf5');
b=ded_read_hdf('pm/test/aino/01/Jdudx.hdf5');
u=squeeze(a.u(20,20,:));
dudx=squeeze(b.dudx(20,20,:));

dudx2=pr_diff(u,c.dJx,1,1,[],-1);
subplot(2,1,1);
plot(c.Jx,u);
subplot(2,1,2);
plot(c.Jx,dudx-dudx2);
size(dudx2)

a=ded_read_hdf('pm/test/aino/01/Ju2.hdf5');
b=ded_read_hdf('pm/test/aino/01/Jdudx6.hdf5');
u=squeeze(a.u(20,20,:));
dudx=squeeze(b.dudx(20,20,:));

dudx2=pr_diff(u,c.dJx,1,1,[],-1);
subplot(2,1,1);
plot(c.Jx,u);
subplot(2,1,2);
plot(c.Jx,dudx-dudx2);
size(dudx2)
%Tz Cheb --sType gc --lbc slip --rbc slip --WallTime 99:00:00 --scheme RK443 --blbc dz --brbc dz --Tx SinCos --Ny 1 --Nx 256 --Nz 96 --Forcing 7 --dtjsw True --fbmult True --dtjadv True --signals True --dtjd2 True --ddiv True --parallel True --db True --dtjd1 True --wwl True --Re 50.0 --Velocity 1.0 --dtjavar 1.0 --hb 1.0 --dtjcheck 30.0 --Width 1.0 --dtjsev 1.0 --dtjay 1.0 --forceb 1.0 --hu 1.0 --Length 8.0 --Scb 1.0 --x5 0.5 --dtstats 0.1 --mindt 0.0001 --Force 50.0 --dT 1.0 --AA 1.5 --dtjayz 1.0 --Height 2.0 --Time 100.0 --dtjb 0.1 --Gravity 3.5 --Angle 0.0 --maxdt 0.1 --dtja 1.0 --Buoyancy 1.0 --AAJ 1.0 --x6 0.3 --AAS 1.0 --x4 0.25 --dtj 1.0 

