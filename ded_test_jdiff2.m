function [rmse maxe]=ded_test_jdiff2

nm='gc/diff/01';
p=ded_read_param(nm);
c=ded_coord(nm);
pr=ded_read_hdf([ded_dedalus_data_dir '/' nm '/parity.hdf5']);
fns=cellstr_ls([ded_dedalus_data_dir '/' nm '/diff/*']);
A=struct();J=struct();dn={};
cc=flip(c.c)';
a=repmat(struct(),1,2);
for j=1:length(fns)
  [dd nn]=fileparts(fns{j});
  b=ded_read_hdf(fns{j});
  m=fieldnames(b);
  if nn(1)=='A'  
    a(1).(nn(2:end))=b.(m{1});
  else          
    a(2).(nn(2:end))=b.(m{1});
  end
  rg=min(find(nn=='d'))+1:max(find(nn=='d'))-1;
  if ~isempty(rg)
    dn{end+1}=nn(rg);
  end
end
a(1).x=c.Ax;

a(1).z=c.Az;
a(1).dd=flip(c.ddA);
a(1).dd=flip(c.ddA);
a(2).dd=flip(c.ddJ);
a(1).X=flip(p.X);
a(2).X=flip(p.X);
a(2).x=c.Jx;
a(2).z=c.Jz;
dn=unique(dn);
rmse=size([2,length(dn),c.dim]);
maxe=rmse;
AA='AJ';
p.T=flip(p.T);
for i=1:2
  for j=1:length(dn)
    pp=flip(pr.(dn{j}));
    figure;clf
    for k=1:c.dim
      switch(p.T{k})
        case 'Cheb'
          f1=ichebdifff(a(i).(dn{j}),k,a(i).X(k)); % Assume spacing is Chebychev points
        case {'SinCos','Fourier'}
          f1=pr_diff(a(i).(dn{j}),a(i).dd(k),k,1,[],pp(k));
        otherwise
          error(sprintf('Unknown coordinates %s',p.T{k}));
      end
      f2=a(i).(['d' dn{j} 'd' cc(k)]);
      rmse(i,j,k)=sqrt(mean((f1(:)-f2(:)).^2));
      maxe(i,j,k)=max(abs(f1(:)-f2(:)));
      switch(k)
        case(1)
          subplot(4,1,3);plot(a(i).z,f1(:,20),a(i).z,f2(:,20));
          title(sprintf('%s d%sd%s parity %d rms %8.1e, max %8.1e',AA(i),dn{j},cc(k),pp(k),rmse(i,j,k),maxe(i,j,k)));
          subplot(4,1,4);plot(a(i).z,f1(:,20)-f2(:,20));
          xlabel(cc(k));
        case(2)
          subplot(4,1,1);plot(a(i).x,f1(20,:),a(i).x,f2(20,:));
          title(sprintf('%s d%sd%s parity %d rms %8.1e, max %8.1e',AA(i),dn{j},cc(k),pp(k),rmse(i,j,k),maxe(i,j,k)));
          subplot(4,1,2);plot(a(i).x,f1(20,:)-f2(20,:));
          xlabel(cc(k));
      end
      drawnow;
    end
  end
end


keyboard;
w1=a(1).w;
w2=a(2).w;

Aw=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Aw.hdf5']);
Jw=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Jw.hdf5']);

w3=icheb_resize(w2,1,size(w1,1));
Jw.w=icheb_resize(Jw.w,1,size(w1,1));
w3=fftresample(w3,size(w1),flip(pr.w));
Jw.w=fftresample(Jw.w,size(w1),flip(pr.w));

figure(1);clf;subplot(2,1,1);plot(c.Ax,w1(20,:),c.Ax,Aw.w(20,:));subplot(2,1,2);plot(c.Ax,w1(20,:)-Aw.w(20,:));
figure(2);clf;subplot(2,1,1);plot(c.Ax,w3(20,:),c.Ax,Jw.w(20,:));subplot(2,1,2);plot(c.Ax,w3(20,:)-Jw.w(20,:));
figure(3);clf;
subplot(3,1,1);plot(c.Ax,w1(20,:)-w3(20,:));
subplot(3,1,2);plot(c.Ax,w1(20,:)-Aw.w(20,:));
subplot(3,1,3);plot(c.Ax,w1(20,:)-Jw.w(20,:));



Aw=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Adwdx.hdf5']);
Jw=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Jdwdx.hdf5']);
w1=a(1).dwdx;
w2=a(2).dwdx;
figure(1);clf;subplot(2,1,1);plot(c.Jx,w2(20,:),c.Jx,Jw.dwdx(20,:));subplot(2,1,2);plot(c.Jx,w2(20,:)-Jw.dwdx(20,:));

%These should agree


w3=icheb_resize(w2,1,size(w1,1));
Jw.dwdx=icheb_resize(Jw.dwdx,1,size(w1,1));
w3=fftresample(w3,size(w1),flip(pr.w));
Jw.dwdx=fftresample(Jw.dwdx,size(w1),flip(pr.w));

figure(2);clf;subplot(2,1,1);plot(c.Ax,w1(20,:),c.Ax,Aw.dwdx(20,:));subplot(2,1,2);plot(c.Ax,w1(20,:)-Aw.dwdx(20,:));
figure(3);clf;subplot(2,1,1);plot(c.Ax,w3(20,:),c.Ax,Jw.dwdx(20,:));subplot(2,1,2);plot(c.Ax,w3(20,:)-Jw.dwdx(20,:));
figure(4);clf;
subplot(4,1,1);plot(c.Ax,w1(20,:)-w3(20,:));
subplot(4,1,2);plot(c.Ax,w1(20,:)-Aw.dwdx(20,:));
subplot(4,1,3);plot(c.Ax,w1(20,:)-Jw.dwdx(20,:));
subplot(4,1,4);plot(c.Ax,w3(20,:)-Jw.dwdx(20,:));

Aw=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Adwdx.hdf5']);
J1=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Jdwdx1.hdf5']);
J2=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Jdwdx2.hdf5']);
J3=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Jdwdx3.hdf5']);
J4=ded_read_hdf([ded_dedalus_data_dir '/' nm '/Jdwdx4.hdf5']);
figure(1);clf;
subplot(4,1,1);plot(c.Jx,J1.dwdx(20,:),c.Jx,J2.dwdx(20,:),c.Jx,J3.dwdx(20,:),c.Jx,J4.dwdx(20,:));
subplot(4,1,2);plot(c.Jx,J1.dwdx(20,:)-J2.dwdx(20,:));
subplot(4,1,3);plot(c.Jx,J1.dwdx(20,:)-J3.dwdx(20,:));
subplot(4,1,4);plot(c.Jx,J1.dwdx(20,:)-J4.dwdx(20,:));

