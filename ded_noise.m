function a=ded_noise(n,t0)
if nargin<1
  t0=[];
end
if isempty(t0)
  t0=0;
end
p=ded_read_param(n);

a=ded_read_g(n,'noise');

a=ded_zgrid(a,2*length(a.z),[],{},{},{},p.H);

if isempty(a)
  return;
end

nmc=intersect(fieldnames(a),{'noisex','noisey','noisez'});

f=find(a.t>t0);
dims=ndims(a.noisef)-1;
switch(dims)
  case 2 
    x=a.noisef(:,:,f);
    fx=sqrt(squeeze(sum(sum(x.^2,1),3)));
    fz=sqrt(squeeze(sum(sum(x.^2,2),3)));
    at=auto_corrn(x,3);
    ax=auto_corrn(x,2);
    az=auto_corrn(x,1);
    x=a.noisef2(:,:,f);
    fx2=sqrt(squeeze(sum(sum(x.^2,1),3)));
    fz2=sqrt(squeeze(sum(sum(x.^2,2),3)));
    at2=auto_corrn(x,3);
    ax2=auto_corrn(x,2);
    az2=auto_corrn(x,1);
    nfp=squeeze(sqrt(mean(mean(a.noisef.^2,1),2)));
    nf2p=squeeze(sqrt(mean(mean(a.noisef2.^2,1),2)));
    msk=a.wnoise>max(a.wnoise(:))/2;
    for j=1:length(nmc)
      x=a.(nmc{j}).^2.*msk;
      nft(:,j)=squeeze(sum(sum(x,1),2))./squeeze(sum(sum(msk,1),2));
      nfx(:,j)=squeeze(sum(sum(x,1),3))./squeeze(sum(sum(msk,1),3));
      nfz(:,j)=squeeze(sum(sum(x,2),3))./squeeze(sum(sum(msk,2),3));
    end
    x=mean(a.wnoise,3).^2;
    nwx=squeeze(sqrt(mean(x,1)));
    nwz=squeeze(sqrt(mean(x,2)));
  case 3
    x=a.noisef(:,:,:,f);
    fx=sqrt(squeeze(sum(sum(sum(x.^2,1),2),4)));
    fy=sqrt(squeeze(sum(sum(sum(x.^2,1),3),4)));
    fz=sqrt(squeeze(sum(sum(sum(x.^2,2),3),4)));
    at=auto_corrn(x,4);
    ax=auto_corrn(x,3);
    ay=auto_corrn(x,2);
    az=auto_corrn(x,1);
    x=a.noisef2(:,:,:,f);
    fx2=sqrt(squeeze(sum(sum(sum(x.^2,1),2),4)));
    fy2=sqrt(squeeze(sum(sum(sum(x.^2,1),3),4)));
    fz2=sqrt(squeeze(sum(sum(sum(x.^2,2),3),4)));
    at2=auto_corrn(x,4);
    ax2=auto_corrn(x,3);
    ay2=auto_corrn(x,2);
    az2=auto_corrn(x,1);
    nfp=squeeze(sqrt(mean(mean(mean(a.noisef.^2,1),2),3)));
    nf2p=squeeze(sqrt(mean(mean(mean(a.noisef2.^2,1),2),3)));
    msk=a.wnoise>max(a.wnoise(:))/2;
    for j=1:length(nmc)
      x=a.(nmc{j}).^2.*msk;
      nft(:,j)=squeeze(sum(sum(sum(x,1),2),3))./squeeze(sum(sum(sum(msk,1),2),3));
      nfx(:,j)=squeeze(sum(sum(sum(x,1),2),4))./squeeze(sum(sum(sum(msk,1),2),4));
      nfy(:,j)=squeeze(sum(sum(sum(x,1),3),4))./squeeze(sum(sum(sum(msk,1),3),4));
      nfz(:,j)=squeeze(sum(sum(sum(x,2),3),4))./squeeze(sum(sum(sum(msk,2),3),4));
    end
    x=mean(a.wnoise,4).^2;
    nwx=squeeze(sqrt(mean(mean(x,1),2)));
    nwy=squeeze(sqrt(mean(mean(x,1),3)));
    nwz=squeeze(sqrt(mean(mean(x,2),3)));
end

nft=sqrt(sum(nft,2));
nfx=sqrt(sum(nfx,2));
nfz=sqrt(sum(nfz,2));

t=a.t(f);

t=t-t(1);
nt=length(t);

disp(sprintf('Looking at noise statistics for %s',n));
disp(sprintf('  %i times dt=%6.4f  %8.6f',nt,mean(diff(t)),std(diff(t))));
meann=sqrt(mean(a.noisef(:).^2));
disp(sprintf('mean noise %8.5f',meann));

sz=size(a.noisef);

figure(1)
subplot(3,1,1);
plot(a.t,nfp/mean(nfp));
ylabel('noisef');
axis('tight');
subplot(3,1,2);
plot(a.t,nf2p/mean(nf2p));
ylabel('noisef2');
axis('tight');
subplot(3,1,3);
plot(a.t,nft,a.t([1 end]),p.noise([1 1]));
ylabel('noise vel');
axis('tight');
xlabel('time');

x=a.x-a.x(1);
z=a.z-a.z(1);

figure(2);clf;
subplot(dims+1,1,1);
plot(t,at/at(1),t,at2/at2(1),t,exp(-t/p.noiseT));
xlabel('t')
axis([0 min(t(end),p.noiseT*5) 0 1]);
subplot(dims+1,1,2);
plot(x,ax/ax(1),x,ax2/ax2(1),x,exp(-(x/p.noiseL).^2));
xlabel('x')
axis([0 min(x(end),p.noiseL*5) 0 1]);
if dims>2
  y=a.y-a.y(1);
  subplot(dims+1,1,3);
  plot(y,ay/ay(1),y,ay2/ay2(1),y,exp(-(y/p.noiseL).^2));
  xlabel('y')
  axis([0 min(y(end),p.noiseL*5) 0 1]);
end
subplot(dims+1,1,dims+1);
plot(z,az/az(1),z,az2/az2(1),z,exp(-(z/p.noiseL).^2));
xlabel('z')
axis([0 min(z(end),p.noiseL*5) 0 1]);

figure(3);clf;
subplot(2,1,1);
plot(a.x,nfx,a.x([1 end]),p.noise([1 1]),a.x,nwx/max(nwx)*p.noise);
axis([a.x(1) a.x(end) 0 inf]);
xlabel('x');

subplot(2,1,2);
plot(a.z,nfz,a.z([1 end]),p.noise([1 1]),a.z,nwz/max(nwz)*p.noise);
axis([a.z(1) a.z(end) 0 inf]);
xlabel('z');

return;




b=zeros(sz);
for j=1:nt
  b(:,:,:,j) = real(ifft(abs(fftn(a.noisef(:,:,:,j))).^2));
end
b=mean(b,4);

figure(3)
plot(a.x,squeeze(mean(mean(b,1),2)),a.y,squeeze(mean(mean(b,1),3)));
keyboard


return

b=zeros(sz);
for j=1:nt
  b(:,:,:,j) = abs(fftn(a.noisef(:,:,:,j))).^2;
end
b=mean(b,4);

dx=double(p.L)/sz(3);
dy=double(p.W)/sz(2);
dz=double(p.H)/sz(1);

w1=dz*(mod((0:sz(1)-1)+floor(sz(1)/2),sz(1))-floor(sz(1)/2));
w2=dy*(mod((0:sz(2)-1)+floor(sz(2)/2),sz(2))-floor(sz(2)/2));
w3=dx*(mod((0:sz(3)-1)+floor(sz(3)/2),sz(3))-floor(sz(3)/2));

w1=abs(fft_modes(sz(1),dz));
w2=abs(fft_modes(sz(2),dy));
w3=abs(fft_modes(sz(3),dx));

maxx=sqrt(max(w1.^2)+max(w2.^2)+max(w3.^2));
N=200;
dx=maxx/(N-1);
w=zeros(1,N);
f=zeros(1,N);
x=zeros(1,N);
for i1=1:sz(1)
  for i2=1:sz(2)
    [ww ff] = jgrid1(sqrt(w1(i1)^2+w2(i2)^2+w3.^2)/dx,squeeze(b(i1,i2,:))',0,N-2,'linear');
    w=w+ww;
    f=f+ff;
  end
end
rg=find(w>0);
w=w(rg);
f=f(rg);
x=(rg-1)*dx;
acx=f./w;%acx=acx/acx(1);
semilogy(x,acx/acx(1),x,exp(-(x/p.noiseL)));
  

bs=mean(abs(b).^2,3);





return;

n=1e6;
x=randn(n,1);


% $$$ sz=size(a.noisef);
% $$$ if strcmp(p.Tz,'Cheb')
% $$$   c=ichebf2c(a.noisef,1);
% $$$   a.z=linspace(0,p.H,p.Nz);
% $$$   N=prod(sz(2:end));
% $$$   c=reshape(c,[sz(1) N]);
% $$$   a.noisef=zeros(p.Nz,N);
% $$$   for i=1:N
% $$$     a.noisef(:,i)=cheval('regular',c(:,i),a.z);
% $$$   end
% $$$   a.noisef=reshape(a.noisef,sz);
% $$$ end
