function ded_gc_tint(nm)
if nargin==0
  nm={'120','122','124','126','128','130','132','134','136','138','140','142','144','146','148',...
      '121','123','125','127','129','131','133','135','137','139','141','143','145','147','149'};
  for j=1:length(nm)
    ded_gc_tint(nm{j});
  end
  return;
end
nm=['gc/' nm];
p=ded_read_param(nm);
a=ded_read_g(nm,'tint');
T=ded_convergence_T(nm);
if isempty(T)
  T=a.t(end-1);
end
f=min(find(a.t>T));
fnm={'ab','ap','au','auu','auw','avv','aw','aww'};
n=length(fnm);
dt=a.t(end)-a.t(f);
zz=linspace(0,p.W,4*length(a.z));
for k=1:n
  n=fnm{k};
  X=(a.(n)(:,:,end)-a.(n)(:,:,f))/dt;
  a.(n)=X;
  b.(n)=interp1(a.z,X,zz,'spline','extrap');
end
b.x=a.x';
b.z=zz';


[dd nm]=fileparts(nm);
fn=['~/Dropbox/Jim-Claudia-GC/mat-files/' nm '.mat'];

save(fn,'b','p');


return
% rsync -vaP hamilton:Dropbox/Jim-Claudia-GC/mat-files ~/Dropbox/Jim-Claudia-GC
typ{end+1}='P';
b.P=b.ap+b.auu/2+b.aww/2;
ntyp=ntyp+1;

figure(1);clf;
ah=jsubplot([1 ntyp],[0.01 0.01],[0.01 0.01],[0.01 0.01]);
for k=1:length(typ)
  axes(ah(k));
  n=typ{k};
  imagesc(b.x,b.z,b.(n));
  colorbar
  h=text(xrg(end)*0.95,p.H*0.95,n,'VerticalAlignment','top','HorizontalAlignment','right');
end
set(ah,'ydir','normal','xtick',[],'ytick',[]);
set(ah,'xlim',xrg,'ylim',[0 p.H],'dataaspect',[1 1 1]);



a.P=a.ap+a.auu/2+a.aww/2+a.avv/2;
a.P=a.P-mean(a.P(:,end));
minP=min(a.P(:));
maxP=max(a.P(:));


nz=length(a.z);
nx=length(a.x);
t1=(0.5:nz-0.5)'/nz;
t2=(0:nz)'/nz;
z1=sin(pi*t/2).^2; 
z2=sin(pi*t2/2).^2; 
plot(t1,a.z,t1,z1,t2,z2);
plot(t,acos(1-2*a.z)/pi-t);
b=real(fft([a.ab(end:-1:1,1);a.ab(:,1)]));


dz = diff([0;a.z;p.H]);
bI  = cumsum(repmat(dz(1:end-1),1,nx).*filter_midpoint(a.ab([1 1:end],:)));
plot(a.z,bI(:,1),a.z,a.P(:,1));


dz=[a.z(1);filter_midpoint(diff(a.z));p.H-a.z(end)];


plot(a.z,a.ab(:,1))



figure(2);clf;
subplot(2,1,1);
h=plot(a.x,a.P([1 end],:));
xlabel('x');
axis([xrg minP maxP]);
legend(h([1 end]),{'bottom','top'});
title('Total pressure');
subplot(2,1,2);
h=plot(a.z,a.P(:,[1 end]));
legend(h([1 end]),{'left','right'});
xlabel('z');
axis([zrg minP maxP]);

for j=1:size(a.P,1)
  plot(a.x,a.P(j,:));
  axis([xrg minP maxP]);
  title(sprintf('z=%6.3f',a.z(j)));
  drawnow
  pause;
end

keyboard
a.ap=a.ap-mean(a.ap(:,end));
minb=min(bI(:));
maxb=max(bI(:));
dp=a.ap+bI;
mindp=min(dp(:));
maxdp=max(dp(:));

for j=1:nx
  subplot(2,1,1);
  plot(a.z,a.ap(:,j),a.z,bI(:,j));
  axis([zrg min(minP,minb) max(maxP,maxb)]);
  subplot(2,1,2);
  plot(a.z,dp(:,j));
  axis([zrg mindp maxdp]);
  title(sprintf('x=%6.3f',a.x(j)));
  drawnow
end

plot(a.x,bI(end,:),p.flux.x,p.flux.B);

plot(p.flux.x,p.flux.P+p.flux.B);

