function ded_gca_2(nm)

c=ded_read_avrg(nm);
p=ded_read_param(nm);
s=ded_read_stats(nm);

%[ts sc] = ded_interp_stats(s,p,c.t);
fg=find(s.t>=c.t(1) & s.t<=c.t(2));
if isfield(s,'g');
  g=mean(s.g(fg));
  disp(std(s.g(fg)));
else
  g=p.g;
end


c=ded_zgrid(c,c.nz,[],{'u','w','au','aw','auu','aww','avv','auw'},{'u','w','au','aw'},{'ab','ap'},p.H);

dp=c.ap(1,:)-c.ap(end,:);

B=c.abI(end,:);
P=c.apI(end,:);

Yp=1.05*max(dp(:));
YB=1.05*max(B(:));
xrg=[0 p.L];

x=c.x;
z=c.z;
dx=x(2)-x(1);
c.audx =pr_diff(c.au,dx,2);
c.awdx =pr_diff(c.aw,dx,2);
c.audxx=pr_diff(c.audx,dx,2);
c.awdxx=pr_diff(c.awdx,dx,2);
c.audxz=pr_diff(c.audz,dx,2);
c.awdxz=pr_diff(c.awdz,dx,2);

%sum(sum(c.audx.*c.awdz))./sqrt(sum(sum(c.audx.^2))*sum(sum(c.awdz.^2)))
div=c.audx+c.awdz;

omega=c.awdx-c.audz;
c.av=0;
uu=max(0,c.auu-c.au.^2);
if isfield(c,'avv')
  vv=max(0,c.avv-c.av.^2);
else
  vv=0;
end
ww=max(0,c.aww-c.aw.^2);
TKE=(uu+vv+ww)/2;
xrg=[0 p.L];

figure(1);clf;
subplot(3,1,1);
plot(x,dp,xrg,p.hu*[1 1],xrg,p.hb*[1 1]);
ylabel('p|^H_0','interpreter','tex')
title('Pressure drop');
axis([xrg 0 Yp]);
subplot(3,1,2);
plot(x,B,xrg,p.hu*[1 1],xrg,p.hb*[1 1]);
ylabel('\int b dz','interpreter','tex')
axis([xrg 0 YB]);
title('Mass hold up');
subplot(3,1,3);
plot(x,dp-B*g,xrg,[0 0])
axis([xrg -inf inf]);
ylabel('p|^H_0-g\int b dz','interpreter','tex');
title('Hydrostatic balance');

figure(2);clf;
imagesc(x,z,omega);
set(gca,'ydir','normal','clim',percentile(omega(:),[1 99]));
title('y vorticity');
colorbar;

figure(3);clf;
imagesc(x,z,div);
set(gca,'ydir','normal','clim',percentile(div(:),[1 99]));
title('Incompressibility');
colorbar;

figure(4);clf;
imagesc(x,z,TKE);
set(gca,'ydir','normal','clim',percentile(TKE(:),[1 99]));
title('Turbulent Kinetic Energy');
colorbar;

figure(5);clf;
imagesc(x,z,c.ab);
set(gca,'ydir','normal','clim',[0 1]);
title('Density');
colorbar;
[cx,cy]=longest_contours(x,z,c.ab,0.05:0.1:0.95,10); for k=1:length(cx);line(cx{k},cy{k},'color',[1 1 1]*.8);end


p=c.ap-repmat(c.ap(1,:),length(c.z),1)+c.abI; % Take of the bottom pressure
%p=c.ap+c.abI;;
prg=percentile(p(:),[1 99]);
figure(6);clf;
imagesc(x,z,p);
set(gca,'ydir','normal','clim',prg);
title('deviation from hydrostatic');
colorbar;
[cx,cy]=longest_contours(x,z,p,linspace(prg(1),prg(2),10),10); for k=1:length(cx);line(cx{k},cy{k},'color',[1 1 1]*.8);end


p=c.ap+c.abI+c.auu/2+c.aww/2;
p=p-mean(p(:,end));

prg=percentile(p(:),[1 99]);
figure(7);clf;
imagesc(x,z,p);
set(gca,'ydir','normal','clim',prg);
title('Total pressure');
colorbar;
[cx,cy]=longest_contours(x,z,p,linspace(prg(1),prg(2),10),10); for k=1:length(cx);line(cx{k},cy{k},'color',[1 1 1]*.8);end


return;
