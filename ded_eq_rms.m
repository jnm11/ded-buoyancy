function ded_eq_rms(rg,aa,lbl,Px,dx,H)

nz = size(aa,1);
nx = size(aa,2);
w  = ichebintw(nz)';
z  = ichebgrid(nz);
if nargin>3
  ded_eq_rms(rg,aa,lbl);
  Ix=pr_int(aa,dx,2,[],Px);
  Ix=Ix-Ix(:,rg(end),:);
  disp('Int x');
  ded_eq_rms(rg,Ix,lbl);
  Iz=ichebintf(aa,1);
  Iz=Iz-Iz(1,:,:);
  disp('Int z');
  ded_eq_rms(rg,Iz,lbl);
  return;
end


% Construct orthonormal polynomial projection
z0 = ones(size(z))/sqrt(2); 
z1 = z*sqrt(3/2);
z2 = z.^2-1/3;
z2 = z2./sqrt(w'*z2.^2);
O=([z0 z1 z2].*w)'*[z0 z1 z2];

aa=cat(3,sum(aa,3),aa);

a.rms   = squeeze(sqrt(mean(mean(aa(:,rg,:).^2,1),2)))';
a.rmst  = squeeze(sqrt(mean(   aa(end,rg,:).^2   ,2)))';
a.rmsb  = squeeze(sqrt(mean(   aa(  1,rg,:).^2   ,2)))';
a.rmsw0 = squeeze(sqrt(mean(sum(w.*z0.*aa(:,rg,:),1).^2,2)))';
a.rmsw1 = squeeze(sqrt(mean(sum(w.*z1.*aa(:,rg,:),1).^2,2)))';
a.rmsw2 = squeeze(sqrt(mean(sum(w.*z2.*aa(:,rg,:),1).^2,2)))';

s=cellsprintf('%7s',lbl);
disp(['   ' [s{:}]]);
disp(['rms ' sprintf('%6.4f ',a.rms  )]);
disp(['top ' sprintf('%6.4f ',a.rmst )]);
disp(['bot ' sprintf('%6.4f ',a.rmsb )]);
disp(['z0  ' sprintf('%6.4f ',a.rmsw0)]);
disp(['z1  ' sprintf('%6.4f ',a.rmsw1)]);
disp(['z2  ' sprintf('%6.4f ',a.rmsw2)]);

