function ded_test_wdu(nm)
nm='pm/test/tokachi/24';
p=ded_read_param(nm);
c=ded_coord(nm);
fyz=ded_get_fn(nm,'force/fyz');
fyz=ded_read_hdf(fyz{1});
fx=ded_get_fn(nm,'force/fyz');
fx=ded_read_hdf(fx{1});

fns=ded_get_fn(nm,'wdu');
%for j=1:length(fns)

for j=1:length(fns)
  a=ded_read_hdf(fns{j});
  if isfield(a,'dd');a.fd=a.dd;end;
  wyz=fyz.wdrn;
  wyz(wyz<0.999)=0;
  w=a.wu.*wyz;
  maxwu=squeeze(max(max(a.wu,[],1),[],2));
  fd=pr_diff(a.fu,c.dAx,3,1,[],-1);
  sw=squeeze(sum(sum(w,1),2));
  C1(:,j) = squeeze(sum(sum(w.*fd.*a.fd,1),2)./sqrt(sum(sum(w.*fd.^2,1),2).*sum(sum(w.*a.fd.^2,1),2)));
  C2(:,j) = sqrt(squeeze(sum(sum(w.*(fd-a.fd).^2,1),2))./sw);
  C2(maxwu<2e-3,j)=0;
  C1(maxwu<2e-3,j)=1;
end

clf;
subplot(2,1,1);plot(c.Ax,C1);
subplot(2,1,2);plot(c.Ax,C2);

return;
plot(c.Ax,squeeze(mean(mean(a.wu,1),2)),c.Ax,squeeze(max(max(abs(a.fu),[],1),[],2)));

figure;
clf;
for j=1:c.NAx
  subplot(2,2,1);mesh(c.Az,c.Ay,fd(:,:,j));            axis('tight');title('fd');
  subplot(2,2,2);mesh(c.Az,c.Ay,a.dd(:,:,j));          axis('tight');title('a.fd');
  subplot(2,2,3);mesh(c.Az,c.Ay,w(:,:,j).*(a.dd(:,:,j)-fd(:,:,j)));axis('tight');title('diff');
  drawnow;
end
