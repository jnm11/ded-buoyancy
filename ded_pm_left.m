function a=ded_pm_left(n)
fns=ded_get_fn(n,'left');
if isempty(fns)
  a=[];
  return;
end


for k=1:length(fns)
  a=ded_read_hdf(fns{k})
  for j=1:size(a.u,4)
    subplot(2,3,1);mesh(a.x,a.y,squeeze(a.u(:,:,:,j)));axis('tight');title(sprintf('t=%6.3f u',a.sim_time(j)));
    subplot(2,3,2);mesh(a.x,a.y,squeeze(a.v(:,:,:,j)));axis('tight');title('v')
    subplot(2,3,3);mesh(a.x,a.y,squeeze(a.w(:,:,:,j)));axis('tight');title('w')
    subplot(2,3,4);mesh(a.x,a.y,squeeze(a.b(:,:,:,j)));axis('tight');title('b')
    subplot(2,3,5);mesh(a.x,a.y,squeeze(a.noisef(:,:,:,j)));axis('tight');title('noisef')
    drawnow;
    pause(0.1);
  end
end


  