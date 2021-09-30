function ded_plot_noise(nm)
fns=ded_get_fn(nm,'noise');
p=ded_read_param(nm);
for j=1:length(fns)
  fn=fns{j};
  a=ded_read_hdf(fn);
  [b f]=max(max(max(a.noisey.^2+a.noisez.^2)));
  b=sqrt(b);
  subplot(2,2,1);
  mesh(a.y,a.z,a.noisey(:,:,f));
  axis([-p.W/2 p.W/2 -p.H/2 p.H/2 -b b]);
  subplot(2,2,2);
  mesh(a.y,a.z,a.noisez(:,:,f));
  axis([-p.W/2 p.W/2 -p.H/2 p.H/2 -b b]);
  subplot(2,2,3);
  mesh(a.y,a.z,a.noise(:,:,f));
  b=max(abs(a.noise(:)));
  axis([-p.W/2 p.W/2 -p.H/2 p.H/2 -b b]);
  drawnow
end

