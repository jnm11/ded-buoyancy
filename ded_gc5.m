function ded_gc5(n)
fn=sprintf('~/data/dedalus/gc5/gc5_s%i/gc5_s%i_p0.h5',n,n);
if ~isfile(fn)
  disp(fn);
  return;
end
disp(fn);
h5disp(fn);
a.x=h5read(fn,'/scales/x/2');
a.y=h5read(fn,'/scales/y/2');
a.dt=h5read(fn,'/scales/timestep');
a.t=h5read(fn,'/scales/sim_time');
a.i=h5read(fn,'/scales/iteration');
a.b=h5read(fn,'/tasks/b');
a.f=h5read(fn,'/tasks/f');

figure(1);
clf;
y=linspace(a.y(1),a.y(end),length(a.y));
for j=1:length(a.t);
  cla
  b=interp1(a.y,a.b(:,:,j),y,'linear');
  imagesc(a.x,y,b);title(sprintf('t=%5.4f, i=%4u',a.t(j),a.i(j)));
  set(gca,'ydir','normal');
% $$$   hold('on');
% $$$   contour(a.x,a.y,a.f(:,:,j),20);
  drawnow;
  pause(0.01);
end
