fns=cellstr_ls('~/data/dedalus/gc-ccle-046-*.mat');
clf;
for j=1:length(fns)
  load(fns{j});
  %imagesc(a.x,a.z,a.b-a.bb);set(gca,'ydir','normal','dataaspect',[1 1 1]);set(gca,'position',[0 0 1 1]);
  %contour(a.x,a.z,a.b,[0.03 0.03]);
  contour(a.x,a.z,a.b,[0.03 0.03]);
  hold('on');
end
set(gca,'ydir','normal','dataaspect',[1 1 1]);set(gca,'position',[0 0 1 1]);
