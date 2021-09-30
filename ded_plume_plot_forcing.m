function a=ded_plume_plot_forcing(nm)
fn=['~/data/dedalus/pm/' nm '/force/force_s1.hdf5'];
a=ded_read_hdf(fn);
if ~isfile(fn)
  disp(sprintf('ded_pm_plot_forcing: %s does not exist',fn));
clf;
ah=jsubplot([5 1],[0.01 0.01],[0.02 0.04],[0.01 0.04]);

if isfield(a,'wz')
  a.ww=a.wz;
  a=rmfield(a,'wz');
end
n={'div','ws','fS','fu','wu','fC','ww','fw','fT'};
n={'fb','fd','fu','wu','wb'};
n=intersect(n,fieldnames(a));
for j=1:length(n)
  x=a.(n{j});
  x=squeeze(mean(x,1));
  axes(ah(j));imagesc(a.z,a.x,x);title(n{j});
end
set(ah,'ydir','normal','dataaspect',[1 1 1],'ytick',[],'xtick',[]);
