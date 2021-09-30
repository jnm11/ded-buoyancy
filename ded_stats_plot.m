function ded_stats_plot(nm)
if iscell(nm) 
  for j=1:length(nm)
    figure(j);
    ded_stats_plot(nm{j});
  end
  return;
end
if isstruct(nm) & length(nm)>1
  for j=1:length(nm)
    figure(j);
    ded_stats_plot(nm(j));
  end
  return;
end
if isstruct(nm)
  a=nm;
else
  a=ded_stats(nm);
end
% $$$ if ~isempty(a.flux)
% $$$   h=plot(a.flux.t([1 end]), a.fXm([1 1]), a.flux.t([1 end]), a.fX([1 1]), 'r-');%,...
if ~isempty(a.stats) 
  trg=[min(a.stats.t) max(a.stats.t)];
  if a.PIDX>0 
    plot(a.stats.t,a.stats.X,trg,a.X([1 1]),trg,a.PIDX([1 1]));
  else
    plot(a.stats.t,a.stats.X);
  end
  axis('tight');aa=axis;
  aa(3:4)=(aa(3)+aa(4) + 1.05*(aa(4)-aa(3))*[-1 1])/2;
  axis(aa);         
  line([a.T a.T],aa(3:4));
  %  legend(h,{'mean thresh','Threshold','mean Moment','Moment'},'location','best')
  xlabel('t');
  title(sprintf('H=%5.1f, U=%5.3f, g=%5.3f, %s %s',a.H,a.U,a.g,a.series,a.msg),'fontsize',8);
end
