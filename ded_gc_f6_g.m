xffunction ded_gc_f6_g
% find the g that will give X=20
clear all
fns={'gc/f6/g/16','gc/f6/g/23','gc/f6/g/22','gc/f6/g/26','gc/f6/g/27','gc/f6/g/15'};
dc;ded_plot_PID(fns,1);
figure;
for j=1:length(fns)
  a=ded_read_stats(fns{j});
  p=ded_read_param(fns{j});
  if isempty(p)
    return;
  end
  X1(j)=a.X(end);
  X2(j)=mean(a.X(round(end/2):end));
  g(j)=p.g;
end
plot(g,X1,'b-s',g,X2,'r-s');

%ded_gca_2('gc/f6/g/18');
