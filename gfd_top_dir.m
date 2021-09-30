function  d=gfd_top_dir() 
fns={'~/Dropbox/GFD-2019-Gravity-Currents'
   'C:\Users\Claudia\Dropbox\GFD 2019 Gravity Currents'};
for j=1:length(fns)
  if exist(fns{j},'dir')
    d=fns{j};
  end
end
