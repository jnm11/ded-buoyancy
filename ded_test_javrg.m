function ded_test_javrg(fns)

if nargin<1
  fns=[];
end

% mpiexec -n 10 ded_gc.py --preset ccle --Nx 4096 --Ny 1 --Nz None -W None -T 20 --Re 500 --lbc slip --ubc slip --reset gc/test/01 


figure(1);
if isempty(nms)
  nms=cellstr_ls([ded_dedalus_data_dir '/gc/test/*']);
  nms=cell_rm_prefix([ded_dedalus_data_dir '/'],nms);
end
if ~iscell(nms)
  nms={nms};
end

for k=1:length(nms)
  nm=nms{k};
  fn=cellstr_ls(sprintf('%s/%s/javrg/javrg-*.hdf5',ded_dedalus_data_dir,nm));
  p=ded_read_param(nm);
  for j=1:length(fns)
    fn=fns{j};
    a=ded_read_hdf(fn);
    b=ded_zgrid(a,[],{'b','S','Wy'},{},{},{},p.H);
    clf;
 
    subplot(3,1,1);
    imagesc(b.x,b.z,b.b/a.dt);
    colorbar
    title(sprintf('%s, t=[%6.2f %6.2f], n=%u',nm,a.t1,a.t2,a.write_num));
    set(gca,'ydir','normal','dataaspectratio',[1 1 1],'position',[0 0.68 1 0.3]);

    subplot(3,1,2);
    imagesc(b.x,b.z,b.S/a.dt);
    colorbar
    set(gca,'ydir','normal','dataaspectratio',[1 1 1],'position',[0 0.35 1 0.3]);

    subplot(3,1,3);
    imagesc(b.x,b.z,b.Wy/a.dt);
    colorbar
    set(gca,'ydir','normal','dataaspectratio',[1 1 1],'position',[0 0.02 1 0.3]);
    drawnow;
  end
end




%mpiexec -n 4 ded_gc.py --presetccle --Nx 1536 --Ny 1 -W None -T 20 --Re 500 --lbc slip --ubc slip --reset gc/test/01 
