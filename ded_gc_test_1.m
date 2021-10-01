


fns={'01','02','03','04','10','11','12','13','14','15','16','19','20','21','22','23'};
for j=1:length(fns);
  try
    p=ded_read_param(nm);
    disp(sprintf('%s, Re=%5.0f, Nx=%4.0f , Ny=%4.0f, Nz=%4.0f, W= %8.6f',p.name,p.Re,p.Nx,p.Ny,p.Nz,p.W));  
    ded_gca_1(nm);
    sfigure(4);clf
    ded_gc_check_bc(nm);
    drawnow;
    nm=['gc/f6/' fns{j}];
    pause;
  end
end
