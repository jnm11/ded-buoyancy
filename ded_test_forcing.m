%mpiexec -n 4 ded_gc.py --dtb 0.1 --dtslice 0.1 --slicex 0 --preset gcf6 --PIDG 1 --Nx 256 --Ny 1 -W None -T 10 --Re 50 --lbc slip --ubc slip gc/test/100

fns=cellstr_ls([ded_dedalus_data_dir '/gc/slice/50']);
clf;



clear('x');
for j=1:length(fns)
  nm=fns{j};
  
  s=ded_read_g(nm,'slicex');
  p=ded_read_param(nm);
  if j==1
    nd=2+isfield(s,'v');
    for k=1:nd+1
      a(k)=subplot(nd+1,1,k);hold('on');
    end
  end  
  
  t=s.t;
  x(j,:)=[p.x0 p.x1 p.x2 p.x3 p.x4 p.x5 p.x6 p.x7];
  nz=size(s.z);
  cu=ichebf2c(s.u,1);cu(nz+1:2*nz,:)=0;
  if isfield(s,'v')
    cv=ichebf2c(s.v,1);cv(nz+1:2*nz,:)=0;
  else
    cv=0*cu;
  end
  cw=ichebf2c(s.w,1);cw(nz+1:2*nz,:)=0;
  cb=ichebf2c(s.b,1);cb(nz+1:2*nz,:)=0;
  
  u=ichebc2f(cu,1)+p.U;
  v=ichebc2f(cv,1);
  w=ichebc2f(cw,1);
  b=ichebc2f(cb,1);
  
  eu=sqrt(ichebDint(u.^2,1));
  ev=sqrt(ichebDint(v.^2,1));
  ew=sqrt(ichebDint(w.^2,1));
  eb=sqrt(ichebDint(b.^2,1));
  axes(a(1));hx(j)=plot(t,eu);
  if nd==3 
    axes(a(nd-1));
    hy(j)=plot(t,ev);
  end
  axes(a(nd));hz(j)=plot(t,ew);
  axes(a(nd+1));hb(j)=plot(t,eb);
end
set(a,'yscale','log');
for k=1:nd+1
  axes(a(k));
  hold('on');
  axis('tight');
end





return;

nm=[ded_dedalus_data_dir '/gc/slice/13'];
s=ded_read_g(nm,'slicex');
p=ded_read_param(nm);

