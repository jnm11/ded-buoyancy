function a=ded_stats(nm,display)
a=[];
if nargin<2
  display=[];
end
if isempty(display)
  display=nargout==0;
end
a=[];
nmo=nm;
if ~isdir(nm)
  nm=[ded_dedalus_data_dir '/' nm];
end

fnd=nm;

if ~isdir(fnd)
  disp(sprintf('ded_stats: Directory "%s" does not exist',nm));
  return;
end

fns=sprintf('%s/stats.mat',fnd);
fnp=sprintf('%s/param.h5',fnd);
fnm=cellstr_ls(sprintf('%s/mom/mom*.hdf5',fnd));
if isempty(fnm)
  fnm=cellstr_ls(sprintf('%s/xyz/xyz*.hdf5',fnd));
end


if isfile(fnm)
  newer=file_nt(fns,fnm,-5*60);
else
  newer=isfile(fns);
end

if all(newer==1) 
  disp(sprintf('ded_stats: Loading "%s"',fns));
  load(fns);
  if ~display
    return;
  end
else
  disp(sprintf('ded_stats: Making "%s"',fns));
  
  a=ded_read_param(fnp);
  
  if isempty(a)
    disp(sprintf('No parameter file for "%s"',nm));
    return
  end
  
  nz=round(a.H/a.L*a.Nx);
  z=filter_midpoint(linspace(0,a.H,nz+1));
  
  a.mom=ded_read_mom(fnd);
  %  if a.name(1)=='0' | all(a.name(1:2)=='50') | all(a.name(1:2)=='51')
  %  a.T=0;
  %else

  a.T=ded_convergence_T(a.mom);
  %end
  if ~isempty(a.mom)
    if any(diff(a.mom.t)<0)
      disp(sprintf('Non-monotonic time in momentum %s',nm));
    end
    b=ded_gc_find_front_mom(a.mom,a.T);
    a.mX  = b.X;
    a.mt  = a.mom.t;
    a.mXm = b.Xm;
  else
    a.mX  = [];
    a.mt  = [];
    a.mXm = NaN;
  end

  a.flux=ded_read_flux(nm);
  if ~isempty(a.flux)
    if any(diff(a.flux.t)<0)
      disp(sprintf('Non-monotonic time in fluxes %s',nm));
    end
    b=ded_gc_find_front(a.flux,a.T);
    a.fX  = b.X;
    a.ft  = a.flux.t;
    a.fXm = b.Xm;
  else
    a.fX  = [];
    a.ft  = [];
    a.fXm = NaN;
  end
  a.stats=ded_read_stats(nm);
  if ~isempty(a.stats)
    f=find(a.stats.t>=a.T);
    if isempty(f)
      f=length(a.stats.t);
    end
    a.g=mean(a.stats.g(f));
    a.sg=std(a.stats.g(f));
    a.X=mean(a.stats.X(f));
    a.sX=std(a.stats.X(f));
    a.U=mean(a.stats.U(f));
    a.sU=std(a.stats.U(f));
  else
    a.g=NaN;
    a.sg=NaN;
    a.X=NaN;
    a.sX=NaN;
    a.U=NaN;
    a.sU=NaN;
  end
  if any(nm=='/')
    [dd nm]=fileparts(nm);
  end
  a.msg=sprintf('nm=%s, t=%6.1f, T=%6.1f, X=%6.3f',nm,a.t, a.T, a.X);
end

if display
  clf;
  subplot(2,1,1);
  h=plot(a.ft([1 end]), a.fXm([1 1]), a.ft,a.fX);
  lg={'mean thresh','Threshold'}
  if ~isempty(a.stats) 
    h(end+1)=line(a.stats.t,a.stats.X,'color',[0 0 0]);
    lg{end+1}='stats';
  end
  axis('tight');aa=axis;aa(3)=aa(3)-0.1;aa(4)=aa(4)+0.1;axis(aa);         
  line([a.T a.T],aa(3:4));
  legend(h,lg,'location','best')
  xlabel('t');
  title(sprintf('Re=%5.0f,H=%3.1f,h=%3.1f,U=%5.3f,PIDX=%4.1f,X=%5.2f',a.Re,a.H,a.h,a.U,a.PIDX,a.fXm));
  subplot(2,1,2);
  h=plot(a.mt([1 end]),  a.mXm([1 1]),  a.mt, a.mX);
  lg={'mean Moment','Moment'};
  if ~isempty(a.stats) 
    h(end+1)=line(a.stats.t,a.stats.X,'color',[0 0 0]);
    lg{end+1}='stats';
  end
  axis('tight');aa=axis;aa(3)=aa(3)-0.1;aa(4)=aa(4)+0.1;axis(aa);         
  line([a.T a.T],aa(3:4));
  legend(h,lg,'location','best')
  xlabel('t');
end
a.forcing = double(a.forcing);
a.Nx = double(a.Nx);
a.Ny = double(a.Ny);
a.Nz = double(a.Nz);

a.U  = double(a.U);
a.U1 = double(a.U1);
a.L  = double(a.L);
a.W  = double(a.W);
a.H  = double(a.H);
a.Re = double(a.Re);
if isfield(a,'PIDX')
  a.PIDX = double(a.PIDX);
else
  a.PIDX=0;
end
a=orderfields(a);
save(fns,'a');


return;
a=ded_stats('121');
