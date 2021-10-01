function a=ded_gc_stats(nm,display)

if nargin<2
  display=[];
end
if isempty(display)
  display=nargout==0;
end
a=[];
fnd=sprintf('~/gc/%s',nm);
fns=sprintf('%s/stats.mat',fnd);
fnp=sprintf('%s/param.h5',fnd);
fnm=cellstr_ls(sprintf('%s/mom/mom*.hdf5',fnd));

if ~isdir(fnd)
  disp(sprintf('ded_gc_stats: Directory "%s" does not exist',nm));
  return;
end

if isempty(fnm)
  disp(sprintf('ded_gc_stats: Non momentum files found for "%s"',nm));
  return;
end

newer=file_nt(fns,fnm);

if all(newer==1) 
  disp(sprintf('ded_gc_stats: Loading "%s"',fns));
  load(fns);
  if ~display
    return;
  end
else
  disp(sprintf('ded_gc_stats: Making "%s"',fns));
  
  a=ded_gc_read_param(fnp);
  if isempty(a)
    disp(sprintf('No parameter file for "%s"',nm));
    return
  end
  
  
  U=-a.U;
  S=1024/a.L;
  HT=min(2,a.H);
  
  nz=round(a.H/a.L*a.Nx);
  z=filter_midpoint(linspace(0,a.H,nz+1));
  
  a.mom=ded_read_mom(nm);
  a.flux=ded_read_flux(nm);
  
  [a.T dt bb]=ded_convergence_T(nm);
  a.X=NaN;
  a.MB=NaN;
  a.MX=NaN;
  a.SX=NaN;
  if ~isempty(a.flux)
    a.maxB =max(abs(a.flux.B(:)));
    a.maxUW=max(abs(a.flux.UW(:)));
    a.maxVV=max(abs(a.flux.VV(:)));
    a.maxP =max(abs(a.flux.P(:)));
    a.maxQ =max(abs(a.flux.Q(:)));
    
    for j=1:size(a.flux.B,2)
      fx=find(a.flux.B(:,j)>a.maxB*0.1);
      a.flux.X(j)=NaN;
      if ~isempty(fx) & max(fx)<=length(a.flux.x)
        a.flux.X(j)=max(a.flux.x(fx));
      end
    end
    a.mom.X = a.mom.MX+sqrt(3)*a.mom.SX;
    
    f=findmin(abs(a.mom.t(end)-a.mom.t-100));
    XX=a.mom.X(f:end);XX=XX(:)';
    NN=length(XX);
    t=linspace(0,1,NN);
    a.T0=a.mom.t(f);
    if NN==1
      a.X0=XX;
    else
      PP=polyfit(t,XX,1);
      a.X0=polyval(PP,0);
    end
    
    f=find(a.mom.t>=min(a.T,a.mom.t));
    a.MB=sum(a.mom.B(f));
    a.MX=sum(a.mom.BX(f))/a.MB;
    a.SX=sqrt(sum(a.mom.BXX(f))/a.MB-a.MX^2);
    a.X = a.MX+sqrt(3)*a.SX;
    %    plot(t,XX,t,polyval(PP,t));
  end
  if ~isfinite(a.X)
    if ~isempty(a.mom.X)
      a.X=a.mom.X(end);
    end
  end
  a.msg=sprintf('nm=%s, t=%6.1f, T=%6.1f, X=%4.1f',nm, a.t, a.T, a.X);
end

if display
  clf;
  h=plot(a.flux.t([1 end]),a.X([1 1]),a.flux.t,a.flux.X,a.mom.t,a.mom.X);
  axis('tight');aa=axis;aa(3)=aa(3)-0.1;aa(4)=aa(4)+0.1;axis(aa);         
  line([a.T a.T],aa(3:4));
  legend(h,{'mean','Threshold','Moment'},'location','best')
  xlabel('t');
  title(sprintf('H=%5.1f,U=%5.3f, %s',a.H,-a.U,a.msg));
end

a=orderfields(a);
save(fns,'a');


return;
a=ded_gc_stats('121');
