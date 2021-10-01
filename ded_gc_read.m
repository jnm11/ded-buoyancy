function a=ded_gc_read(nm,cmd)
if nargin<2
  cmd={};
end
a.param=ded_read_param(nm);
if isempty(a.param)
  disp(sprintf('ded_gc_read: Simulation does not exist %s',nm));
  a=[];
  return;
end
a.gc = ded_read_2d(nm,'gc');
%a.y  = ded_read_2d(nm,'y');
a.y  = ded_read_avrg(nm);  
nz=round(a.param.H/a.param.L*a.param.Nx);
a.y=ded_zgrid(a.y,nz,[],{'au','aw','auu','aww','avv','auw'},{'au','aw'},{'ab','ap'},a.param.H);
a.L.x=a.y.x;
a.L.z=a.y.z;
a.R.x=a.y.x;
a.R.z=a.y.z;
fnm=setdiff(fieldnames(a.y),{'x','z','t','b','bI','p','pI','pdzz','u','uI','udz','v','vI','vdz','w','wI','wdz'});
dt=diff(a.y.t);
for j=1:length(fnm)
  a.L.(fnm{j}(2:end))=a.y.(fnm{j})(1,:)/dt;
  a.R.(fnm{j}(2:end))=a.y.(fnm{j})(end,:)/dt;
end

%a.L  = ded_read_2d(nm,'left');
%a.R  = ded_read_2d(nm,'right');
a.yz = ded_read_2d(nm,'yz');
a.s  = ded_read_stats(nm);

if isempty(a.L) & ~isempty(a.gc)
  fg=fieldnames(a.gc); 
  for k=1:length(fg)
    b=fg{k};
    if b(1)=='L'
      a.L.(b(2:end))=a.gc.(b);
      a.gc=rmfield(a.gc,b);
    elseif b(1)=='R'
      a.R.(b(2:end))=a.gc.(b);
      a.gc=rmfield(a.gc,b);
    end
  end
  a.L.t=a.gc.t;
  a.L.x=a.gc.x;
  a.R.t=a.gc.t;
  a.R.x=a.gc.x;
end
if isempty(a.yz) & ~isempty(a.gc)
  fg=fieldnames(a.gc); 
  for k=1:length(fg)
    b=fg{k};
    if length(b)>3
      if all(b(1:3)=='Iyz')
        a.yz.(b(4:end))=a.gc.(b);
        a.gc=rmfield(a.gc,b);
      end
    end
  end
  a.yz.t=a.gc.t;
  a.yz.x=a.gc.x;
end
f={'y','L','R','yz'};
for j=1:length(f)
  if isempty(a.(f{j}))
    a=rmfield(a,f{j});
  else
    fg=setdiff(fieldnames(a.(f{j})),{'x','y','z','sz','nx','ny','nz','nt'});
    for k=1:length(fg)
      x=a.(f{j}).(fg{k});
      tz=size(x);
      if ~isfield(a.(f{j}),'nt')
        continue;
      end
      nd=max(find(tz==a.(f{j}).nt));
      if isempty(nd)
        continue;
      end
      switch(cmd)
        case 'meant'
          a.(f{j}).(fg{k})=mean(a.(f{j}).(fg{k}),nd);
        case 'finalt'
          switch(nd)
            case 1
              a.(f{j}).(fg{k})=a.(f{j}).(fg{k})(end);
            case 2
              a.(f{j}).(fg{k})=a.(f{j}).(fg{k})(:,end);
            case 3
              a.(f{j}).(fg{k})=a.(f{j}).(fg{k})(:,:,end);
            case 4
              a.(f{j}).(fg{k})=a.(f{j}).(fg{k})(:,:,:,end);
          end
      end
    end
  end
end

if ~isempty(a.gc)
  a.x=a.gc.x;
  a.t=a.gc.t;
elseif isfield(a,'y')
  a.x=a.y.x;
  a.t=a.y.t;
end
if isfield(a,'y')
  a.z=a.y.z;
  if isfield(a.y,{'Ex','Ey','Ez'});
    a.y.E=a.y.Ex+a.y.Ey+a.y.Ez;
    a.y.OW=a.y.S-a.y.E;
  end
end
if ~isfield(a.y,'E')
  dx=a.x(2)-a.x(1);
  dimx=ndims(a.y.au);
  a.y.audx=pr_diff(a.y.au,dx,dimx);
  a.y.awdx=pr_diff(a.y.aw,dx,dimx);
  a.y.E=(2*a.y.audx.^2+2*a.y.awdz.^2+(a.y.audz+a.y.awdx).^2)/a.param.Re;
end



if isfield(a.s,'g')
  a.g=mean(a.s.g(a.s.t>=a.t(1)))
else
  a.g=a.param.g;
end

%function c=NSx(c) % Calculate Navier Stokes terms in x direction
