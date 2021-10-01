function a=ded_gc_read_param(num)

fns{1}=num;
fns{2}=[num '/param.h5'];;
fns{3}=['~/data/dedalus/gc/' num '/param.h5'];
fns{4}=[num '/param.h5'];
fn=[];
for j=1:length(fns)
  if isfile(fns{j})
    fn=fns{j};
    dd=fileparts(fn);
    [dd num]=fileparts(dd);
    break;
  end
end
if ~isfile(fn)
  disp(sprintf('No match found for "%s"',num));
   a=[];
  return;
end

if isempty(dd)
  dd='./';
end

fnt=sprintf('%s/%s/time',dd,num);
fns=sprintf('%s/%s/status',dd,num);
fnj=sprintf('%s/%s/jobid',dd,num);

%mm={'H','L','Nx','Ny','Nz','Re','Sc','T','U','U1','W','average','h','lbc','name','ubc','parallel'};
if isfile(fn)
  i=h5info(fn);
  mm={i.Datasets.Name};
  for j=1:length(mm)
    m=mm{j};
    x=h5read(fn,['/' m]);
    if iscell(x) & length(x)==1
      a.(m)=x{1};
    elseif isempty(x)
      %
    elseif isnumeric(x)
      a.(m)=double(x);
    end
  end
end

m={'parallel','signals','ck_time','forcing','Smooth','average','Wx','Wz','x0','x1','x2','x3'};
for j=1:length(m)
  if ~isfield(a,m{j})
    a.(m{j})=NaN;
  end
end
a.forcing(~isfinite(a.forcing))=1;

if isfile(fnt)
  for j=1:10
    fp=fopen(fnt);
    a.t=fscanf(fp,'%f');
    fclose(fp);
    if ~isempty(a.t)
      break;
    end
    pause(0.1);
  end
  if isempty(a.t)
    a.t=NaN;
  end
else
  a.t=NaN;
end

if isfile(fnj)
  fp=fopen(fnj);
  a.jobid=fscanf(fp,'%f');
  fclose(fp);
else
  a.jobid=[];
end
if isfile(fns)
  fp=fopen(fns);
  a.status=fscanf(fp,'%s');
  fclose(fp);
else
  a.status=[];
end

a=orderfields(a);

return
