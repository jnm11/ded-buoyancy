function a=ded_read_param(n)

DD=ded_dedalus_data_dir;
quiet=1;
a=[];
if iscell(n)
  for j=1:length(n)
    aa=ded_read_param(n{j});
    if j==1
      a=aa;
      continue;
    end
    if isempty(aa)
      disp(sprintf('ded_read_param: empty parameter file %s',n{j}));
    else
      [dd nm1]=fileparts(n{j});
      [dd nm2]=fileparts(dd);
      [dd nm3]=fileparts(dd);
      a=struct_array_append(a,aa,[nm3 '/' nm2],quiet);
    end
  end
  return;
end


fn1=n;
fn2=sprintf('%s/param.h5',n);
fn3=sprintf('%s/%s/param.h5',DD,n);
fn4=sprintf('%s/results/%s/param.mat',DD,n);

if isfile(fn1)
  fn=fn1;
elseif isfile(fn2)
  fn=fn2;
elseif isfile(fn3)
  fn=fn3;
elseif isfile(fn4)
  a=load(fn4);
  a=a.p;
  return;
else
  disp(sprintf('ded_read_param: cannot find parameter file for "%s"',n));
  return;
end

try
  i=h5info(fn);
catch
  disp(sprintf('ded_read_param: h5info(''%s'') failed ',fn));
  return;
end
mm={i.Datasets.Name};
[aa dmp]=unix(sprintf('h5dump %s',fn));

for j=1:length(mm)
  m=mm{j};
  %if ~isempty(strmatch(m,{'pmzt'}))
  %  keyboard
  %end
  switch(i.Datasets(j).Datatype.Class)
    case 'H5T_ENUM'
      f = strfind(dmp,m);
      [tk tt]=strtok(dmp(f:end));
      while ~strcmp(tk,'DATA') & length(tt)>0
        [tk tt]=strtok(tt);
      end
      while ~strcmp(tk,'TRUE') & ~strcmp(tk,'FALSE') & length(tt)>0
        [tk tt]=strtok(tt);
      end
      a.(m) = strcmp(tk,'TRUE');
    otherwise
      try
        a.(m)=h5read(fn,['/' m]);
      catch
        disp(sprintf('ded_read_param: Failed to read field %s from %s',m,fn));
        a=[];
        return;
      end
  end  
  if iscell(a.(m)) & length(a.(m))==1
    a.(m)=a.(m){1};
  end
  if isempty(a.(m))
    keyboard
  end
end
if ~isfield(a,'sType')
  a.sType='gc';
end
% $$$ mf={'signals','parallel','forcing','ck_time','Wx','Wy','Wz','x0','x1','x2','x3','Smooth'};
% $$$ for j=1:length(mf)
% $$$   if ~isfield(a,mf{j});
% $$$     a.(mf{j})=NaN;
% $$$   end
% $$$ end
switch(a.sType)
  case 'gc'
    if isfield(a,'U')
      a.U=abs(a.U);
      if ~isfield(a,'V')
        a.V=a.U;
      end
    end
end


[dd nm]=fileparts(fn);
fnv=sprintf('%s/version',dd);
if isfile(fnv)
  fp=fopen(fnv);
  a.version=fgetl(fp);
  fclose(fp);
end

fnt=sprintf('%s/time',dd);
fns=sprintf('%s/status',dd);
fnj=sprintf('%s/jobid',dd);
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

f=findstr(dd,'dedalus');
if isempty(f)
  a.name=dd;
else
  a.name=dd(f+8:end);
end

f=find(a.name=='/');
switch(length(f))
  case 0
    a.series='';
    a.num='';
  case 1
    a.series='';
    a.num=a.name(f(1)+1:end);
  otherwise
    a.series=a.name(f(1)+1:f(2)-1);
    a.num=a.name(f(2)+1:end);
end

if ~isfield(a,'Scb') & isfield(a,'Sc')
  a.Scb=a.Sc;
end

a=orderfields(a);

nms=fieldnames(a);
for j=1:length(nms)
  nm=nms{j};
  if isnumeric(a.(nm));
    a.(nm)=double(a.(nm));
  end
end
if  isfield(a,'AAj'); a.AAJ=a.AAj;a=rmfield(a,'AAj');end;
if ~isfield(a,'AAJ'); a.AAJ=1;end;
if ~isfield(a,'AAS'); a.AAS=1;end;

if ~isfield(a,'hb') & isfield(a,'h')
  a.hb=a.h;
end
if ~isfield(a,'hu') & isfield(a,'h')
  a.hb=a.h;
end

% Rename fields
r={{'Length',   'L'}
   {'Height',   'H'}
   {'Width',    'W'}
   {'Gravity',  'g'}
   {'Buoyancy', 'B'}
   {'Sediment', 'S'}
   {'diffdom',  'ddddiv'}
   {'Velocity', 'V'}};
for j=1:length(r)
  b=r{j}{1};
  c=r{j}{2};
  if isfield(a,b) 
    a.(c)=a.(b);
    a=rmfield(a,b);
    %disp(sprintf('Renaming %s: %s to %s',a.nm,b,c));
  end
end
r={'topdudx','h','wdiffdom'};
a.T={};
if isfield(a,'Tx'); a.T{end+1}=a.Tx;end
if isfield(a,'Ty'); a.T{end+1}=a.Ty;end
if isfield(a,'Tz'); a.T{end+1}=a.Tz;end
if a.Ny>1
  a.X   = [a.L a.W a.H];
else
  a.X   = [a.L a.H];
end

for j=1:length(r)
  if isfield(a,r{j}) 
    a=rmfield(a,r{j});
  end
end

return

b={'dbD','ddD','dbz','dblogh','ddlogh','dbdirect'};
for j=1:length(b)
  if isfield(a,b{j}) 
    for k=1:length(a)
      a(k).(b{j})=logical(a(k).(b{j}));
    end
  end
end
