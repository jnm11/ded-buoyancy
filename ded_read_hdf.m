function a=ded_read_hdf(fn)
%a=ded_read_hdf('force/force_s1.hdf5');

if iscell(fn)
  for j=1:length(fn)
    a(j)=ded_read_hdf(fn{j});
  end
  return;
end
if ~isfile(fn)
  fn2=[ded_dedalus_data_dir '/' fn];
  if ~isfile(fn2)
    disp(sprintf('ded_read_hdf "%s" is not a file',fn));
    a=[];
    return;
  else
    fn=fn2;
  end
end

try
  b=h5info(fn);
catch
  disp(sprintf('ded_read_hdf: "%s" not currently readable as an hdf5 file',fn));
  a=[];
  return;
end
%h5disp(fn);

a=rec([],b,fn);
 
if length(b.Groups)==0
  %  keyboard
  %a=b;
  return;
end

for j=1:length(b.Groups)
  g=b.Groups(j);
  gn=g.Name;
  gnm=gn(gn~='/');
  a.(gnm)=rec(gn,g,fn);
  a=cs(a,gnm,'tasks');
  a=cs(a,gnm,'scales');
end

s='xyz';
for j=1:length(s)
  if isfield(a,s(j))
    if isnumeric(a.(s(j)))
      a.(s(j))=a.(s(j));
    else
      nn=fieldnames(a.(s(j)));
      for jj=1:length(nn)
        if isfloat(a.(s(j)).(nn{jj}))
          break;
        end
      end
      if ~isfloat(a.(s(j)).(nn{jj}))
        disp(a.(s(j)));
        error('ded_read_hdf no float scales');
      end
      a.(s(j))=a.(s(j)).(nn{jj});
    end
  end
end

function a=cs(a,m,n)
if strcmp(m,n)
  b=fieldnames(a.(m));
  for k=1:length(b)
    a.(b{k})=a.(m).(b{k});
  end
  a=rmfield(a,m);
end


function a=rec(nm,g,fn)
for k=1:length(g.Datasets)
  d=g.Datasets(k).Name;
  [nd nr]=sscanf(d,'%f');
  if nr==0
    dd=fxnm(d);
  else
    if length(g.Datasets(k).Attributes)>2
      dd=g.Datasets(k).Attributes(2).Value(1);
      a.(dd).scale=nd;
    else
      disp(sprintf('ded_read_hdf: failed to read "%s"',fn));
      disp(b);
      a=[];
    end
  end
  a.(dd)=h5read(fn,[nm '/' d]);
  A=g.Datasets(k).Attributes;
  for j=1:length(A)
    anm=[dd '_' A(j).Name];
    a.(anm)=A(j).Value;
  end
end
for k=1:length(g.Groups)
  gn=g.Groups(k).Name;
  d=fxnm(gn(max(find(gn=='/'))+1:end));
  aa=rec(gn,g.Groups(k),fn);
  if ~isempty(aa)
    a.(d)=aa;
  end
end
if ~exist('a','var')
  a=[];
end

function d=fxnm(c)
if(any(c(1)=='0123456789'))
  d=['a' c];
else
  d=c;
end

return

