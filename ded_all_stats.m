function [a A nms]=ded_all_stats(dd)


if iscell(dd)
  fns=dd;
else
  fns=wordexp(['~/data/dedalus/' dd '/param.h5']);
  if length(fns)==1
    fn2=wordexp(['~/data/dedalus/' dd '/*/param.h5']);
    if length(fn2)>1
      fns=fn2;
    end
  end
end
if isempty(fns)
  disp('ded_all_stats: No files found');
  a=[];A=[];nms={};
  return;
end
if ~isfile(fns{1})
  disp('ded_all_stats: No files found');
  a=[];A=[];nms={};
  return;
end

fns=cellstrremove(fns,'/param.h5');
n=length(fns);
kk=zeros(1,n);
k=0;
for jj=1:3
  for j=1:n
    if kk(j)==1
      break;
    end
    try
      aa=ded_stats(fns{j});
      kk(j)=1;
    catch
      disp(sprintf('ded_all_stats: failed %s',fns{j}));
      continue
    end
    if ~isempty(aa)
      k=k+1;
      nms{k}=fns{j};
      if k==1
        a=orderfields(aa);
      else
        a=struct_array_append(a,aa);
      end
    end
  end
end

A=[[a.forcing]' [a.Nx]'  [a.L]'  [a.H]' [a.W]' [a.U]' [a.U1]' [a.Re]' [a.PIDX]'];
[A f]=sortrows(round(1e5*A));
a=a(f);
nms={nms{f}};
A=A/1e5;
