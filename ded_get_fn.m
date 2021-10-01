function fns=ded_get_fn(n,m,w,mm)
if nargin<3
  w=0;
end
if nargin<4
  mm=[];
end
if isempty(mm)
  mm=m;
end
if isdir(n)
  fsa=sprintf('%s/%s/%s*.h*5',n,m,mm);
  fs2=sprintf('%s/%s.h5',     n,mm);
  fs3=sprintf('%s/%s.hdf5',   n,mm);
else
  fsa=sprintf('%s/%s/%s/%s*.h*5',ded_dedalus_data_dir,n,m,mm);
  fs2=sprintf('%s/%s/%s.h5',     ded_dedalus_data_dir,n,mm);
  fs3=sprintf('%s/%s/%s.hdf5',   ded_dedalus_data_dir,n,mm);
end
if isfile(fs2)
  fns={fs2};
elseif isfile(fs3)
  fns={fs3};
else
  fns=cellstr_ls(fsa);
end
if isempty(fns)
  if w
    disp(sprintf('ded_get_fn: No files of type "%s" found for "%s"',m,n));
  end
  fns={};
  return;
end
nfn=length(fns);
k=zeros(nfn,1);
if nfn>1
  for i=1:nfn
    [dd nm]=fileparts(fns{i});
    kk=str2num(nm(length(m)+3:end));
    if isempty(kk)
      kk=i;
    end
    k(i)=kk;
  end
  [k f]=sort(k);
  fns={fns{f}};
end

return;

% $$$ 
% $$$ while 1
% $$$   fn1=sprintf('%s/%s/%s_s%i.h5',n,m,m,k);
% $$$   fn2=sprintf('%s/%s/%s_s%i/%s_s%i_p0.h5',n,m,m,k,m,k);
% $$$   fn3=sprintf('%s/%s/%s_s%i.hdf5',n,m,m,k);
% $$$   fn4=sprintf('%s/%s/%s_s%i/%s_s%i_p0.hdf5',n,m,m,k,m,k);
% $$$   if isfile(fn1)
% $$$     fns{j}=fn1;
% $$$     j=j+1;
% $$$   elseif isfile(fn2)
% $$$     fns{j}=fn2;
% $$$     j=j+1;
% $$$   elseif isfile(fn3)
% $$$     fns{j}=fn3;
% $$$     j=j+1;
% $$$   elseif isfile(fn4)
% $$$     j=j+1;
% $$$     fns{j}=fn4;
% $$$   end
% $$$   k=k+1;
% $$$   if k>j+50
% $$$     break;
% $$$   end
% $$$ end
% $$$ if j==1
% $$$   fns={};
% $$$ end
