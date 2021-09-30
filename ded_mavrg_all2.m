nms=ded_mavrg_nms();

for j=1:length(nms)
  nm=nms{j};
  fno=sprintf('%s/%s.mat',ded_dedalus_data_dir,nm);
  if isfile(fno)
    continue;
  end
  a=[];
  N=0;
  switch(nm)
    case 'gc/ccle/046'
      t1=5;
      t2=18;
    otherwise
      t1=18;
      t2=42;
  end
  
  fns=cellstr_ls(sprintf('%s/%s-*.mat',ded_dedalus_data_dir,nm));
  a=[]
  for k=1:length(fns)
    b=load(fns{k});
    if a.t1<t1 | a.t2>t2
      continue;
    end
    if isempty(a)
      a=b.a;
      fnm=fieldnames(a);
      f=[];
      for jj=1:length(fnm)
        if ndims(a.(fnm{jj}))==2
          f(end+1)=jj;
        end
      end
      fnm={fnm{f}};
      fnm=celldiff(fnm,{'T'});
      continue;
    end
    for jj=1:length(fnm)
      a.(fnm{jj})=a.(fnm{jj})+b.a.(fnm{jj})*b.a.T;
    end
    a.T=a.T+b.a.T;
    a.t1=min(a.t1,b.a.t1);
    a.t2=max(a.t2,b.a.t2);
  end
  if isempty(a)
    continue;
  end
  for k=1:length(fnm)
    a.(fnm{k})=a.(fnm{k})./a.T;
  end
  save(fno,'a');
end


