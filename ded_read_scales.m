function [a typ tnm stnm]=ded_read_scales(fn,a)
if nargin<2
  a=struct;
end
b=h5info(fn);
if length(b.Groups)>0
  n1g=length(b.Groups(1).Groups);
  k=0;
  for j=1:n1g
    ds=b.Groups(1).Groups(j).Datasets;
    if ~isempty(ds)
      k=k+1;
      n1=b.Groups(1).Groups(j).Name;
      n2=b.Groups(1).Groups(j).Datasets.Name;
      nm=n1(end);
      a.(nm)=h5read(fn,[n1 '/' n2]);
      a.(['n' nm])=length(a.(nm));
      a.sz(k)=length(a.(nm));
    end
  end
  stnm='/scales/sim_time';
else
  typ=setdiff({b.Datasets.Name},{'x','y','z','dt','i1','i2','t1','t2'});
  tnm='';
  stnm='/t1';
end


if length(b.Groups)>1
  tnm=b.Groups(2).Name;
  typ={b.Groups(2).Datasets.Name};
end

if isfield(b.Datasets,'Name')
  if sum(cellstrfind({b.Datasets.Name},'x'))>0
    a.x=h5read(fn,'/x');
  end
  if sum(cellstrfind({b.Datasets.Name},'y'))>0
    a.y=h5read(fn,'/y');
  end
  if sum(cellstrfind({b.Datasets.Name},'z'))>0
    a.z=h5read(fn,'/z');
  end
end

if ~isfield(a,'sz')
  a.sz=[];
  if isfield(a,'x')
    a.sz(end+1)=length(a.x);
  end
  if isfield(a,'y')
    a.sz(end+1)=length(a.y);
  end
  if isfield(a,'z')
    a.sz(end+1)=length(a.z);
  end
end


