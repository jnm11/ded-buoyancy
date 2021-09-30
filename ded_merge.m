function a=ded_merge(a)
if length(a)>1
  s={'x','z','y','nm'};
  flnm=setdiff(fieldnames(a),s);
  for j=1:length(flnm)
    nd=1;
    for k=1:length(a)
      nd=max([nd,find(size(a(k).(flnm{j}))>1)]);
    end
    aa.(flnm{j})=cat(nd,a.(flnm{j}));
  end
  for j=1:length(s)
    if isfield(a,s{j})
      aa.(s{j})=a(1).(s{j});
    end
  end
  a=aa;
end
return
