function ded_find_non_monotonic

fns=wordexp('~/data/dedalus/gc/*/xyz');
f=findstr(fns{1},'dedalus/');
fns=cellstrremoveprefix(fns,fns{1}(1:f+7));
fns=cellstrremove(fns,'/xyz');

for j=1:length(fns)
  a=ded_read_g(fns{j},'xyz');
  if ~isfield(a,'t')
    continue
  end
  f=find(diff(a.t)<0);
  if isempty(f)
    continue;
  end
  fprintf(1,'ded_find_non_monotonic: %s',fns{j});
  disp(['  ' sprintf('%6.2f ',a.t(f))]);
end

