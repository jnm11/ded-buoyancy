function nms=ded_parse_nms(nms)
if ~iscell(nms)
  nms={nms};%ded_regexp_names(nms);
end
dd=[ded_dedalus_data_dir '/'];
for j=1:length(nms)
  fns{j}=wordexp([dd nms{j} '/param.h5']);
end
nms=sort(cat(1,fns{:}));
nms=cellstrremoveprefix(nms,dd);
nms=cellstrremovesuffix(nms,'/param.h5');

