function nma=ded_regexp_names(nm)
DD=[ded_dedalus_data_dir '/'];
if ~iscell(nm)
  if any(nm=='*' | nm=='?')
    nms=cellstr_ls([DD nm ],'ls -d');
  elseif ~isfile(nm)
    nms={[DD nm]};
  elseif ~isfile(nm)
    nms={fileparts(nm)};
  else
    disp(sprintf('ded_regexp_names: "%s" is not a directory, file or regexp',nm{j}));
  end
end
k=0;
nma={};
for i=1:length(nms)
  if isdir(nms{i})
    k=k+1;
    nma{k}=nms{i};
  else
    disp(sprintf('ded_regexp_names: "%s" is not a directory',nms{i}));
  end
end
if isempty(nma)
  disp(sprintf('ded_regexp_names: No simulations found for "%s"',nm));
end
