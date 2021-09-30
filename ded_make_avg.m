function [fnmat a]=ded_make_avg(nm,typ,ctp,force)

%ded_make_avg(cellstr_ls('gc/f6/g/*/*/*/[23]*/status',[],'dir'),'ay');
if nargin<3; ctp=struct;end;
if nargin<4; force=false;end;

if iscell(nm)
  for j=1:length(nm)
      fnmat{j}=ded_make_avg(nm{j},typ,ctp);
      if nargout>1
        aa=load(fnmat{j});
        a(j)=aa.a;
      end
  end
  return;
end

fnmat = sprintf('%s/results/%s/%s.mat',ded_dedalus_data_dir,nm,typ);
fnay  = ded_get_fn(nm,typ);
if isempty(fnay); 
  disp(sprintf('ded_make_avg: no files "%s" %s',nm,typ));
% $$$   if isfile(fnmat)
% $$$     unix(sprintf('/bin/rm -f %s',fnmat'));
% $$$     disp(sprintf('/bin/rm -f %s',fnmat'));
% $$$   end
   return; 
end
if file_nt(fnmat,fnay{end}) & ~force
  disp(sprintf('ded_make_avg: up to date "%s"',nm));
  return;
end

T=ded_convergence_T(nm,ctp);
trg=[T inf];
if typ(1)=='a'
  a=ded_read_javrg(nm,typ,trg);
else
  a=ded_read_g(nm,typ,'avrg');
end
mkdirifnotexist(fnmat);
save(fnmat,'a');
