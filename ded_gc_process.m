function ded_gc_process

fns=cellstr_ls('~/gc/1*/yavg',[],'dirname');
dd=[fileparts(fileparts(fileparts(fns{1}))) '/'];
fns=cellstrremove(cellstrremoveprefix(fns,dd),'/yavg');
n=length(fns);

for j=1:n
  [T(j) dt(j) b(j)]=ded_convergence_T(fns{j});
  pp=ded_read_param(fns{j});
  if j>1
    df=setxor(fieldnames(pp),fieldnames(p));
    if ~isempty(df)
      disp(df)
      keyboard;
    end
  end
  p(j)=pp;
end

Tr=2000;

X=[b.MB]+sqrt(3)*[b.SB];

fp=fopen('~/misc/alt-gc','w');
fprintf(fp,'cd %s\n',dd);
for j=1:n
  Xmin  = 7;
  Xmax  = p(j).L-0.5;
  nm    = p(j).name;
  t     = p(j).t;
  XT(j) = ((X(j)<Xmin | X(j)>Xmax)  & t>200) | t>Tr+T(j);
  if str2num(nm)<100
    continue;
  end
  switch(p(j).status)
    case {'Running','Aborted','Submitted'}
      if XT(j)
        disp(sprintf('%s %s but X not %.3f<%.3f<%.3f',p(j).status,p(j).name,Xmin,X(j),Xmax));
        s=sprintf('ded_setparam.py %s/param.h5 T %.0f',nm,floor(t));disp(s);fprintf(fp,[s '\n']);
        s=sprintf('echo Terminated > %s/status',nm);disp(s);fprintf(fp,[s '\n']);
        s=sprintf('cat  %s/jobid >> %s/jobid.log',nm,nm);disp(s);fprintf(fp,[s '\n']);
        s=sprintf('/bin/rm -f %s/jobid',nm);disp(s);fprintf(fp,[s '\n']);
        switch(p(j).status)
          case 'Running'
            s=sprintf('touch %s/abort',nm);disp(s);fprintf(fp,[s '\n']);
          case 'Aborted'
          case 'Submitted'
            s=sprintf('scancel %u',p(j).jobid);disp(s);fprintf(fp,[s '\n']);
        end
      else
        switch(p(j).status)
          case 'Running'
          case 'Aborted'
            s=sprintf('dedresub  %s',nm);disp(s);fprintf(fp,[s '\n']);
          case 'Submitted'
        end
      end
    case 'Terminated'
      if ~XT(j)
        if isfinite(T(j))
          nT=T(j)+200;
        else
          nT=t+300;
        end
        nT = 50*ceil(nT/50);
        s=sprintf('ded_setparam.py %s/param.h5 T %.0f',p(j).name,nT);disp(s);fprintf(fp,[s '\n']);
        s=sprintf('echo Aborted > %s/status',nm);disp(s);fprintf(fp,[s '\n']);
        s=sprintf('dedresub  %s',nm);disp(s);fprintf(fp,[s '\n']);
      end
    case 'Suspended'
    otherwise
      disp(sprintf('Unknown status %s %s',nm,p(j).status));
  end
end 
fclose(fp);
% $$$ f=findmax(dt);
% $$$ nm=nms{f};
% $$$ [T dt]=ded_gc_convergence_T(nm);


