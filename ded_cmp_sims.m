function a=ded_cmp_sims(nm,nmr)
% Compare a series of simulations
%ded_cmp_sims('gc/5*');
if nargin<2
  nmr={'PIDDD','PIDIT','PIDS1','PIDS2','PIDST','PIDT','ubc','lbc','rbc','gmax','gmin','Umax','Umin'...
       'blbc','brbc','rbc','slbc','srbc','AAj','ck_time','time',...
       'dtb','dtcheck','dtforce','dtjWx','dtjWy','dtjWz','dtja','dtjavrg','dtjayz',...
       'dtjb','dtjcheck','dtjd','dtjp','dtjr','dtjs','dtju','dtjv','dtjw','dtjy','dtxyz','dtyz','m1','m2'};
end
fns=ded_parse_nms(nm);
if length(fns)==0
  a=[];
  disp(['ded_cmp_sims.m: no matching simulation ' m]);
  return;
end
if length(fns)==1
  a=[];
  disp(['ded_cmp_sims.m: only one matching simulation ' fns{1}]);
  return;
end

a=ded_read_param(fns);
if isempty(a)
  return;
end
fna=fieldnames(a);
f=intersect(fna,{'cm','PIDG','signals','parallel','avrg','fpidu','fpidv','fpidw','fpidb','fpids','clipB','clipS','fbmult','fsmult','fbmin','fsmin','fbmax','fsmax','Conservative','db','oddg','dtjadv','dtjsw','dtjd1','dtjd2','pmzr','pmzt','pmzrt','topdrt','topd','topr','topdd','topu','ful','fur','wul','wvl','wur','wvr','wwl','wwr','rescale','fixI','pmIu','ddiv','ddivu','ddivq','ddivscl','intdx','divxp','divxI','topq','topdivr','ddddiv'});
for j=1:length(f)
  for jj=1:length(a)
    if isfinite(a(jj).(f{j}))
      a(jj).(f{j})=logical(a(jj).(f{j}));
    end
  end
end

f=intersect(fna,{'B','S','Skg','Ski','Skp','pmsl','pmss','SV'});
for j=1:length(f)
  ff=~isfinite([a.(f{j})]);
  if all(ff)
    a=rmfield(a,f{j});
  else
    ff=find(ff);
    for jj=1:length(ff)
      a(ff(jj)).(f{j})=0;
    end
  end
end

f=intersect(fieldnames(a),{'Scs','Scb','Scd'});
for j=1:length(f)
  ff=~isfinite([a.(f{j})]);
  if all(ff)
    a=rmfield(a,f{j});
  else
    ff=find(ff);
    for jj=1:length(ff)
      a(ff(jj)).(f{j})=1;
    end
  end
end



s=ded_read_stats({a.name});
for j=1:length(s)
  a(j).Xrms = 0;
  a(j).gstd = 0;
  a(j).Urms = 0;  
  
  n=max(round(length(s(j).t)/2),min(find(s(j).t)>=s(j).t(end)-max(20,s(j).t(end)*0.8)));
  if isfield(s,'g')
    if all(isfinite(s(j).g))
      a(j).gstd = std(s(j).g(n:end));
      a(j).g    = mean(s(j).g(n:end));
    end
  end
  if isfield(s,'X')
    if all(isfinite(s(j).X))
      a(j).Xrms = std(s(j).X(n:end));
      a(j).X    = mean(s(j).X(n:end)); 
      if isfield(a,'PIDX')
        if isfinite(a(j).PIDX)
          a(j).Xrms = sqrt(mean((s(j).X(n:end)-a(j).PIDX).^2));
          a(j).X = a(j).X -a(j).PIDX;
        end
      end
    end
  end
  if isfield(s,'U')
    if all(isfinite(s(j).U))
      a(j).Ustd = std(s(j).U(n:end));
      a(j).U    = mean(s(j).U(n:end));
    end
  end
  if ~isfinite(a(j).Xrms)
    a(j).Xrms=0;
  end
end

fff={'b','s','d'};
for j=1:length(fff)
  ffS=['Sc' fff{j}];
  ffP=['Pe' fff{j}];
  if isfield(a,ffS)
    for j=1:length(a)
      a(j).(ffP)=round(a(j).Re*a(j).(ffS)*1000)/1000;
    end
    a=rmfield(a,ffS);
  end
end

nmr=intersect(fieldnames(a),nmr);
a=rmfield(a,nmr);

[x yy]=struct_cmp_array(a);

n=length(x);
fprintf('Differences\n');
nm=fieldnames(x);
nm=setdiff(nm,{'jobid','name'});
nm={'name',nm{:}};
m=length(nm);
nmst=intersect(nm,{'Re','g','Scb','Nx'});
xx=zeros(length(x),length(nmst));
for k=1:length(nmst)
  xx(:,k)=[x.(nmst{k})];
end
%[xx f]=sortrows(xx);
%[xx f]=sortrows({x.name});
%x=x(f);

for k=1:m
  for j=1:length(x)
    y=x(j).(nm{k});
    if isfinite(y) & ~isempty(y);break;end
  end
  if islogical(y)
    for kk=1:length(x)
      if x(kk).(nm{k})==0
        x(kk).(nm{k})='F';
      elseif x(kk).(nm{k})==1
        x(kk).(nm{k})='T';
      else
        x(kk).(nm{k})='';
      end
    end
    y='T';
    maxx=1;
  elseif isnumeric(y)
    xx=sort([x(:).(nm{k})]);
    xx=xx(isfinite(xx));
    if length(xx)>1
      dx=diff(xx);
    else
      dx=xx/10;
    end
    dx(dx==0)=1;
    dx=ceil(min(4,max(0,1-log10(min(dx)))));
    maxx=log10(max(abs(xx)));
    maxx=sign(maxx)*ceil(0.01+abs(maxx));
    dx=max(dx,-maxx);
    maxx=max(0,maxx)+dx+2*(dx>0);
  elseif ischar(y)
    xx=[];
    for j=1:length(x)
      xx(j)=length(x(j).(nm{k}));
    end
    maxx=max(xx);
  end
  maxx=1+max([maxx length(nm{k})]);
  if ischar(y)
    fms{k}=sprintf('%ds',maxx);
  elseif isnumeric(y)
    fms{k}=sprintf('%d.%df',maxx,dx);
  else
    fms{k}='6.1g';
  end
  switch(nm{k})
    case {'Nx','Ny','Nz'}
      fms{k}='4.0f';
    case {'Ski','Skg','Skp'}
      fms{k}='7.5f';
    case {'PIDI'}
      fms{k}='6.4f';
    case {'PIDP'}
      fms{k}='6.3f';
    case {'PIDD'}
      fms{k}='6.2f';
     case {'Wx','Wy','Wz','hu','hb','wu','wb'}
      fms{k}='4.2f';
     case {'U1'}
      fms{k}='4.2f';
     case {'Angle','PIDX','dtjyz','dtjay','dbT','ddT'}
      fms{k}='5.1f';
    case {'dblogT','ddlogT'}
      fms{k}='6.1f';
    case {'B'}
      fms{k}='6.4f';
    case {'S'}
      fms{k}='+7.4f';
    case {'g'}
      fms{k}='6.4f';
     case {'dUX'}
      fms{k}='4.1f';
   case {'Qu'}
      fms{k}='5.3f';
  end
  fms{k}=['%' fms{k} ' '];
  fw{k}=length(sprintf(fms{k},y))-1;
  
  fmt{k}=['%' num2str(fw{k}) 's '];
  fprintf(fmt{k},nm{k});
end
fprintf('\n');
%disp(fms);
%disp(fmt);

for j=1:n
  for k=1:m
    if ~isempty(x(j).(nm{k})) & isfinite(x(j).(nm{k}))
      fprintf(fms{k},x(j).(nm{k}))
    else
      f=min(find(fms{k}>='0' & fms{k}<='9'));
      bw=str2double(fms{k}(f));
      fprintf(repmat(' ',1,1+bw));
    end
  end
  fprintf('\n');
end

return

%nmr={'PIDT','PIDW','Sc','Smooth','T','Tz','U1','VMEM','W','Wx','Wy','Wz','XS','XT','ck_time','dtavrg','dtb','dtgc','dtleft','dtmomb','dtnoise','dtp','dtpm','dtright','dtu','dtv','dtw','dtx','dtxy','dtxyz','dtxz','dty','dtyz','dtz','g','h','noiseL','noiseT','parallel','q','signals','status','x0','x1','x2','x3','x4','x5','x6'};
%nmr={'PIDDD','PIDIT'};
%nmr=intersect(fieldnames(y),nmr);
%y=rmfield(y,nmr);
%disp('');
%disp('Similarities');
%disp(yy);

  nmr={'PIDDD','PIDIT','PIDS1','PIDS2','PIDST','PIDT','ubc','lbc','rbc','gmax','gmin','Umax','Umin'...
       'blbc','brbc','rbc','slbc','srbc','AAj','ck_time','time',...
       'dtb','dtcheck','dtforce','dtjWx','dtjWy','dtjWz','dtja','dtjavrg','dtjayz',...
       'dtjb','dtjcheck','dtjd','dtjp','dtjr','dtjs','dtju', ...
       'dtjv','dtjw','dtjy','dtxyz','dtyz','m1','m2',...
       'AB','PIDD','PIDG','PIDI','PIDP','PIDX',...
      'xa','xb','x0','x1','x2','x3','x4','x5','x6','x7',...
       'clipB','clipS','clipu','dnoise','dtavrg','dtfpid','dtgc','dtjay', ...
       'dtjxyz','dtjyz','dtleft','dtmomb','dtnoise','dtp','dtpm', ...
       'dtright', 'dtslice','dtstats','dtu','dtv','dtw','dtx', ...
       'dtxy','dtxz', 'dty','dtz','force','forcing','fpid', ...
       'fpidb','fpids', 'fpidu','fpidv','fpidw','hu','inlet','inoise', ...
       'maxdt','mindt','noise','noiseL','noiseT','noised','num', ...
       'q','radius', 'scheme','version','xn1','xn2','xn3','xn4',...
       'Time','U','V','VMEM','W','Wx','Wy','Wz','alpha','avrg','bV','Angle','Peb','AA'};
  ded_cmp_sims('gc/ccle/*',nmr); 
  