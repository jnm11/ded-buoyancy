function [t nmc dt]=ded_get_times(nmb)

if ~iscell(nmb);
  nmb={nmb};
end	

t=[];
dt=[];
tn =1;
jj=[];
nmc={};
for k=1:length(nmb)
  fn=nmb{k};
  %b=h5info(fn);
  tt=[];
  tt2=NaN;
  if isempty(tt)
    try
      tt=h5read(fn,'/t1');
      tt2=h5read(fn,'/t2');
    end
  end
  if isempty(tt)
    try
      tt=h5read(fn,'/t');
    end
  end
  if isempty(tt)
    try
      tt=h5read(fn,'/scales/sim_time');
    end
  end
  if length(tt)>0
    t=[t;double(tt)];
    dt=[dt;double(tt2)-double(tt)];
    jj(end+1)=k;
    for kkk=1:length(tt)
      nmc{end+1}=fn;
    end
  else
    disp(sprintf('ded_fix_time: %s contains no times',fn));
  end
end
nmb={nmb{jj}};

