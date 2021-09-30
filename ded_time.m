function [t,fns]=ded_time(nm,typ,display) 
%ded_time(nm,'b');
if nargin<3
  display=0;
end
nmb=ded_get_fn(nm,typ,1);
kk=[];
tn =1;
n=length(nmb);
f=[];
for k=1:length(nmb)
  fn=nmb{k};
  try
    ttt=h5read(fn,'/scales/sim_time');
    f(end+1)=k;
  catch
    disp(fn);
    continue
  end
  tt{k}=ttt;
  kk(end+(1:length(tt{k})))=k;
  if length(tt{k})~=1 
    tn=0;
  end
end
t=[tt{f}];
nmb={nmb{f}};

p=ded_read_param(nm);
s=round(t/p.dtb);
e = abs(t/p.dtb-s);
f = find(e>0.005);
nt=1:length(t);
plot(nt,t,'s',nt,s*p.dtb);

if ~isempty(f)
  r1={nmb{f}};
  r2=cellstrprefix([ded_dedalus_data_dir 'gc', cellstrremoveprefix(r1,[ded_dedalus_data_dir 'pgc'));
  for j=1:length(f)
    cmd=sprintf('/bin/rm -f %s',r1{j});
    unix(cmd);
    cmd=sprintf('/bin/rm -f %s',r2{j});
    disp(cmd);   
  end
end


return


if any(diff(t)<=0)
  disp('Non increasing t');
  %  disp(t);
  if tn==1
    f=min(find(diff(t)<=0));
    f=find(t(1:f)>=t(f+1));
    for j=1:length(f)
      k=f(j);
      disp(sprintf('/bin/rm -f %s %7.2f',nmb{k},t(k)));
      unix(sprintf('/bin/rm -f %s\n',nmb{k}));
      unix(sprintf('echo /bin/rm -f %s >> ~/misc/ded-3d-rm',nmb{k}));
    end
    t(f)=[];
    nmb(f)=[];
  else
    %keyboard;
  end
  if p.cont
    unix(sprintf('/bin/rm -f %s/*.png',dd'));
    p.cont=0;
  end
  disp('ded_3d: should be rerun');
  return;
end

return

