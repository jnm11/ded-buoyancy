function a=ded_pm_psd(nm,rgt)

fnpsd=[ded_dedalus_data_dir '/' nm '/psd.mat'];
if isfile(fnpsd)
  disp(sprintf('ded_pm_psd: %s exists',fnpsd));
  load(fnpsd);
  return;
end
a=[];
c=ded_coord(nm);
if isempty(c); return; end
p=ded_read_param(nm);
P=[strcmp(p.Tz,'SinCos') strcmp(p.Ty,'SinCos')];

fns=ded_get_fn(nm,'b');
nb=0;
pb=0;
for j=1:length(fns)
  disp(sprintf('ded_pm_psd: %s',fns{j}));
  b=ded_read_hdf(fns{j});
  if isempty(b); continue; end;
  if b.t<rgt(1) | b.t>rgt(2); continue; end
  for i=1:2
    if P(i)
      f=squeeze(mean(abs(fft(cat(i,b.b,flip(b.b,i)),[],1)).^2,3-i));
    else
      f=squeeze(mean(abs(fft(b.b,[],i)).^2,3-i));
    end
    pb=pb+f(1:end/2,:);          
  end
  nb=nb+1;
end
pb=pb/nb;

fns=ded_get_fn(nm,'ps');
ns=0;
ps=0;
for j=1:length(fns)
  disp(sprintf('ded_pm_psd: %s',fns{j}));
  b=ded_read_hdf(fns{j});
  if isempty(b); continue; end;
  if b.t<rgt(1) | b.t>rgt(2); continue; end
  for i=1:2
    if P(i)
      f=squeeze(mean(abs(fft(cat(i,b.b,flip(b.b,i)),[],1)).^2,3-i));
    else
      f=squeeze(mean(abs(fft(b.b,[],i)).^2,3-i));
    end
    ps=ps+f(1:end/2,:);          
  end
  ns=ns+1;
end
ps=ps/ns;

a.nm=nm;
a.ps=ps;
a.pb=pb;
a.ps=ps;
a.ns=ns;
a.nb=nb;
a.x=c.Jx;
a.k=c.dJy*(0:size(pb,1)-1);
save(fnpsd,'a');
return;
d={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','50','51','54','55','56','57','58'};
for j=1:length(d);
  nm=['pm/f7/e/' d{j}];
  rgt=[-inf inf];
  a=ded_pm_psd(nm,rgt);
end



