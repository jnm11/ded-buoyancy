function a=ded_gc_steady(nm,trg,fx,fu,typ,LL)

d=ded_dedalus_data_dir;
dd=[d '/results/' nm];
fnsteady = [dd '/steady.mat'];
fncheb   = [dd '/steady-cheb.mat'];

fns=ded_get_fn(nm,typ);
[tt nmb]=ded_get_times(fns);
[tt f]=sort(tt);

f=f(find(tt>=trg(1) & tt<=trg(2)));
if isempty(f)
  disp(sprintf('ded_gc_steady: No times %s %4.1f %4.1f',nm,trg(1),trg(2)));
  a=[];
  return
end

nmb={nmb{f}};
if any(file_nt(nmb,fncheb)) | ~isfile(fncheb)
  a=ded_mavrg2(nm,trg,fx,fu,typ,LL);
  if ~isempty(a)
    save(fncheb,'a');
  end
end

if file_nt(fncheb,fnsteady) | ~isfile(fnsteady)
  load(fncheb);
  p=ded_read_param(nm);
  b=ded_zgrid(a,2*p.Nz,{},[],[],[],p.H,1);
  b.t1=a.t1;
  b.t2=a.t2;
  a=b;
  save(fnsteady,'a');
end
load(fnsteady);

