fnmat='~/mat-ccle-046';
nm='gc/ccle/046';

fntime=[fnmat '/046-time.mat'];

p=ded_read_param(nm);
W  = p.W;
H  = p.H;
Nz = p.Nz;

save([fnmat '/046-param.mat'],'p');
s=ded_read_stats(nm);
save([fnmat '/046-stats.mat'],'s');

dt=max(diff(sort(s.t)));
[w X t]=jgrid(s.t',s.X',dt,'cubic');
X([1 end])=[];
t([1 end])=[];

t1=5;
t2=33.5;
f= find(t>=t1 & t<=t2);
p=polyfit(t(f),X(f),2);
dp=p(1:2).*[2 1];


mkdir(fnmat);


fx= @(t) polyval(p,t);
fu= @(t) polyval(dp,t);

save([fnmat '/046-traj.mat'],'t','X','fx','fu');

LL=[-20 5];
trg=[16 33.5];
%trg=[28 28.1];
typ='y';
a=ded_mavrg2(nm,trg,fx,fu,typ,LL);
save([fnmat '/046-steady-cheb.mat'],'a');
b=ded_zgrid(a,2*Nz,{},[],[],[],H,1);
b.t1=a.t1;
b.t2=a.t2;
a=b;
save([fnmat '/046-steady.mat'],'a');





if isfile(fntime)
  load(fntime);
else
  clear('b');
  w  = ichebintw(Nz,H);
  wz = a.z'.*w;
  
  fns=ded_get_fn(nm,'y');
  for k=1:length(fns)
    disp(sprintf('%u/%u',k,length(fns)));
    a=ded_read_hdf(fns{k});
    
    
    dx=a.x(2)-a.x(1);
    
    a.b  = a.b/W;
    a.bb = a.bb/W;
    a.bu = a.bu/W;
    a.bv = a.bv/W;
    a.bw = a.bw/W;
    a.p  = a.p/W;
    a.pp = a.pp/W;
    a.u  = a.u/W;
    a.uu = a.uu/W;
    a.uv = a.uv/W;
    a.uw = a.uw/W;
    a.v  = a.v/W;
    a.vv = a.vv/W;
    a.vw = a.vw/W;
    a.w  = a.w/W;
    a.ww = a.ww/W;
    
    b.u(k)=dx*sum(w*a.u);
    b.w(k)=dx*sum(w*a.w);
    b.v(k)=dx*sum(w*a.v);
    
    b.uu(k)=dx*sum(w*a.uu);
    b.vv(k)=dx*sum(w*a.vv);
    b.ww(k)=dx*sum(w*a.ww);
    b.uv(k)=dx*sum(w*a.uv);
    b.uw(k)=dx*sum(w*a.uw);
    b.vw(k)=dx*sum(w*a.vw);
    
    b.Auu(k)=dx*sum(w*(a.u.*a.u));
    b.Avv(k)=dx*sum(w*(a.v.*a.v));
    b.Aww(k)=dx*sum(w*(a.w.*a.w));
    b.Auv(k)=dx*sum(w*(a.u.*a.v));
    b.Auw(k)=dx*sum(w*(a.u.*a.w));
    b.Avw(k)=dx*sum(w*(a.v.*a.w));
    
    b.bu(k)=dx*sum(w*(a.b.*a.u));
    b.bv(k)=dx*sum(w*(a.b.*a.v));
    b.bw(k)=dx*sum(w*(a.b.*a.w));
    
    b.buu(k)=dx*sum(w*(a.b.*a.uu));
    b.bvv(k)=dx*sum(w*(a.b.*a.vv));
    b.bww(k)=dx*sum(w*(a.b.*a.ww));
    b.buv(k)=dx*sum(w*(a.b.*a.uv));
    b.buw(k)=dx*sum(w*(a.b.*a.uw));
    b.bvw(k)=dx*sum(w*(a.b.*a.vw));
    
    b.b(k)  = dx*sum(w*a.b);
    b.bz(k) = dx*sum(wz*a.b);
    b.t(k)=a.t1;
    [ddd fn]=fileparts(fns{k});
    b.nm{k}=fn;
  end
  a=b;
  save(fntime,'a');
end
