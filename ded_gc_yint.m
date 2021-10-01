function a=ded_gc_yint(nm,dd)

fntime   = [dd '/time.mat'];

fns=ded_get_fn(nm,'y');
if isfile(fntime)
  f=find(file_nt(fns,fntime));
  if isempty(f)
    load(fntime);
    return;
  fns={fns{f}};
  end
end

p=ded_read_param(nm);

W=p.W;

if isfile(fntime)
  load(fntime)
else
  a.t=[];
end
for k=1:length(fns)
  b=ded_read_hdf(fns{k});
  if isempty(b)
    continue;
  end
  w=ichebintw(length(b.z),p.H);
  dx=b.x(2)-b.x(1);
  wz = b.z'.*w;  
  if isfield(b,'b' ); b.b  = squeeze(b.b/W);  end
  if isfield(b,'bb'); b.bb = squeeze(b.bb/W); end
  if isfield(b,'bu'); b.bu = squeeze(b.bu/W); end
  if isfield(b,'bv'); b.bv = squeeze(b.bv/W); end
  if isfield(b,'bw'); b.bw = squeeze(b.bw/W); end
  if isfield(b,'p');  b.p  = squeeze(b.p/W);  end
  if isfield(b,'pp'); b.pp = squeeze(b.pp/W); end
  if isfield(b,'u');  b.u  = squeeze(b.u/W);  end
  if isfield(b,'uu'); b.uu = squeeze(b.uu/W); end
  if isfield(b,'uv'); b.uv = squeeze(b.uv/W); end
  if isfield(b,'uw'); b.uw = squeeze(b.uw/W); end
  if isfield(b,'v');  b.v  = squeeze(b.v/W);  end
  if isfield(b,'vv'); b.vv = squeeze(b.vv/W); end
  if isfield(b,'vw'); b.vw = squeeze(b.vw/W); end
  if isfield(b,'w');  b.w  = squeeze(b.w/W);  end
  if isfield(b,'ww'); b.ww = squeeze(b.ww/W); end
  
  if isfield(b,'sim_time') b.t=b.sim_time; end;
  if isfield(b,'t1')       b.t=b.t1; end;
  
  for i=1:length(b.t)
    disp(sprintf('ded_gc_yint: %s %4u %4u %7.3f',nm,i,k,b.t(i)))
    j=find(a.t==b.t(i));
    if isempty(j)
      j=length(a.t)+1;
    end

    if isfield(b,'u');  a.u(j)=dx*sum(w*b.u(:,:,i)); end
    if isfield(b,'v');  a.v(j)=dx*sum(w*b.v(:,:,i)); end
    if isfield(b,'w');  a.w(j)=dx*sum(w*b.w(:,:,i)); end
    
    if isfield(b,'uu'); a.uu(j)=dx*sum(w*b.uu(:,:,i)); end
    if isfield(b,'vv'); a.vv(j)=dx*sum(w*b.vv(:,:,i)); end
    if isfield(b,'ww'); a.ww(j)=dx*sum(w*b.ww(:,:,i)); end
    if isfield(b,'uv'); a.uv(j)=dx*sum(w*b.uv(:,:,i)); end
    if isfield(b,'uw'); a.uw(j)=dx*sum(w*b.uw(:,:,i)); end
    if isfield(b,'vw'); a.vw(j)=dx*sum(w*b.vw(:,:,i)); end
    
    if isfield(b,'u') & isfield(b,'u');  a.Auu(j)=dx*sum(w*(b.u(:,:,i).*b.u(:,:,i))); end
    if isfield(b,'v') & isfield(b,'v');  a.Avv(j)=dx*sum(w*(b.v(:,:,i).*b.v(:,:,i))); end
    if isfield(b,'w') & isfield(b,'w');  a.Aww(j)=dx*sum(w*(b.w(:,:,i).*b.w(:,:,i))); end
    if isfield(b,'u') & isfield(b,'v');  a.Auv(j)=dx*sum(w*(b.u(:,:,i).*b.v(:,:,i))); end
    if isfield(b,'u') & isfield(b,'w');  a.Auw(j)=dx*sum(w*(b.u(:,:,i).*b.w(:,:,i))); end
    if isfield(b,'v') & isfield(b,'w');  a.Avw(j)=dx*sum(w*(b.v(:,:,i).*b.w(:,:,i))); end
    
    if isfield(b,'u');  a.bu(j)=dx*sum(w*(b.b(:,:,i).*b.u(:,:,i))); end
    if isfield(b,'v');  a.bv(j)=dx*sum(w*(b.b(:,:,i).*b.v(:,:,i))); end
    if isfield(b,'w');  a.bw(j)=dx*sum(w*(b.b(:,:,i).*b.w(:,:,i))); end
    
    if isfield(b,'uu'); a.buu(j)=dx*sum(w*(b.b(:,:,i).*b.uu(:,:,i))); end
    if isfield(b,'vv'); a.bvv(j)=dx*sum(w*(b.b(:,:,i).*b.vv(:,:,i))); end
    if isfield(b,'ww'); a.bww(j)=dx*sum(w*(b.b(:,:,i).*b.ww(:,:,i))); end
    if isfield(b,'uv'); a.buv(j)=dx*sum(w*(b.b(:,:,i).*b.uv(:,:,i))); end
    if isfield(b,'uw'); a.buw(j)=dx*sum(w*(b.b(:,:,i).*b.uw(:,:,i))); end
    if isfield(b,'vw'); a.bvw(j)=dx*sum(w*(b.b(:,:,i).*b.vw(:,:,i))); end
    
    a.b(j)  =  dx*sum(w*b.b(:,:,i));
    a.bz(j) =  dx*sum(wz*b.b(:,:,i));
    a.E1(j) = -dx*sum(w*(max(0,b.b(:,:,i)).*log(max(realmin,b.b(:,:,i)))));
    a.E2(j) =  dx*sum(w*(max(0,b.b(:,:,i).*(1-b.b(:,:,i)))));
    if a.E1(j)>1e3
      keyboard;
    end
    if isfield(b,'Eb')
      a.E3(j) =  dx*sum(w*b.Eb);
    else
      a.E3(j) = NaN;
    end
    a.t(j)=b.t(i);
    [ddd fn]=fileparts(fns{k});
    a.nm{j}=fn;
  end
end
save(fntime,'a');

