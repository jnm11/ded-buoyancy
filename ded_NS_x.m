function a=ded_NS_x(nm,typ,trg)


if typ(1)=='a'
  a=ded_read_javrg(nm,typ,trg);
elseif strcmp(typ,'final')
  fns=ded_get_fn(nm,'final',[],'state');
  a=ded_read_hdf(fns{end});
end
parity=ded_read_parity(nm);
coord=ded_coord(nm);
param=ded_read_param(nm);
a=ded_add_quad(a,{'bu','bw','uu','uw','ww'});
a=ded_add_diff(a,param,parity,coord,{'b','u','w','p','uu','ww','uw','bu','bw'});
a=ded_add_quad(a,{'ududx','wdudz','udwdx','wdwdz','udbdx','wdbdz'});
a=ded_add_diff(a,param,parity,coord,{'dudx','dudz','dwdx','dwdz','dbdx','dbdz'});

nn=setdiff({'ududx','wdudz','udwdx','wdwdz'},fieldnames(a));
for j=1:length(nn)
  n=nn{j};
  n1=n(1);
  n2=n(2:end);
  a.(n)=a.(n1).*a.(n2);
end
nz=size(a.u,1);
w=ichebintw(nz)/2;
[w1 w2]=ichebendsw(nz);

nn={'udbdx','wdbdz','dbudx','dbwdz','duudx','duwdz','ududx','wdudz','duwdx','dwwdz','udwdx','wdwdz','dpdx','dpdz','ddudxdx','ddudzdz','ddwdxdx','ddwdzdz','ddbdxdx','ddbdzdz','b'};
a.b=a.b*param.g;
a.ddudxdx=a.ddudxdx/param.Re;
a.ddudzdz=a.ddudzdz/param.Re;
a.ddwdxdx=a.ddwdxdx/param.Re;
a.ddwdzdz=a.ddwdzdz/param.Re;
a.ddbdxdx=a.ddbdxdx/param.Re/param.Scb;
a.ddbdzdz=a.ddbdzdz/param.Re/param.Scb;
z=coord.Jz/param.H;
for j=1:length(nn)
  b.rms.(nn{j})=sqrt(w*a.(nn{j}).^2);
  b.w1.(nn{j})=w1*a.(nn{j});
  b.w2.(nn{j})=w2*a.(nn{j});
  b.m0.(nn{j})=w *a.(nn{j});
  b.m1.(nn{j})=w*(z.^1.*a.(nn{j}));
  b.m2.(nn{j})=w*(z.^2.*a.(nn{j}));
end

f=find(coord.Jx>1 & coord.Jx<30);
x=coord.Jx(f);
tt={'rms','m0','m1','m2','w1','w2'};
nb={'dbudx','dbwdz','udbdx','wdbdz','ddbdxdx','ddbdzdz'};
nx={'duudx','duwdz',        'wdudz','ddudxdx','ddudzdz'    , 'dpdx'};
nz={'duwdx','dwwdz','udwdx',        'ddwdxdx','ddwdzdz','b'} 'dpdz';

dc;
for i=1:3
  switch(i)
    case(1)
      nn=nx;
      T='u';
    case(2)
      nn=nz;
      T='v';
    case(3)
      nn=nb;
      T='b';
  end
  for k=1:length(tt)
    mmax=0;
    mmin=0;
    figure;
    clf;
    ah=jsubplot([1 length(nn)],[0.1 0.05],[0.01 0.01],[0.01 0.05]);
    for j=1:length(nn)
      axes(ah(j))
      plot(x,b.(tt{k}).(nn{j})(f));
      axis('tight');
      mmax=max(mmax,max(b.(tt{k}).(nn{j})(f)));
      mmin=min(mmin,min(b.(tt{k}).(nn{j})(f)));
      ylabel(nn{j});
    end
    axes(ah(1));
    title([T ' ' tt{k}]);
    set(ah,'ylim',[mmin mmax]);
    set(ah(2:end),'xticklabels',[]);
  end
end
