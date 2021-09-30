nms=ded_mavrg_nms();

for j=1:length(nms)
  nm=nms{j};
  disp(nm);
  nnm=nm;
  nnm(nm=='/')='-';
  p=ded_read_param(nm);
  b=ded_read_stats(nm);
  dt=max(diff(sort(b.t)));
  [w X t]=jgrid(b.t',b.X',dt,'cubic');
  X=filter_lowpass(X,4);
  
  switch(nm)
    case 'gc/ccle/046'
      t1=5;
      t2=28;
      dt=0.5;
    otherwise
      t1=18;
      t2=42;
      dt=1;
  end
  f=find(t>=t1 & t<=t2);
  if length(f)>4
    while (1)
      px=polyfit(t(f),X(f),2);
      fX=polyval(px,t(f));
      gg=find(abs(fX-X(f))>0.01);
      if isempty(gg)
        break;
      else
        f(gg)=[];
      end
    end
    pv=poly_diff(px,1);
    fx = @(t) polyval(px,t);
    fu = @(t) polyval(pv,t);
    figure;
    subplot(2,1,1);plot(t,X,t,fx(t));axis('tight');ylabel('x');
    subplot(2,1,2);plot(filter_midpoint(t),diff(X)./diff(t),t,fu(t));axis('tight');ylabel('u');
    drawnow;
  end
  tt = dt*(0:1+ceil(max(b.t/dt)))
  for k=1:length(tt)-1
    nmo=sprintf('%s/%s-%02i.mat',ded_dedalus_data_dir,nnm,k);
    disp(nmo);
    if isfile(nmo)
      continue;
    end
    t1=tt(k);
    t2=tt(k+1);;
    if all(b.t<t1)
      continue;
    end
    a=ded_mavrg2(nm,[t1 t2],fx,fu,'y',[-20 5]);
    a=ded_zgrid(a,length(a.z)*2,[],{},{},{},p.H);
    a.t1=t1;
    a.t2=t2;
    a.px=px;
    a.bt=t;
    a.bX=X;
    a.V=fu(t);
    a.X=fx(t);
    a.param=p;
    save(nmo,'a');
    disp(sprintf('Writing %s',nmo));
  end
end

