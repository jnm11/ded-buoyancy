function [b c]=ded_gc_f7_g_classify(fns,xmin,display)
j=0;
for jj=1:length(fns)
  nm=fns{jj};
  
  fn=[ded_dedalus_data_dir '/results/' nm];
  if isfile(fn)
    disp(sprintf('ded_gc_f7_g_classify: loading %s',fn));
    pr=load(fn);pr=pr.a;
    fn=fileparts(fn);
  else
    fnp=[fn '/profile-ray.mat'];
    if isfile(fnp);
      pr=load(fnp);pr=pr.a;
    else
      disp(sprintf('ded_gc_f7_g_classify: No file "%s"',fnp));
      continue;
    end
  end
  j=j+1;
  p=load([fn '/param.mat']);p=p.p;
  
  if isfield(pr,'x')
    x=pr.x;
  else
    c=load([fn '/coord.mat']);
    x=c.c.Jx;
    clear('c');
  end
  dx=x(2)-x(1);
  if length(x)==length(pr.bm0)
    xmin=max(xmin,min(x(pr.bm0>0)));
  end
  fr=ded_gc_find_head(x,pr.bm0,xmin,display);
  xmax=p.L-xmin;

  du=max(0,pr.u1-pr.u2);
  db=max(0,pr.b1-pr.b2);
  dudz=du./pr.uw/sqrt(pi);
  dbdz=db./pr.bw/sqrt(pi);
  R=pr.uw./pr.bw;
  Re=p.Re*du.*pr.bh;
  Rib=p.g*db.*pr.uh./du.^2;
  J=dbdz./dudz.^2*p.g;
  
  if isfield(pr,'bDI');pr.SBI=pr.bDI;end;
  if isfield(pr,'SI');pr.SRI=pr.SI;end;
  if ~isfield(pr,'SBI');pr.SBI=NaN*pr.bm0;end;
  if ~isfield(pr,'SRI');pr.SRI=NaN*pr.bm0;end;

  f=find(x>xmin&x<=xmax);
  xf=x(f);
  if display
    figure;clf;
    subplot(5,1,1);plot(xf,pr.SBI(f),xf,pr.dbs(f));ylabel('SB');
    title(nm);
    subplot(5,1,2);plot(xf,pr.SRI(f),xf,pr.dus(f));ylabel('SR');
    subplot(5,1,3);plot(xf,pr.Fuw(f),xf,pr.Fww(f));ylabel('uw');
    subplot(5,1,4);plot(xf,pr.Fbu(f),xf,pr.Fbw(f));ylabel('buw');
    subplot(5,1,5);plot(xf,pr.Fbb(f));ylabel('bb');
    drawnow
  end
  for k=1:5
    switch(k)
      case(1)
        f=find(x>xmin&x<=fr.Xf);
      case(2)
        f=find(x>fr.Xt&x<=fr.Xf);
      case(3)
        f=find(x>xmin&x<=fr.Xt);
      case(4)
        f=find((fr.Xt-x)/fr.ht>=1 &  (fr.Xt-x)/fr.ht<=10);
      case(5)
        f=find(x>=25 & x<=30);
    end
    if fr.Xt>fr.Xf
      disp('de_gc_f7_g_classify: tail and front out of order');
    end
    b(k).x1(j)     = min(x(f)/fr.ht);
    b(k).x2(j)     = max(x(f)/fr.ht);
    b(k).Xf(j)     = fr.Xf;
    b(k).Xt(j)     = fr.Xt;
    b(k).Xh(j)     = fr.Xh;
    b(k).X0(j)     = fr.X0;
    b(k).ht(j)     = fr.ht;
    b(k).hh(j)     = fr.hh;
    b(k).h0(j)     = fr.h0;
    b(k).alpha(j)  = fr.alpha;
    b(k).alpha1(j) = (fr.h0-fr.ht)/(fr.Xt-fr.X0);
    b(k).SBI(j)  = max(0,sum(pr.SBI(f))*dx);
    b(k).SRI(j)  = max(0,sum(pr.SRI(f))*dx);
    b(k).SB(j)   = max(0,sum(pr.SBI(f)-pr.dbs(f))*dx);
    b(k).SR(j)   = max(0,sum(pr.SRI(f)-pr.dus(f))*dx);
    b(k).Fuu(j)  = sum(pr.Fuu(f))*dx;
    b(k).Fuw(j)  = sum(pr.Fuw(f))*dx;
    b(k).Fww(j)  = sum(pr.Fww(f))*dx;
    b(k).Fbu(j)  = sum(pr.Fbu(f))*dx;
    b(k).Fbw(j)  = sum(pr.Fbw(f))*dx;
    b(k).Fbb(j)  = sum(pr.Fbb(f))*dx;
    b(k).nRe(j)   = p.Re;  
    b(k).nPe(j)   = p.Re*p.Scb;
    b(k).typ{j}  = ded_gc_f7_g_classification(nm);
    b(k).typ2{j} = [ded_gc_f7_g_classification(nm) ' ' num2str(round(p.Re*p.Scb))];
    b(k).nm{j}   = nm;  
    b(k).R(j)    = mean(R(f));
    b(k).Re(j)   = mean(Re(f));
    b(k).Rib(j)  = mean(Rib(f));
    b(k).J(j)    = mean(J(f));
    b(k).sR(j)   = std(R(f));
    b(k).sRe(j)  = std(Re(f));
    b(k).sRib(j) = std(Rib(f));
    b(k).sJ(j)   = std(J(f));

    b(k).typ{j}  = b(k).typ{j} 
    c(k).Re(j)   = p.Re;
    c(k).uw(j)   = mean(pr.uw(f));
    c(k).uh(j)   = mean(pr.uh(f));
    c(k).u1(j)   = mean(pr.u1(f));
    c(k).u2(j)   = mean(pr.u2(f));
    c(k).bw(j)   = mean(pr.bw(f));
    c(k).bh(j)   = mean(pr.bh(f));
    c(k).b1(j)   = mean(pr.b1(f));
    c(k).b2(j)   = mean(pr.b2(f));
  end
end
return

