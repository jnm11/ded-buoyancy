function a=ded_gc_f7_g_RJ(fns,xmin,loc,display)
j=0;
a=struct;
for jj=1:length(fns)
  fn=[ded_dedalus_data_dir '/results/' fns{jj}];
  if isfile(fn)
    disp(sprintf('ded_gc_f7_g_RJ: loading %s',fn));
    pr=load(fn);pr=pr.a;
    fn=fileparts(fn);
  else 
    fnp=[fn '/profile-ray.mat'];
    if isfile(fnp);
      pr=load(fnp);pr=pr.a;
    else
      disp(sprintf('ded_gc_f7_g_RJ: No file "%s"',fnp));
      continue;
    end
  end
  if ~isfield(pr,'bI');continue;end
  j=j+1;
    
  nm=pr.nm;
  p=load([fn '/param.mat']);p=p.p;
  %c=load([fn '/coord.mat']);c=c.c;
  fr=ded_gc_find_head(pr.x,pr.bI,xmin,display);

  x=pr.x(:)';
  switch(loc)
    case('front')
      X0=fr.Xf;
    case('tail')
      X0=fr.Xt;
  end
  
  f=find(x<=X0 & x>xmin);
  x=(X0-x)/fr.ht;

  R=pr.uw./pr.bw;
  du=max(0,pr.u1-pr.u2);
  db=max(0,pr.b1-pr.b2);
  dudz=du./pr.uw/sqrt(pi);
  dbdz=db./pr.bw/sqrt(pi);
  Re=p.Re*du.*pr.bh;
  Rib=p.g*db.*pr.uh./du.^2;
  J=dbdz./dudz.^2*p.g;
  Fr=sqrt(max(0,(pr.uuI-p.H)./(p.g*pr.bz)));
  aa.Fbw=pr.Fbw(f);
  aa.Tbw=pr.Tbw(f);
  aa.x=x(f);
  aa.x=x(f);
  aa.R=R(f);
  aa.Re=Re(f);
  aa.Rib=Rib(f);
  aa.J=J(f);
  aa.Fr=Fr(f);
  aa.b=pr.bI(f);
  aa.nRe=p.Re; 
  aa.nPe=round(p.Scb*p.Re);
  aa.g = p.g;
  aa.nm=nm;
  aa.typ= ded_gc_f7_g_classification(nm);
  aa.loc=loc;
  aa.Xf=(X0-fr.Xf)/fr.ht;
  aa.Xh=(X0-fr.Xh)/fr.ht;
  aa.Xt=(X0-fr.Xt)/fr.ht;
  a(j)=aa;
end


return;

figure(1);clf;hold('on');
for j=1:length(a);plot(a(j).R,a(j).J);end
set(gca,'xscale','log','yscale','log');
xlabel('R');ylabel('J');

figure(2);clf;hold('on');
for j=1:length(a);plot(a(j).R,a(j).Re);end
set(gca,'xscale','log','yscale','log');
xlabel('R');ylabel('Re');

figure(3);clf;hold('on');
for j=1:length(a);plot(a(j).R,a(j).Rib);end
set(gca,'xscale','log','yscale','log');
xlabel('R');ylabel('Rib');

figure(4);clf;hold('on');
for j=1:length(a);plot(a(j).x,a(j).J);end
set(gca,'xscale','linear','yscale','log');
xlabel('x');ylabel('J');


X=[0 25 50 75 100];
col=[[1 0 0];[0 0 1];[1 1 0];[1 0 1];[0 0 1];[0 0 0]];
X=linspace(0,100,200);
col=jet(length(X));
mk='s<>^vpdo*';
k=50;plot(a(k).x,a(k).R)
figure(1);gfd_plot_along_x(a,X,'R','J',col,mk);

figure(2);gfd_plot_along_x(a,X,'R','Re',col,mk);

figure(3);gfd_plot_along_x(a,X,'R','Rib',col,mk);

figure(4);gfd_plot_along_x(a,X,'J','Fbw',col,mk);set(gca,'yscale','linear');


figure(5);gfd_plot_along_x(a,X,'R','Fbw',col,mk);set(gca,'yscale','linear');



figure(2);gfd_plot_along_x(a,X,'Re','Fr',col,mk);set(gca,'yscale','linear');



f=find([a.nPe]==8000);% & [a.nRe]<=1000);
                      %gfd_plot_along_x(a(f),X,'R','J',col,mk);

clf;
for k=1:length(f)
  j=f(k);
  subplot(3,1,1);  plot(a(j).x,a(j).J);hold('on');ylabel('J');
  subplot(3,1,2);  plot(a(j).x,a(j).R);hold('on');ylabel('R');
  subplot(3,1,3);  plot(a(j).R,a(j).J);hold('on'); xlabel('R');ylabel('J');
end
