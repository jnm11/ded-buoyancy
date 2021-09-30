function ded_test_inertial(num)
nm=['gc/inertial/' num];
fns=ded_get_fn(nm,'s');
fnu=ded_get_fn(nm,'u');
fnw=ded_get_fn(nm,'w');
fnp=ded_get_fn(nm,'p');
fu=ded_hdf(nm,'fu','force');if ~isempty(fu);fu=fu.fu;end;
fw=ded_hdf(nm,'fw','force');if ~isempty(fw);fw=fw.fw;end;
c=ded_coord(nm);
if ~isempty(fu) & ~isempty(fw); 
  [fp fq fS fW fD]=pressure_fft(fu',fw',c.dJx,c.dJz);
  %plot(fq(:)-(fS(:)-fW(:)),'.');
end

p=ded_read_param(nm);

eu=[];cu=[];
ew=[];cw=[];
ep=[];cp=[];
eq=[];cq=[];
t=[];

x=c.Jx;
z=c.Jz;

clf;
if isempty(fns)
  ah=jsubplot([1 2],[0.05 0.05],[0.05 0.05],[0.01 0.05]);
  hp=ah(1);
  hq=ah(2);
else
  ah=jsubplot([2 2],[0.05 0.05],[0.05 0.05],[0.01 0.05]);
  hs=ah(1,1);
  hh=ah(2,1);
  hp=ah(1,2);
  hq=ah(2,2);
  
end
th=gu_title('');

if ~isempty(fns)
  axes(hs);cla;
  ih1=imagesc(x,z,x*z');
  title('s');
  hold('on');
  ch1=[];
  colorbar('EastOutside')
end

axes(hp);cla;
ih2=imagesc(x,z,x*z');
title('p');
hold('on');
ch2=[];
colorbar('EastOutside')


axes(hq);cla;
ih3=imagesc(x,z,x*z');
title('q');
hold('on');
ch3=[];
colorbar('EastOutside')

set(ah,'ytick',[]);
set(ah(:,1),'xtick',[]);

for j=1:length(fns)
  as=ded_read_hdf(fns{j});
  s.mean(j) = mean(as.s(:));
  s.std(j)  = std(as.s(:));
  s.min(j)  = min(as.s(:));
  s.max(j)  = max(as.s(:));
  s.t(j)    = double(as.t);
end
shist=linspace(min(s.min),max(s.max),100);

n=max([length(fns),length(fnu),length(fnw),length(fnp)]);
for j=1:n
  
  if j<=length(fns);as=ded_read_hdf(fns{j});else;as=[];end
  if j<=length(fnu);au=ded_read_hdf(fnu{j});else;au=[];end
  if j<=length(fnw);aw=ded_read_hdf(fnw{j});else;aw=[];end
  if j<=length(fnp);ap=ded_read_hdf(fnp{j});else;ap=[];end
  
  P=ap.p-min(ap.p(:));
  
  mu=sqrt(mean(mean(au.u.^2+aw.w.^2)));
  [f g w d sr rr q]=helmholtz_fft(au.u',aw.w',c.dJx,c.dJz);
  [pp pq ps pr pd]=pressure_fft(au.u',aw.w',c.dJx,c.dJz);
  %pe(:,j)=[max(abs(pq(:)-q(:))) max(abs(ps(:)-sr(:))) max(abs(pr(:)-rr(:))) max(abs(pd(:)-d(:)))];
    
  if ~isempty(as);cc=corrcoef(q(:),as.s(:));else;cc=NaN;end;

  
  cnts=contourc(au.x,au.z,f',20);
  set(th,'string',sprintf('t=%5.2f, c=%6.3f, u=%6.3f',au.t,cc,mu));

  if ~isempty(as);
    axes(hs);
    set(ih1,'cdata',as.s');
    delete(ch1);
    ch1=plot_contours(cnts);
    set(gca,'clim',[shist(1) shist(end)]);
    set(ch1,'color',[1 1 1]);

    axes(hh);
    hist(as.s(:),shist);
  
  end
  
% $$$   clf;
% $$$   quiver(c.Jx,c.Jz,au.u,aw.w)
% $$$   hold('on');
% $$$   plot_contours(cnts);
% $$$   
  
  axes(hp);
  set(ih2,'cdata',P);
  delete(ch2);
  ch2=plot_contours(cnts);
  set(gca,'clim',[min(P(:)) eps+max(P(:))]);
  set(ch2,'color',[1 1 1]);
  
  axes(hq);
  set(ih3,'cdata',q);
  delete(ch3);
  ch3=plot_contours(cnts);
  set(gca,'clim',[min(q(:)) eps+max(q(:))]);
  set(ch3,'color',[1 1 1]);

  %drawnow;
  if ~isempty(fu);[eu(end+1) cu(end+1)]=errcorr(au.u,fu); end
  if ~isempty(fw);[ew(end+1) cw(end+1)]=errcorr(aw.w,fw); end
  if ~isempty(fp);[ep(end+1) cp(end+1)]=errcorr(ap.p,fp); end
  if ~isempty(fq);[eq(end+1) cq(end+1)]=errcorr(   q,fq); end
  t(end+1)=au.t;
end

figure(2);clf;
ah=jsubplot([1 4],[0.2 0.05],[0.05 0.05],[0.01 0.05]);
if ~isempty(eu);axes(ah(1));plot(t,log10(eu));ylabel('u');gu_setalim;end
if ~isempty(ew);axes(ah(2));plot(t,log10(ew));ylabel('w');gu_setalim;end
if ~isempty(ep);axes(ah(3));plot(t,log10(ep));ylabel('p');gu_setalim;end
if ~isempty(eq);axes(ah(4));plot(t,log10(eq));ylabel('q');gu_setalim;end
set(ah(1:3),'xtick',[]);

figure(3);clf;
ah=jsubplot([1 4],[0.2 0.05],[0.05 0.05],[0.01 0.05]);
if ~isempty(cu);axes(ah(1));plot(t,cu);ylabel('u');gu_setalim;end
if ~isempty(cw);axes(ah(2));plot(t,cw);ylabel('w');gu_setalim;end
if ~isempty(cp);axes(ah(3));plot(t,cp);ylabel('p');gu_setalim;end
if ~isempty(cq);axes(ah(4));plot(t,cq);ylabel('q');gu_setalim;end
set(ah(1:3),'xtick',[]);



function [e c cx cy]=errcorr(x,y)
sxx=mean(x(:).^2);
syy=mean(y(:).^2);
sxy=mean(x(:).*y(:));
e=sqrt(mean((x(:)-y(:)).^2));
c=sxy/sqrt(sxx*syy);
cx=sxy/sxx;
cy=sxy/syy;







return;
if isfield(p,'fuspecial')
  U=p.U;
  H=p.H;
  L=p.L;
  %  switch(p.fuspecial)
  %  case 'sinsin'
      kx=2*pi/L;
      kz=2*pi/H;
      kk=sqrt((kx^2+kz^2)/2);
      [xx zz]=ndgrid(kx*c.Jx,kz*c.Jz);
      psi = U*sin(xx).*sin(zz)/kk;
      fu = +U*sin(xx).*cos(zz)*kz/kk;
      fw = -U*cos(xx).*sin(zz)*kx/kk;
      fp = U^2*(kz^2*cos(2*xx) + kx^2*cos(2*zz))/kk/4;
      fq = U^2*(     cos(2*xx) +      cos(2*zz))*kx^2*kz^2/kk^2;
      fww= 2*U*sin(xx).*sin(zz)*kk;
      fss= U^2/kk^2*(2*kx^2*kz^2*cos(xx).^2.*cos(zz).^2+sin(xx).^2.*sin(zz).^2*(kx^2-kz^2)^2/2);
      [f g ww d ss rr q]=helmholtz_fft(fu,fw,c.dJx,c.dJz);
      [fe fc fcx fcy]=errcorr(f,psi);
      [qe qc qcx qcy]=errcorr(q,fq);
      [we wc wcx wcy]=errcorr(ww,fww);
      [re rc rcx rcy]=errcorr(rr,fww.^2/2);
      [se sc scx scy]=errcorr(ss,fss);
      disp([fe qe we re se]);
      fu=fu';
      fw=fw';
      fp=fp';
      fq=fq';
      %  end
end
