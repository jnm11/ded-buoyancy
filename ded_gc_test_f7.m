function p=ded_gc_test_f7(nm,trg)
if nargin<2
  trg=[0 inf];
end
c=ded_coord(nm);
a=ded_read_g(nm,'sev',[],trg);
d=ded_read_g(nm,'avar',[],trg);
b=ded_read_javrg(nm,'a',trg);
bb=ded_read_javrg(nm,'a',trg,'separate');
if isempty(b)
  b=ded_read_javrg(nm,'ay',trg);
  bb=ded_read_javrg(nm,'ay',trg,'separate');
end
p=ded_read_param(nm);
s=ded_read_stats(nm);
w=ded_hdf(nm,'fx','force');
rfu=ded_hdf(nm,'rfu','force');
a.z=c.Az;
a.x=c.Ax;
b.x=c.Jx;
b.z=c.Jz;
dimx=c.dim;

dc;
% sdl is velocity sample region
% sbl is density sample  region
% wdl is density forced  region

f=intersect({'uul','uur','vvl','vvr','wwl','wwr'},fieldnames(p));
for j=1:length(f)
  if p.(f{j})
    disp(sprintf('%s True',f{j}));
  else
    disp(sprintf('%s False',f{j}));
  end
end
if isfield(w,'sdl');nsdl=length(w.sdl);xsdl=c.Ax(1:nsdl);Isdl=sum(w.sdl*c.dAx);ssdl=sprintf('d sample %4.2f',Isdl);end;
if isfield(w,'wdl');nwdl=length(w.wdl);xwdl=c.Ax(1:nwdl);Iwdl=sum(w.wdl*c.dAx);swdl=sprintf('d source %4.2f',Iwdl);end;
if isfield(w,'wul');nwul=length(w.wul);xwul=c.Ax(1:nwul);Iwul=sum(w.wul*c.dAx);swul=sprintf('u source %4.2f',Iwul);end;
if isfield(w,'sbl');nsbl=length(w.sbl);xsbl=c.Ax(1:nsbl);Isbl=sum(w.sbl*c.dAx);ssbl=sprintf('b sample %4.2f',Isbl);end;
if isfield(w,'wbl');nwbl=length(w.wbl);xwbl=c.Ax(1:nwbl);Iwbl=max(w.wbl);      swbl=sprintf('b source %4.2f',Iwbl);end;

if isfield(w,'wdr');nwdr=length(w.wdr);xwdr=c.Ax(end+1-nwdr:end);Iwdr=sum(w.wdr*c.dAx);swdr=sprintf('d source %4.2f',Iwdr);end;
if isfield(w,'wur');nwur=length(w.wur);xwur=c.Ax(end+1-nwur:end);Iwur=sum(w.wur*c.dAx);swur=sprintf('u source %4.2f',Iwur);end;
if isfield(w,'wbr');nwbr=length(w.wbr);xwbr=c.Ax(end+1-nwbr:end);Iwbr=max(w.wbr);swbr=sprintf('b source %4.2f',Iwbr);end;

tol=1e-6;

if isfield(w,'wdl') ; if isfield(w,'wdls') xwdl=c.Ax(w.wdls+1) ; else f=(w.wdl>tol*max(w.wdl));xwdl=c.Ax(f);w.wdl=w.wdl(f);end; end
if isfield(w,'sdl') ; if isfield(w,'sdls') xsdl=c.Ax(w.sdls+1) ; else f=(w.sdl>tol*max(w.sdl));xsdl=c.Ax(f);w.sdl=w.sdl(f);end; end
if isfield(w,'wdr') ; if isfield(w,'wdrs') xwdr=c.Ax(w.wdrs+1) ; else f=(w.wdr>tol*max(w.wdr));xwdr=c.Ax(f);w.wdr=w.wdr(f);end; end
if isfield(w,'sdr') ; if isfield(w,'sdrs') xsdr=c.Ax(w.sdrs+1) ; else f=(w.sdr>tol*max(w.sdr));xsdr=c.Ax(f);w.sdr=w.sdr(f);end; end
                                                                                                                       
if isfield(w,'wbl') ; if isfield(w,'wbls') xwbl=c.Ax(w.wbls+1) ; else f=(w.wbl>tol*max(w.wbl));xwbl=c.Ax(f);w.wbl=w.wbl(f);end; end
if isfield(w,'sbl') ; if isfield(w,'sbls') xsbl=c.Ax(w.sbls+1) ; else f=(w.sbl>tol*max(w.sbl));xsbl=c.Ax(f);w.sbl=w.sbl(f);end; end
if isfield(w,'wbr') ; if isfield(w,'wbrs') xwbr=c.Ax(w.wbrs+1) ; else f=(w.wbr>tol*max(w.wbr));xwbr=c.Ax(f);w.wbr=w.wbr(f);end; end
if isfield(w,'sbr') ; if isfield(w,'sbrs') xsbr=c.Ax(w.sbrs+1) ; else f=(w.sbr>tol*max(w.sbr));xsbr=c.Ax(f);w.sbr=w.sbr(f);end; end
                                                                                                                       
if isfield(w,'wul') ; if isfield(w,'wuls') xwul=c.Ax(w.wuls+1) ; else f=(w.wul>tol*max(w.wul));xwul=c.Ax(f);w.wul=w.wul(f);end; end
if isfield(w,'sul') ; if isfield(w,'suls') xsul=c.Ax(w.suls+1) ; else f=(w.sul>tol*max(w.sul));xsul=c.Ax(f);w.sul=w.sul(f);end; end
if isfield(w,'wur') ; if isfield(w,'wurs') xwur=c.Ax(w.wurs+1) ; else f=(w.wur>tol*max(w.wur));xwur=c.Ax(f);w.wur=w.wur(f);end; end
if isfield(w,'sur') ; if isfield(w,'surs') xsur=c.Ax(w.surs+1) ; else f=(w.sur>tol*max(w.sur));xsur=c.Ax(f);w.sur=w.sur(f);end; end
dc






nsfd=+isfield(w,'wdl')+isfield(w,'sdl')+isfield(w,'wdr')+isfield(w,'sdr');
nsfb=+isfield(w,'wbl')+isfield(w,'sbl')+isfield(w,'wbr')+isfield(w,'sbr');
nsfu=+isfield(w,'wul')+isfield(w,'sul')+isfield(w,'wur')+isfield(w,'sur');



figure;j=0;
if isfield(w,'wdl') j=j+1;subplot(nsfd,1,j);plot(xwdl,w.wdl);ylabel('wdl');title(sprintf('l d source %4.2f',sum(w.wdl*c.dAx)));end;
if isfield(w,'sdl') j=j+1;subplot(nsfd,1,j);plot(xsdl,w.sdl);ylabel('sdl');title(sprintf('l d sample %4.2f',sum(w.sdl*c.dAx)));end;
if isfield(w,'wdr') j=j+1;subplot(nsfd,1,j);plot(xwdr,w.wdr);ylabel('wdr');title(sprintf('r d source %4.2f',sum(w.wdr*c.dAx)));end;
if isfield(w,'sdr') j=j+1;subplot(nsfd,1,j);plot(xsdr,w.sdr);ylabel('sdr');title(sprintf('r d sample %4.2f',sum(w.sdr*c.dAx)));end;
xlabel('x');

figure;j=0;
if isfield(w,'wbl') j=j+1;subplot(nsfb,1,j);plot(xwbl,w.wbl);ylabel('wbl');title(sprintf('l b source %4.2f',sum(w.wbl*c.dAx)));end;
if isfield(w,'sbl') j=j+1;subplot(nsfb,1,j);plot(xsbl,w.sbl);ylabel('sbl');title(sprintf('l b sample %4.2f',sum(w.sbl*c.dAx)));end;
if isfield(w,'wbr') j=j+1;subplot(nsfb,1,j);plot(xwbr,w.wbr);ylabel('wbr');title(sprintf('r b source %4.2f',sum(w.wbr*c.dAx)));end;
if isfield(w,'sbr') j=j+1;subplot(nsfb,1,j);plot(xsbr,w.sbr);ylabel('sbr');title(sprintf('r b sample %4.2f',sum(w.sbr*c.dAx)));end;
xlabel('x');

figure;j=0;
if isfield(w,'wul') j=j+1;subplot(nsfu,1,j);plot(xwul,w.wul);ylabel('wul');title(sprintf('l u source %4.2f',sum(w.wul*c.dAx)));end;
if isfield(w,'sul') j=j+1;subplot(nsfu,1,j);plot(xsul,w.sul);ylabel('sul');title(sprintf('l u sample %4.2f',sum(w.sul*c.dAx)));end;
if isfield(w,'wur') j=j+1;subplot(nsfu,1,j);plot(xwur,w.wur);ylabel('wur');title(sprintf('r u source %4.2f',sum(w.wur*c.dAx)));end;
if isfield(w,'sur') j=j+1;subplot(nsfu,1,j);plot(xsur,w.sur);ylabel('sur');title(sprintf('r u sample %4.2f',sum(w.sur*c.dAx)));end;
xlabel('x');



n=2*p.dd+isfield(a,'db')+isfield(d,'db');
k=0;
Ibw  = ichebintw(c.NAz);
Iaw  = ichebintw(c.NJz);
dim=(p.Nx>1)+(p.Ny>1)+(p.Nz>1);

Nxx=round(c.Nx*p.AA/p.AAJ);

if isfield(b,'u')
  sz=size(b.u);
  sz(end)=Nxx;
  if 0 % Check parity resampling
    b=ded_read_javrg(nm,'a',trg);
    d1=b.d(1,:);d2=fftresample(b.d,sz,2,[0  1]);d2=d2(2,:);
    u1=b.u(1,:);u2=fftresample(b.u,sz,2,[0 -1]);u2=u2(2,:);
    subplot(4,1,1);plot(c.Jx,d1,'s-',c.Ax,d2,'s-');axis([0 0.5 -inf inf])
    subplot(4,1,2);plot(c.Jx,d1,'s-',c.Ax,d2,'s-');axis([p.L-0.5 p.L -inf inf])
    subplot(4,1,3);plot(c.Jx,u1,'s-',c.Ax,u2,'s-');axis([0 0.5 -inf inf])
    subplot(4,1,4);plot(c.Jx,u1,'s-',c.Ax,u2,'s-');axis([p.L-0.5 p.L -inf inf])
  end
end


ff={'u','v','w','d','b','s','p','uu'};
P= [-1   1   1   1   1   1   1    1];
if Nxx~=c.Nx
  for k=1:length(ff)
    f=ff{k};
    if isfield(b,f);
      PP=zeros(1,1,dim);
      PP(dimx)=P(k);
      b.(f)=fftresample(b.(f),sz,PP);
      for j=1:length(bb)
        bb(j).(f)=fftresample(bb(j).(f),sz,PP);
      end
    end
  end
end

sz=ones(1,dim);
sz(end)=Nxx;

sdl=NaN;
sbl=NaN;
wdl=NaN;

if isfield(w,'sdl');sdl=zeros(sz);sdl(1:length(w.sdl))=w.sdl*c.dAx;end;
if isfield(w,'sbl');sbl=zeros(sz);sbl(1:length(w.sbl))=w.sbl*c.dAx;end;
if isfield(w,'wdl');wdl=zeros(sz);wdl(1:length(w.wdl))=w.wdl*c.dAx;end;

if isfield(b,'b')
  sb    = sum(b.b.*sbl,c.dimx);
  wpl   = sum(b.p.*wdl,c.dimx);
  wpr   = sum(b.p.*wdl,c.dimx);
  uz    = sum(b.w.*sdl,c.dimx);
  uxs   = sum(b.u.*sdl,c.dimx);
  uxw   = sum(b.u.*wdl,c.dimx);

  dw    = ichebdifff(uz,1,p.H);
  buxst = Iaw*squeeze(          sum(cat(3,bb.u).*sdl,c.dimx));
  dblt  = Iaw*squeeze(          sum(cat(3,bb.p).*wdl,c.dimx));
  dbrt  = Iaw*squeeze(          sum(cat(3,bb.p).*wdl,c.dimx));
  buzt  = sqrt(Iaw*squeeze(           sum(cat(3,bb.w).*sdl,c.dimx)).^2);
  bdwt  = sqrt(Iaw*squeeze(ichebdifff(sum(cat(3,bb.w).*sdl,c.dimx),1,p.H)).^2);
  
  
  Nl=length(w.wdl);rgl=1:Nl;         xl=c.Ax(rgl);
  Nr=length(w.wdr);rgr=Nxx-Nr+1:Nxx; xr=c.Ax(rgr);
  
  dp = pr_diff(b.p,c.dAx, dimx,[],[],1);
  du = pr_diff(b.u,c.dAx, dimx,[],[],-1);
  dpxI  = Iaw*dp;
  duuxI = Iaw*(b.u.*du);
  pI  = Iaw*b.p;
  uuI = Iaw*b.uu;
  bgI = p.g*Iaw*b.b;
  bI  = pr_int(bgI,c.dAx, dimx,[],1);
  
  figure;clf;
  subplot(4,1,1);plot(xl,uuI(rgl)/2,'s-',xl,-pI(rgl)+pI(1),  's-');axis([0 xl(end) -inf inf]);
  subplot(4,1,2);plot(xr,uuI(rgr)/2,'s-',xr,-pI(rgr)+pI(end),'s-');axis([xr(1) p.L -inf inf]);
  subplot(4,1,3);plot(xl,uuI(rgl)/2+pI(rgl)-pI(1)  ,'s-');axis([0 xl(end) -inf inf]);
  subplot(4,1,4);plot(xr,uuI(rgr)/2+pI(rgr)-pI(end),'s-');axis([xr(1) p.L -inf inf]);
  gu_title('Integrated uu p balance'); 
  
  figure;clf;
  subplot(4,1,1);plot(xl,duuxI(rgl),'s-',xl,-dpxI(rgl),'s-');axis([0 xl(end) -inf inf]);
  subplot(4,1,2);plot(xr,duuxI(rgr),'s-',xr,-dpxI(rgr),'s-');axis([xr(1) p.L -inf inf]);
  subplot(4,1,3);plot(xl,duuxI(rgl)+dpxI(rgl),         's-');axis([0 xl(end) -inf inf]);
  subplot(4,1,4);plot(xr,duuxI(rgr)+dpxI(rgr),         's-');axis([xr(1) p.L -inf inf]);
  gu_title('duux dpdx balance'); 
  
end

%subplot(1,n,1);plot(uz,b.z);xlabel('w');axis('tight');
  %subplot(1,n,2);plot(uxs,b.z,uxw,b.z);xlabel('u');axis('tight');
  %title(sprintf('%5.3f %5.3f',-Ibw*uxs,-Ibw*uxw));
  %subplot(1,n,3);plot(dw,b.z);xlabel('dw');axis('tight');
  %k=k+3;

addt=repmat(NaN,length(a.t),1);
dwdzt=addt;
adpt=NaN;
dpt=NaN;



n=isfield(a,'db')+isfield(a,'dd')+isfield(d,'db')+isfield(d,'ww')+isfield(d,'dwdz');
figure;clf;
k=0;

if isfield(a,'db');  k=k+1;subplot(1,n,k);plot(a.db,  a.z);title('sev'); xlabel('db');  axis([0 1 0 p.H]);end
if isfield(d,'db');  k=k+1;subplot(1,n,k);plot(d.db,  a.z);title('avar');xlabel('db');  axis([0 1 0 p.H]);end
if isfield(a,'dd');  k=k+1;subplot(1,n,k);plot(a.dd,  a.z);title('sev'); xlabel('dd');  axis([-inf inf 0 p.H]);axis('tight');end
if isfield(d,'dwdz');k=k+1;subplot(1,n,k);plot(d.dwdz,a.z);title('avar');xlabel('dwdz');axis([-inf inf 0 p.H]);axis('tight');end
if isfield(d,'ww');  k=k+1;subplot(1,n,k);plot(d.ww,  a.z);title('avar');xlabel('ww');  axis([-inf inf 0 p.H]);axis('tight');end

figure;clf;
k=0;

if isfield(a,'db');  k=k+1;subplot(n,1,k);plot(a.t,Ibw*a.db  );ylabel('sev  db');  axis([-inf inf 0 1 ]);end
if isfield(d,'db');  k=k+1;subplot(n,1,k);plot(d.t,Ibw*d.db  );ylabel('avar db');  axis([-inf inf 0 1 ]);end
if isfield(a,'dd');  k=k+1;subplot(n,1,k);plot(a.t,sqrt(Ibw*(a.dd-p.V).^2/p.H)  );ylabel('rms dd');  axis([-inf inf -inf inf]);axis('tight');end
if isfield(d,'dwdz');k=k+1;subplot(n,1,k);plot(d.t,sqrt(Ibw*d.dwdz.^2/p.H));ylabel('rms dwdz');axis([-inf inf -inf inf]);axis('tight');end
if isfield(d,'ww');  k=k+1;subplot(n,1,k);plot(d.t,Ibw*d.ww  );ylabel('avar ww');  axis([-inf inf -inf inf]);axis('tight');end
xlabel('t');

% Plot w dw u du near the end wall
figure;clf;
N=3;
ww=b.w(:,1:N);
dw=ichebdifff(ww,c.dimz,p.H,1);
uu=b.u(:,1:N);
du=ichebdifff(uu,c.dimz,p.H,1);
lb={'w','dw','u','du'};
for k=1:4
  subplot(1,4,k);
  switch(k)
    case(1)
      f=ww;
    case(2)
      f=dw;
    case(3)
      f=uu;
    case(4)
      f=du;
  end
  maxf=cumsum([0 max(f)-min(f)]);
  maxf(end)=[];
  xw=maxf/2;
  h=plot(f+xw,c.Jz);
  set(gca,'xtick',xw,'xticklabels',cellsprintf('%d',1:N));
  for j=1:N
    set(h(j),'color',[1 1 1]*(j-1)/(N-1));
  end
  line([xw;xw],[0 p.H],'color',[0.7 0.7 0.7]);
  xlabel('w');
  axis('tight');
  title(sprintf('%7.3f %7.3f',min(f(:)),max(f(:))));
  xlabel(lb{k});
  axis('tight');
end
gu_title('Left profiles'); 

if isfield(s,'dd')
  figure;
  f= find(s.t>=trg(1) & s.t<=trg(2));
  t=s.t(f);
  for k=1:3
    switch(k)
      case 1
        g=s.dd(f,:);
        l='dd';
      case 2
        if ~isfield(s,'dw');continue;end
        g=s.dw(f,:);
        l='dw';
      case 3
        if isfield(p,'PIDX')
          g=s.X(f,:)-p.PIDX;
        else
          g=s.X(f,:)-mean(s.X(f,:));
        end
        l='X';
    end
    subplot(3,1,k);
    plot(t,g);
    line(t,t*0,'color',[0.7 0.7 0.7]);
    axis('tight');
    ylabel(l);
    set(gca,'ytick',[])
    title(sprintf('%7.0e %7.0e',min(g(:)),max(g(:))));
    if k==3;
      xlabel('t');
    else
      set(gca,'xtick',[]);
    end
  end
  gu_title('Time series'); 
end

if isfield(w,'wudr')
  if isfield(w,'wdrs')
    fx=w.wdrs+1;
    wudr=w.wudr;
  else
    x1=min(c.Ax(abs(w.wudr)>0.01*max(abs(w.wudr))));
    fx=find(c.Ax>=x1);
    wudr=w.wudr(fx);
  end
  figure;
  x=c.Ax(fx);
  
  figure;
  title('wudr')
  if isfield(b,'d')
    subplot(4,1,1);
    dd=sum(b.u(:,fx).*b.d(:,fx).*wudr',2)/sum(wudr.^2);  
    mean(a.dd,2);
    plot(x,b.u(:,fx).*b.d(:,fx)-dd*wudr','b')
    axis('tight');
    ylabel('ududx-dd*f');
  end
  subplot(4,1,2);
  plot([x;p.L],[wudr;0],'bs-')
  axis('tight');
  ylabel('f');
  subplot(4,1,3);cla;
  plot(x,b.p(:,fx(end))-b.p(end,fx),'g');
  axis('tight');
  ylabel('p');
  subplot(4,1,4);cla;
  plot(x,b.u(:,fx));
  axis('tight');
  ylabel('u');
end

if isfield(w,'wudl')
  if isfield(w,'wdls')
    fx=w.wdls+1;
    wudl=w.wudl;
  else
    x1=max(c.Ax(abs(w.wudl)>0.01*max(abs(w.wudl))));
    fx=find(c.Ax<=x1);
    wudl=w.wudl(fx);
  end
  figure;
  x=c.Ax(fx);
  figure;
  title('wudl')
  if isfield(b,'d')
    dd=sum(b.u(:,fx).*b.d(:,fx).*wudl',2)/sum(wudl.^2);  
    subplot(4,1,1);
    mean(a.dd,2);
    plot(x,b.u(:,fx).*b.d(:,fx)-dd*wudl','b')
    axis('tight');
    ylabel('ududx-dd*f');
  end
  subplot(4,1,2);
  plot([x;p.L],[wudl;0],'bs-')
  axis('tight');
  ylabel('f');
  subplot(4,1,3);cla;
  plot(x,b.p(:,fx(end))-b.p(end,fx),'g');
  axis('tight');
  ylabel('p');
  subplot(4,1,4);cla;
  plot(x,b.u(:,fx));
  axis('tight');
  ylabel('u');
end


