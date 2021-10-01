function ded_gc_X_U(n,Xtyp,Ytyp)

if nargin==0
  figure(1);clf;ded_gc_X_U('6*','U','X');
  figure(2);clf;ded_gc_X_U('5*','Re','U');
  figure(3);clf;ded_gc_X_U('5*','Re','X');
  figure(4);clf;ded_gc_X_U('7*','Re','U');
  figure(5);clf;ded_gc_X_U('7*','Re','X');
  return;
end

if nargin<1
  n='';
end
if nargin<2
  Xtyp='';
end
if nargin<3
  Ytyp='';
end

if isempty(Xtyp)
  Xtyp='U';
end
if isempty(Ytyp)
  Ytyp='X';
end
if isempty(n)
  n='*';
end

a=ded_all_stats('gc',n);


HX=unique([[a.H]' [a.PIDX]'],'rows');

n=size(HX,1);
Hu=HX(:,1);
Xu=HX(:,2);


for j=1:n
  figure(j);
  clf;hold('on');
  f=find(Hu(j)==[a.H] & Xu(j)==[a.PIDX]);
  b=a(f);
  
  switch(Xtyp)
    case 'U'
      x=abs([b.U]);
    case 'Re'
      x=[b.Re];
  end

  [x f]=sort(x);

  switch(Ytyp)
    case 'X'
      y=[b(f).fXm]; 
    case 'U'
      y=[b(f).U]; 
  end  
  GS=cellsprintf('%4.1f',[b(f).Nx].*[b(f).Ny].*[b(f).Nz]/1e6);
  
  T=[b.T];
  pp.lt='-';
  [h lh uc]=groupplot(x,y,GS,pp);
  if length(lh)>1
    legend(lh,uc,'location','best');
  end
  
  %  g1= isfinite(T);
  %  g2=~isfinite(T);
  %  plot(U(g1), Xm(g1),'markerfacecolor',[0 0 0],'markeredgecolor',[0 0 1],'linestyle','none');
  %  plot(U(g2), Xm(g2),'markerfacecolor','none','markeredgecolor',[1 0 0],'linestyle','none');
  %  legend(h,{'moments','threshold'});
  xlabel(Xtyp);
  ylabel(Ytyp);
  title(sprintf('H=%4.1f, PIDX %4.1f',Hu(j),Xu(j)));
  axis('tight');
  aa=axis;
  mU=aa(1);xU=aa(2);dU=0.05*(xU-mU);
  mX=aa(3);xX=aa(4);dX=0.05*(xX-mX);
  rg =[mU-dU xU+dU mX-dX xX+dX]; 
  axis(rg);
end
