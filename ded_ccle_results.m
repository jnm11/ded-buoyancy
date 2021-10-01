nm={'gc/ccle/020',...
    'gc/ccle/021',...
    'gc/ccle/022',...
    'gc/ccle/023',...
    'gc/ccle/024',...
    'gc/ccle/025',...
    'gc/ccle/046'};   % v=@(t) polyval([-2.69e-4 0.4267],t);

%    'gc/ccle/011',...
%    'gc/ccle/012',...
%    'gc/ccle/013',...
%    'gc/ccle/014',...
%    'gc/ccle/016',...
%    'gc/ccle/017',...
%    'gc/ccle/018',...
%    'gc/ccle/020'};

a=ded_read_stats(nm);
p=ded_read_param(nm);


preprint([6 4],8);
colour_lines;
clf;
k=4;
n=1;
tol1=1e-4;
tol2=0.05;
display=0;
for j=1:length(a)
  t=a(j).t;
  X=a(j).X;
  v=(X(k+1:end)-X(1:end-k))./(t(k+1:end)-t(1:end-k));
  tv=(t(k+1:end)+t(1:end-k))/2;
  for jj=1:2
    if jj==2
      figure(1);
    else
      figure(j+1);clf;
    end
    subplot(2,1,1);
    hold('on');
    h(j)=plot(t,X);
    subplot(2,1,2);
    hold('on');
    [c t1 t2 z]=fit_poly_clip(tv,v,n,tol1,tol2,display);
    plot(tv,v,tv,z);
    axis('tight');
    xlabel('t');
    ylabel('v');
    l{j}=sprintf('Nx=%u',p(j).Nx);
    if jj==1
      subplot(2,1,1);
      title(sprintf('%s %5.2f %5.2f %8.2e %6.4f',nm{j},t1,t2,c(1),c(2)));
    end
  end
end
figure(1)
subplot(2,1,1);
legend(h,l,'location','NW');
axis('tight');
ylabel('X');

dd='~/collaborations/Claudia-GC/';
print('-depsc2','-loose',[dd 'ccle-v.eps']);
