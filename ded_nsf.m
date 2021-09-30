

nm='pm/025';
p=ded_read_param(nm);
load('~/data/dedalus/pm/025/a/a-r-00010.mat');


nx=length(a.x);

s=max(realmin,a.s./min(a.s));
b=max(realmin,a.b./max(a.b));


preprint([6 3],12);colour_lines;
k=50;
h=semilogy(a.r,s(k,:),a.r,b(k,:));
legend(h,{'sediment','bouyancy'},'location','se');
axis([0 inf 1e-4 1.1]);
title(sprintf('x=%6.2f',a.x(k)));
drawnow;
xlabel('r');
ylabel('s,b');
print('-depsc2','~/SB.eps');


a=ded_pm_r('pm/025',4);


