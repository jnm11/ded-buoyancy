function a=ded_bg_plot(nm)
%~/python/beta/ded_gc.py --Nx 1024 -g 0 -L 10 --xa 1 --xb 2 --dtu 0.1 --Re 1000 bg/01 
%nm='bg/01';
p=ded_read_param(nm);
a=ded_read_g(nm,'u');

minu=min(a.u(:))-0.01;
maxu=max(a.u(:))+0.01;

a.dx=(a.x(2)-a.x(1));
a.tv=sum(abs(diff([a.u;a.u(1,:)],1,1)),1)*a.dx;
a.iv=sum(a.u,1)*a.dx;
figure(1);clf;
disp(max(diff(a.tv)));

subplot(3,1,1);plot(a.t,a.tv);
ylabel('TV');
subplot(3,1,2);plot(filter_midpoint(a.t),diff(a.tv),filter_midpoint(a.t),diff(a.iv));
subplot(3,1,3);
for j=1:length(a.t)
  plot(a.x,a.u(:,j));
  axis([0 p.L minu maxu]);
  drawnow;
end


