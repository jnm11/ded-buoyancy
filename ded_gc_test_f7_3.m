function ded_gc_test_f7_3(nm)
p=ded_read_param(nm);
c=ded_coord(nm);
w=ichebintw(c.NAz);

a=ded_read_g(nm,'sev');
b=ded_read_g(nm,'avar');
mindd=min(a.dd);
maxdd=max(a.dd);
rgdd=maxdd-mindd;

if ~isfield(a,'dd')
  s=ded_read_hdf([nm '/force/fx.hdf5']);
  a.dd=s.repmat(GUR,length(a.t));
end


sdd=(a.dd-mindd)./rgdd;

figure(1);
subplot(4,1,1);plot(c.Az,a.dd);line(c.Az,a.dd(:,end),'linewidth',2);ylabel('dd');
subplot(4,1,2);semilogy(a.t,max(eps,sum(abs(diff(a.dd)))./rgdd-1));ylabel('TVdd');
subplot(4,1,3);plot(a.t,w*sdd);ylabel('hu');
subplot(4,1,4);plot(a.t,w*(sdd.*(1-sdd)));ylabel('wu');
%subplot(4,1,4);plot(a.t,a.dd(1,:));ylabel('U0');
figure(2);
subplot(4,1,1);plot(c.Az,a.db);line(c.Az,a.db(:,end),'linewidth',2);ylabel('db');
subplot(4,1,2);semilogy(a.t,max(eps,sum(abs(diff(a.db)))./(max(a.db)-min(a.db))-1));ylabel('TVdb');
subplot(4,1,3);plot(a.t,w*a.db);ylabel('hb');
subplot(4,1,4);plot(a.t,w*(a.db.*a.dd));ylabel('bu');
figure(3);
subplot(4,1,1);plot(c.Az,b.dwdz);line(c.Az,b.dwdz(:,end),'linewidth',2);ylabel('dwdz');title(sprintf('%5.1e',max(abs(b.dwdz(:,end)))));
subplot(4,1,2);semilogy(c.Az,sqrt(max(0,b.ww)));line(c.Az,sqrt(max(0,b.ww(:,end))),'linewidth',2);ylabel('ww');set(gca,'ytick',10.^[-5:0]);title(sprintf('%5.1e',sqrt(max(b.ww(:,end)))));
%subplot(3,1,2);plot(b.t,sum(abs(b.dwdz)));ylabel('TVdwdz');
subplot(4,1,3);semilogy(b.t,sqrt(w*b.dwdz.^2),b.t,w*abs(b.dwdz));ylabel('Idwdz');
subplot(4,1,4);semilogy(b.t,sqrt(w*b.ww));ylabel('Iww');
figure(4);
ded_plot_X(nm);
% $$$ figure(5);
% $$$ ded_gc_2db(nm,{'b','psi'},'ay')
