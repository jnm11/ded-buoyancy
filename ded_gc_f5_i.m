function ded_gc_f5_i(a);
nm=['gc/f5/i/' a];
c=ded_coord(nm);
p=ded_read_param(nm);
fn=ded_get_fn(nm,'b');
a=ded_read_hdf(fn{end});

f=max(find(a.b(1,:)>0.1)):length(c.x);
x=min(c.x(f(diff(a.b(1,f))>0)));

bs1=2*sqrt(mean(a.b(1,c.x<p.x3).^2));
bs2=max(abs(a.b(1,c.x>x)));
bs=max(bs1,bs2);
figure(1);clf;
subplot(4,1,1);plot(c.x,a.b(1,:)-1,'s-');axis([p.x3 p.x3+0.5 -bs bs]);title(sprintf('%6.1e',bs));
subplot(4,1,2);plot(c.x,a.b(1,:),'s-');axis([0 p.x3*1.1 -bs bs]);
subplot(4,1,4);plot(c.x,a.b(1,:)-1,'s-');axis([x-0.5 x -bs +bs]);
subplot(4,1,3);plot(c.x,a.b(1,:),'s-');axis([x p.L -bs bs]);
nm=['~/' nm];
wb=ded_read_hdf([nm '/force/wb.hdf5']);
wu=ded_read_hdf([nm '/force/wu.hdf5']);
fb=ded_read_hdf([nm '/force/fb.hdf5']);

figure(2);clf;
subplot(3,1,1);plot(c.x,fb.fb(1,:),'s-',c.x,wb.wb(1,:),'^-',c.x,wu.wu(1,:),'v-');axis([0 3 0-1e-2 1+1e-2]);
subplot(3,1,2);plot(c.x,fb.fb(1,:),'s-',c.x,wb.wb(1,:),'^-',c.x,wu.wu(1,:),'v-');axis([0 3 0-1e-4 0+1e-4]);
subplot(3,1,3);plot(c.x,fb.fb(1,:),'s-',c.x,wb.wb(1,:),'^-',c.x,wu.wu(1,:),'v-');axis([0 3 1-1e-4 1+1e-4]);

