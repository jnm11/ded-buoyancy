
%a=ded_read_hdf('gc/f6/f/25/final/state-00012.hdf5');
%a=ded_read_javrg('gc/f6/f/25','a');
c=ded_coord('gc/f6/f/25');
p=squeeze(a.p(1,:,:))/a.dt;
u=squeeze(a.u(1,:,:))/a.dt;
v=squeeze(a.v(1,:,:))/a.dt;
w=squeeze(a.w(1,:,:))/a.dt;
b=squeeze(a.b(1,:,:))/a.dt;
p=p-mean(p(:,end));
P=p+(u.^2+v.^2+w.^2)/2;
P=P-mean(P(:,end));
plot(c.x,mean(p,1),'-');

U=mean(u,1);
V=mean(v,1);
W=mean(w,1);
P=mean(p,1)+U.^2/2;
P=P-P(end);
plot(c.x,P,c.x,V.^2/2,c.x,W.^2/2);

P=squeeze(mean((a.p+(a.uu+a.vv+a.ww)/2)/a.dt,2));
P=P-P(1,end);



a=ded_read_hdf('gc/f6/f/25/a/a-00009.hdf5');
save('gc/f6/f/25/a/a-00009.mat', '-v7.3','a');

fn='~/gc/f6/f/20/a/a-00009.hdf5'
t1=h5read(fn,'/t1');
t2=h5read(fn,'/t2');
u=h5read(fn,'/u');
b=h5read(fn,'/b');
p=h5read(fn,'/p');
t2=h5read(fn,'/t2');
