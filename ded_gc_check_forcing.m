function ded_gc_check_forcing(nm)

nm='gc/test/04';
a=ded_read_g(nm,'avrg')
p=ded_read_param(a);
nd=ndims(a.u)-1;
if nd==2
  c.u=squeeze(a.u(:,end,:));
  c.w=squeeze(a.w(:,end,:));
  c.b=squeeze(a.b(:,end,:));
else
  c.u=squeeze(a.u(:,:,end,:));
  c.v=squeeze(a.v(:,:,end,:));
  c.w=squeeze(a.w(:,:,end,:));
  c.b=squeeze(a.b(:,:,end,:));
  c.y=a.y;
end
c.z=a.z;
c.x=a.x(end);
c=ded_zgrid(c,2*a.nz,[],{},{},{});
subplot(nd+1,1,1);plot(c.z,c.b);ylabel('b');axis('tight');
subplot(nd+1,1,2);plot(c.z,c.u);ylabel('u');axis('tight');
subplot(nd+1,1,nd+1);plot(c.z,c.w);ylabel('w');axis('tight');
if nd==3
  subplot(nd+1,1,3);plot(c.z,c.v);ylabel('v');axis('tight');
end
dz=c.z(2)-c.z(1);
if nd==2
  e.u=sum((c.u-p.U).^2,1)*dz;
  e.w=sum(e.w.^2,1)*dz;
  e.b=sum(e.b.^2,1)*dz;
else
  dy=a.y(2)-a.y(1);
  e.u=squeeze(sum(sum((c.u-p.U).^2,1),2)*dz*dy);
  e.w=squeeze(sum(sum(e.w.^2,1),2)*dz*dy);
  e.v=squeeze(sum(sum(e.v.^2,1),2)*dz*dy);
  e.b=squeeze(sum(sum(e.b.^2,1),2)*dz*dy);
end

