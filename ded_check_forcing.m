function ded_check_forcing(nm)
nm='gc/qgcf6';
s=ded_read_stats(nm);
p=ded_read_param(nm); 
u=ded_read_g(nm,'u'); 
v=ded_read_g(nm,'v'); 
w=ded_read_g(nm,'w'); 
t=u.t;
x=u.x;
y=u.y;
z=u.z;
u=u.u(:,:,1,:)+p.U;
v=v.v(:,:,1,:);
w=w.w(:,:,1,:);
E=squeeze(sqrt(sum(sum(u.^2+v.^2+w.^2,1),2)/(size(u,1)*size(u,2))));
[ts s r]=ded_interp_stats(nm,t);
plot(t,E,t,r.UE);

