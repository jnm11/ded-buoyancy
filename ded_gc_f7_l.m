rsync -uavp --progress asahi:gc/f7/l/05000 --exclude final  --exclude b --exclude force --exclude a --exclude "initial*" ~/gc/f7/l
rsync -uavp --progress asahi:gc/f7/l/05000 --exclude final  --exclude b --exclude force --exclude a --exclude y --exclude ay --exclude yz --exclude ayz --exclude "initial*" ~/gc/f7/l

nm='gc/f7/l/05000/2250';
pm=ded_read_param(nm);
c=ded_coord(nm);
[a fld]=ded_read_javrg(nm,'ay',[20 inf]);
for j=1:length(fld)
    a.(fld{j})=a.(fld{j})/pm.W;
end
[a b]=ded_gc_diff_int(a,p,c.dJx);
Ib=max(0,a.bIz(end,rgx));

uu=a.uu([1 end],rgx);
vv=a.vv([1 end],rgx);
ww=a.ww([1 end],rgx);
uv=a.uu([1 end],rgx);
uw=a.vv([1 end],rgx);
vw=a.ww([1 end],rgx);
p =a.p([1 end],rgx);


rgx=find(c.Jx>0.2 & c.Jx<pm.L-0.2);
x=c.Jx(rgx);
rgb=max(a.b(:,rgx))>1e-3;
rgb(max(find(rgb)):length(rgb))=0;
rgb(1:min(find(rgb))-1:1)=1;
Fr=pm.U./sqrt(max(0,pm.g*Ib));
subplot(2,1,1);
plot(x,Fr,x,Ib);
axis([0 pm.L 0 1]);
P = p+uu/2;
P = p+uu + a.duwdzIx([1 end],rgx);
P=P-P(:,end);
subplot(2,1,2);
h=plot(x,P,x,a.dudx([1 end],rgx)/pm.Re);
h=legend(h,'bottom','top');




plot(x,b.cpx(end,rgx))



h=plot(x,uu,x,vv,x,ww,x,uv,x,uw,x,vw);
legend(h,{'uu','vv','ww','uv','uw','vw'});

