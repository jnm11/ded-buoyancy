%rsync -vap --progress tttuck:pm/f7/e/5* ~/pm/f7/e --exclude force --exclude [absu] --exclude final --exclude sev --exclude yz --exclude avar 
%rsync -vap --progress tttuck:pm/f7/e/57/s/s-00112.hdf5 ~/pm/f7/e/57
nm='pm/f7/e/57';
trg=[0 inf];
p=ded_read_param(nm);
s=ded_read_stats(nm);
c=ded_coord(nm);
yz=ded_read_javrg(nm,'ayz',trg);
 

btol=max(yz.b(:))/1000;
stol=max(yz.s(:))/1000;
utol=max(yz.u(:))/1000;
bbtol=max(yz.bb(:))/1000;
sstol=max(yz.ss(:))/1000;
uutol=max(yz.uu(:))/1000;

yz.rb(yz.b<btol)=NaN;
yz.rs(yz.s<stol)=NaN;
yz.ru(yz.u<utol)=NaN;
rb=max(0,yz.b)./sqrt(max(bbtol,yz.bb));
rs=max(0,yz.s)./sqrt(max(sstol,yz.ss));
ru=max(0,yz.u)./sqrt(max(uutol,yz.uu));
Ab=max(0,yz.b).^2./max(bbtol,yz.bb);
As=max(0,yz.s).^2./max(sstol,yz.ss);
Au=max(0,yz.u).^2./max(uutol,yz.uu);
rb(yz.b<btol)=NaN;
rs(yz.s<stol)=NaN;

plot(c.Jx,yz.rb,c.Jx,yz.rs,c.Jx,yz.ru,c.Jx,rb,c.Jx,rs,c.Jx,ru);
axis([0 p.L 0 inf]);


subplot(2,1,1);
plot(c.Jx,yz.rb./rb,c.Jx,yz.rs./rs,c.Jx,yz.ru./ru);
axis([0 p.L 0.6 0.7]);

subplot(2,1,2);
plot(c.Jx,yz.Ab./Ab,c.Jx,yz.As./As,c.Jx,yz.Au./Au);
axis([0 p.L 0.45 0.6]);

s=ded_read_hdf('~/pm/f7/e/57/s/s-00112.hdf5');

clf;
hist(s.s(:),100)
set(gca,'yscale',log);

[w f xg yg zg]=jgrid(x,y,dx,typ,pr,iminx,imaxx,nlz);

mins=min(s.s(:));
maxs=max(s.s(:));
h=linspace(mins,maxs,100);
[a b]=jhist(h,s.s(:));



p=percentile(s.s(:),[0.1 99.9]);


for j=1:size(s.s,3)
  imagesc(s.s(:,:,j));
  set(gca,'clim',p);
  drawnow;
end
