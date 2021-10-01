
jnm@tokachi f]$ rsync -vap tttuck:gc/f6/f ~/gc/f6 --exclude final --exclude a --exclude b

rsync -vap tuck:gc/f6/f/[23]* ~/gc/f6/f --exclude final --exclude a --exclude b


nm='gc/f6/f/20';
c=ded_coord(nm);
cd(['~/' nm '/force']);

wu=ded_read_hdf('wu.hdf5');wu=squeeze(mean(wu.wu,2));
wb=ded_read_hdf('wb.hdf5');wb=squeeze(mean(wb.wb,2));
fu=ded_read_hdf('fu.hdf5');fu=squeeze(mean(fu.fu,2));
fw=ded_read_hdf('fw.hdf5');fw=squeeze(mean(fw.fw,2));
fb=ded_read_hdf('fb.hdf5');fb=squeeze(mean(fb.fb,2));
fx=ded_read_hdf('fx.hdf5');
w=ichebintw(size(wu,1));

subplot(3,2,1);imagesc(c.Ax,c.Az,wu);title('wu');set(gca,'ydir','normal');
subplot(3,2,2);imagesc(c.Ax,c.Az,wb);title('wb');set(gca,'ydir','normal');
subplot(3,2,3);imagesc(c.Ax,c.Az,fu);title('fu');set(gca,'ydir','normal');
subplot(3,2,4);imagesc(c.Ax,c.Az,fb);title('fb');set(gca,'ydir','normal');
subplot(3,2,5);imagesc(c.Ax,c.Az,fw);title('fw');set(gca,'ydir','normal');
subplot(3,2,6);plot(c.Az,fu(:,1:10:end));

ay=ded_read_hdf('ay/ay-00018.hdf5');



rgx=find(c.x>30 &c.x<45);
imagesc(c.x(rgx),c.z,ay.SR(:,rgx));set(gca,'ydir','normal');
rsync -vap poseidon:gc/f6/f/[23]* ~/gc/f6/f --delete --exclude final --exclude b --exclude a

