%function
%rsync -uavp --progress tttuck:pm/f7/e/25 --exclude final --exclude b --exclude a --exclude force ~/pm/f7/e
%rsync -uavP tttuck:pm/f7/e/25/a/a-00014.hdf5 ~/pm/f7/e/25/a
%mkdir -p ~/pm/f7/e/27 ~/pm/f7/e/27/s ~/pm/f7/e/27/u 
%rsync -uavp --progress tttuck:pm/f7/e/27 --exclude final --exclude b --exclude s --exclude u --exclude avar --exclude sev --exclude a --exclude force --exclude yz --exclude ayz ~/pm/f7/e
%rsync -uavP tttuck:pm/f7/e/27/s/s-00207.hdf5 ~/pm/f7/e/27/s
%rsync -uavP tttuck:pm/f7/e/27/u/u-00207.hdf5 ~/pm/f7/e/27/u

nm='pm/f7/e/25';
ded_get_fns

fna='~/pm/f7/e/25/a/a-00014.hdf5';
p=ded_read_param(nm);
b=h5info(fna);
d={'u','uu','s','ss','su','b','bb','bu'};
dt=h5read(fna,'/dt');
for j=1:length(d)
  dd=d{j};
  a.(dd)=h5read(fna,['/' dd])/dt;
end
cs=(a.su-a.s.*a.u)./sqrt(a.ss.*a.uu);
cb=(a.bu-a.b.*a.u)./sqrt(a.bb.*a.uu);

subplot(1,2,1);hist(cs(a.s>1e-2),100)
subplot(1,2,2);hist(cb(a.b>1e-2),100)

[a fld b]=ded_read_javrg(nm,'a',[],'combine',{});




