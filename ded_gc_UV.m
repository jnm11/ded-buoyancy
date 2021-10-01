function a=ded_gc_UV(nm)
%rsync -vap hamilton:gc/f7/md/00100/27000 ~/gc/f7/md/00100 --exclude a --exclude b


% $$$ fns=cellstr_ls('~/gc/f7/*/*/*/sev/sev*.hdf5');
% $$$ 
% $$$ for j=1:length(fns)
% $$$   a=ded_read_hdf(fns{j});
% $$$   if isfield(a,'bp')
% $$$     unix(sprintf('/bin/rm -f %s',fns{j}));
% $$$   end
% $$$ end
% $$$ 
% $$$   nm=fns{j}(18:end);
% $$$  

if nargin==0
  nm='ma/*/*';
end

cd('~/');
fns=cellstr_ls([nm '/param.h5'],[],'dir');

clear a
j=0;
for jj=1:length(fns)
  nm=fns{jj};
  q=ded_read_param(nm);
  c=ded_coord(nm);
  s=ded_read_g(nm,'sev');
  if isempty(s);continue;end;
  j=j+1;
  if isfield(s,'divz')
    divz=s.divz(:,end);
  else
    y=ded_read_hdf([ded_dedalus_data_dir '/' nm '/force/fx.hdf5']);
    divz=y.GUR;
  end
  db=s.db(:,end);

  w=ichebintw(c.NAz)*q.H/2;
  [up un] = fit_serf(c.Az,divz,w);
  [bp bn] = fit_serf(c.Az,db,  w);
  a(j).u1 = up(1);  
  a(j).u2 = up(2);
  a(j).uh = up(3);  
  a(j).uw = up(4);
  a(j).b1 = bp(1);  
  a(j).b2 = bp(2);
  a(j).bh = bp(3);  
  a(j).bw = bp(4);
  
  subplot(2,1,1);plot(c.Az,divz,c.Az,un);xlabel('x');ylabel('u');
  subplot(2,1,2);plot(c.Az,db,  c.Az,bn);xlabel('x');ylabel('u');
  title(nm);
  a(j).Re=q.Re;
  a(j).g=q.g;
  a(j).nm=nm;
end
clf
subplot(4,2,1);plot([a.Re],[a.u1],'.');ylabel('u1');
subplot(4,2,3);plot([a.Re],[a.u2],'.');ylabel('u2');
subplot(4,2,5);plot([a.Re],[a.uh],'.');ylabel('uh');
subplot(4,2,7);plot([a.Re],[a.uw],'.');ylabel('uw');
subplot(4,2,2);plot([a.Re],[a.b1],'.');ylabel('b1');
subplot(4,2,4);plot([a.Re],[a.b2],'.');ylabel('b2');
subplot(4,2,6);plot([a.Re],[a.bh],'.');ylabel('bh');
subplot(4,2,8);plot([a.Re],[a.bw],'.');ylabel('bw');

return

  
  w=ichebintw(size(y,1));
U=0*c.Az;
U(c.Az<=mean(h))=mean(U1);
U(c.Az>mean(h))=mean(U2);
plot(c.Az,y,c.Az,U);

p=polyfit(c.Az,mean(y,2),10);
plot(c.Az,y,c.Az,polyval(p,c.Az));

[p d]=fit_erf(c.Az,mean(y,2));
plot(c.Az,y,c.Az,d);

U1=p(1)-p(2);
U2=p(1)+p(2);
h=p(4);
w=1/p(3);
