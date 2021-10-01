function ded_gc_results
  dc;figure;[c f pp]=ded_gc_fit_g_Re('gc/f7/ma/*/[23]*',8,'rat34');
  
  
  dc;figure;[c f pp]=ded_gc_fit_g_Re('gc/f7/ma/*/*',8,'poly2');
  dc;figure;[c f pp]=ded_gc_fit_g_Re('gc/f7/mb/*/*',8,'poly2');

  dc;figure;[c f pp]=ded_gc_fit_g_Re('gc/f7/ma/*/*',8,'rat23');

  nm={'gc/f7/m/*/*','gc/f7/jb/*/*'};
  ded_gc_fit_g_Re(nm,8);
  ded_cmp_sims(nm);
  ded_plot_X(nm,-2);

  ded_plot_X('gc/f7/m/*/*',2);
  ded_time_profiles('gc/gc2d7n/73',3:18,50:1:90);

dc;figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/jb/*/*','X','Re');figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/jb/*/*','g','Re');

dc;figure;[c f pp]=ded_gc_fit_Z_Y_X('gc/f7/jb/*/*','Re','X','g');h=legend;delete(h);

ded_gc_fit_Z_Y_X('gc/f7/m/*/*','Re','X','g');axis('tight');



dc;figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/jb/*/*','X','Re');figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/jb/*/*','g','Re');
  

  dc;figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/l/05000/2*','X','g','quadratic');pp(end)=pp(end)-8;g=roots(pp);g=min(g);hold('on');plot(g,8,'s');title(sprintf('g=%8.5f',g));

dc;figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/l/0**00/2*','X','Re','quadratic');

ded_plot_X('gc/f7/l/0**00/2*',2


dc;

%g=2.282
%g=2.284
%g=2.283
%g=2.282
%g=2.286
%g=2.291
%g=2.300
%g-2.29
%g=2.295
%g=2.2966
%g=2.2957
nms=ded_regexp_names('gc/f7/m/*');
dc;[T dt]=ded_convergence_T(nms,[],1);
f=find(dt>500);
for j=1:length(f)
  k=f(j);
  nm=nms{k};
  p=ded_read_param(nm);
  disp(sprintf('%s: %10s T=%6.1f dt=%6.1f',nm,p.status,T(k),dt(k)))
  if strcmp(p.status,'Running')
    unix(sprintf('touch %s/abort',nm'));
  end
  if strcmp(p.status,'Aborted') | strcmp(p.status,'Submitted') 
    unix(sprintf('echo Terminated > %s/status',nm'));
  end
end

dc;
for j=2:5;
  figure;
  clf;
  ded_gc_fit_Y_X(sprintf('gc/f7/m/%02u000/*',j),'X','g','quadratic');
end




ded_gc_fit_Y_X('gc/f7/l/05000/*','X','g');axis('tight');


  
  pp(end)=pp(end)-8;  g=roots(pp);  g=min(g(g>2 & g<3));  hold('on');  plot(g,8,'s');  title(sprintf('g=%8.5f',g));
  
  
  f=find(abs(err)<0.05);
  
  figure;ded_gc_fit_Y_X('gc/f7/ia/*','g','Re','rat34');title('gc/f7/ia');

  figure;ded_gc_fit_Y_X('gc/f7/j/*','g','Re','rat34');title('gc/f7/j');
  ded_plot_X('gc/f7/ib/*',2,[],1);
end

figure;ded_gc_fit_Y_X('gc/f7/ib/*','X');title('gc/f7/ib');


figure(2);ded_plot_X('gc/f7/m/*');


figure(3);ded_gc_fit_Z_Y_X('gc/f7/[ij]*/*','nmb','g','Re');title('gc/f7/ia gc/f7/j');




figure(4);clf;ded_gc_fit_Y_X('gc/f7/l/*/*','X','W','linear');title('gc/f7/l X');
figure(5);clf;ded_gc_fit_Y_X('gc/f7/l/*/*','g','W','linear');title('gc/f7/l g');
figure(6);clf;ded_gc_fit_Y_X('gc/f7/l/*/*','X','g','linear');title('gc/f7/l Xg');


figure(4);clf;ded_gc_fit_Y_X('gc/f7/l/*/2*','X','g','linear');


figure(4);clf;ded_gc_fit_Y_X('gc/f7/M/*','X','Peb','linear');
title('gc/f7/M X');
figure(5);clf;ded_gc_fit_Y_X('gc/f7/M/*','g','Peb','linear');
title('gc/f7/M g');


figure(3);dc;ded_gc_fit_Z_Y_X('gc/f7/L/*/*','nmb','g','L');

figure(3);ded_gc_fit_Z_Y_X('gc/f7/L/*/*','Nz','g','L');


ded_cmp_sims('gc/f7/L/01/*')


dc;ded_plot_X('gc/f7/ia/*',2);
ded_plot_X('gc/f7/j/*',2);
ded_plot_X('gc/f7/l/*/*',2);

ded_plot_X('gc/f7/l/*/*',2);


ded_plot_X('gc/f7/L/01/*',2);


return;




ded_gc_fit_Y_X('gc/test/1*','g','L','linear');



dc;
ded_plot_X('gc/f7/g/*/*');
figure;ded_gc_fit_Z_Y_X('gc/f7/b/*/*','g','X','Re');axis([0 3000 4 9]);
figure;ded_gc_fit_Z_Y_X('gc/f7/b/*/*','Re','X','g');axis([1.6   2.1    4  9]);
 
nm='gc/gc2d7n/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.094 0.290 0.477 1.102 0.916 0.392 0.726 7.200]);


ded_gc_fit_Y_X('gc/f7/a/*','X','Re','linear');


axis([0 5000 4 9]);



nm='gc/gc2d7n/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','d2pow');figure(2);ded_gc_fit_Y_X(c,'X','Re');
nm='gc/gc2d7[an]/*';dc;figure(1);

ded_gc_fit_Y_X('gc/f6/f/2*','X','g','linear');
ded_gc_fit_Y_X('gc/f6/f/3*','X','g','linear');
ded_gc_fit_Y_X('gc/f6/f/0[23]*','X','g','linear');

ded_gc_fit_Y_X('gc/f6/f/03*','X','g','d2pow');
ded_gc_fit_Y_X('gc/f6/f/2*','X','g','linear');
ded_gc_fit_Y_X('gc/f7/*','X','g','linear');
ded_gc_fit_Y_X('gc/f7/g/*/*','X','g','linear');

ded_gc_fit_Z_Y_X('gc/f7/g/*/*','Re','X','g');



nm='gc/f6/f/03*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','d2pow');figure(2);ded_gc_fit_Y_X(c,'X','Re');

ded_plot_PID('gc/gc2d7n/*',-1,'Running');

clear('h');
for j=1:10
  h(j)=figure(j);clf;
end
for j=1:9
 ded_plot_PID({sprintf('gc/gc2d7n/%02u',j),sprintf('gc/gc2d7n/%02u',30+j)},h(j));
end


  ded_plot_PID('gc/gc2d7n/*',1);

nm='gc/f6/l/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'X','g','ipow');
ded_plot_PID('gc/f6/l/*',1);


nm='gc/gc2d7a/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.0286    0.3275    0.4550    1.1139    1.1942    0.2669    8.7332   29.9380]);
nm='gc/gc2d7k/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[0.8192    0.5125    0.3476    1.2600    0.0015    1.3331    0.3150    1.1832]);
nm='gc/gc2d7l/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.094 0.290 0.477 1.102 0.916 0.392 0.726 7.200]);
nm='gc/gc2d7n/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.094 0.290 0.477 1.102 0.916 0.392 0.726 7.200]);
nm='gc/gc2d7[lma]/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.0222    0.3348    0.4484    1.1035    0.7249    0.5848    0.5869   25.8326]);
nm='gc/gc2d7[am]/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.0222    0.3348    0.4484    1.1035    0.7249    0.5848    0.5869   25.8326]);
nm='gc/gc2d7[ak]/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[0.8073    0.5234    0.3408    1.0993    0.4343    0.8704    0.3708   34.5642]);
nm='gc/gc2d7[ae]/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.0077 0.3410 0.4463 1.1848 0.4748 1.7868 0.0926 26.8840]);
nm='gc/gc2d7[ha]/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[0.9657 0.3599 0.4412 1.0848 0.4821 1.7949 0.1486 26.8840]);
ded_plot_PID('gc/gc2d7a/*',1);
ded_plot_PID('gc/gc2d7e/*',1);
ded_plot_PID('gc/gc2d7h/*',1);
nm='gc/gc2d7h/*';dc;[c f pp]=ded_gc_fit_Y_X(nm,'g','Re','d2pow',[1.0220    0.3315    0.4527    1.2058    1.1574    0.2144    3.9324    6.3757]);

nm='gc/f6/f/0*';figure(1);ded_gc_fit_Y_X(nm,'X','Nx');
nm='gc/f6/f/2*';figure(2);ded_gc_fit_Y_X(nm,'X','Nx');
nm='gc/f6/f/3*';figure(3);ded_gc_fit_Y_X(nm,'X','Nx');

nm='gc/f6/f/*';figure(1);ded_gc_fit_Y_X(nm,'X','Nx');
ded_plot_PID('gc/f6/f/[023]*',1);

clf;figure(2);ded_gc_fit_Y_X('gc/f3/*','X','U');

ded_plot_PID('gc/f6/f/0*',1);

nms='[pa]gc/qgc[257]*/*';
nms='*gc/qgc*a/*';
nms='*gc/qgc*/*';

nms='pgc/qgc7a/*';
dbstop if error
dc;
nms='pgc/qgc7a/*';
figure(1);c=ded_gc_fit_Y_X(nms,'g','Re','dpow');
figure(2);ded_gc_fit_Y_X(c,'X','Re');
figure(3);ded_gc_fit_Y_X(c,'R','Re');
figure(4);ded_gc_fit_Y_X(c,'sX','Re');
figure(5);ded_gc_fit_Y_X(c,'sg','Re');
figure(6);ded_gc_fit_Y_X(c,'sR','Re');

nm='gc/f3/0*'
figure(1);c=ded_gc_fit_Y_X(nm,'X','U','dpow');
nm='gc/gc2d7[a]/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);

nm='gc/gc2d7[a]/*';dc;figure(1);
nm={'gc/gc2d7a/01','gc/gc2d7a/02','gc/gc2d7a/03','gc/gc2d7a/04','gc/gc2d7a/05','gc/gc2d7a/06','gc/gc2d7a/07','gc/gc2d7a/08'};
[c,f]=ded_gc_fit_Y_X(nm,'g','Re','dpow');
f(1e3*(1:19))

nm='cgc/f63/*';figure(1);c=ded_gc_fit_Y_X(nm,'X','Re');figure(2);ded_plot_PID(nm);ded_cmp_sims(nm);
ded_print(nm);
nm='pgc/qgc7d/*';figure(1);c=ded_gc_fit_Y_X(nm,'X','Re');figure(2);ded_plot_PID(nm);ded_cmp_sims(nm);
ded_print(nm);

nm='gc/f6/a/*';figure(1);c=ded_gc_fit_Y_X(nm,'X','Re');figure(2);ded_plot_PID(nm);ded_cmp_sims(nm);
nm='gc/f6/b/*';figure(1);c=ded_gc_fit_Y_X(nm,'X','Re');figure(2);ded_plot_PID(nm);ded_cmp_sims(nm);
nm='gc/f6/a/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);
nm='gc/f6/b/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);
nm='gc/f5/c/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);
nm='gc/f5/a/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);


ded_print(nm);

nm='gc/2df3/*';figure(1);c=ded_gc_fit_Y_X(nm,'U','g','ipow');figure(2);ded_gc_fit_Y_X(c,'X','g');figure(3);ded_plot_PID(nm);ded_cmp_sims(nm);
ded_print(nm);

nm='gc/f3/0*'
figure(1);c=ded_gc_fit_Y_X('pgc/qgc7*/*','g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID('pgc/qgc7*/*');

nm='gc/2df3/*'figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);

figure(1);c=ded_gc_fit_Y_X('gc/gc2d*/*','g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID('gc/gc2d*/*');



figure(1);clf;c=ded_gc_fit_Y_X('gc/gc2d7a/*','g','Re','dpow');figure(2);clf;a=array_select(c,1:11);dc;ded_gc_fit_Y_X(a,'g','Re','dpow');
figure(1);clf;c=ded_gc_fit_Y_X('gc/gc2d7[ah]/*','g','Re','dpow');figure(2);clf;a=array_select(c,13:20);ded_gc_fit_Y_X(a,'g','Re','dpow');
[c f pp]=ded_gc_fit_Y_X(c,'g','Re','d2pow',[1.012  0.337 0.449 1.186 0.302 1.618 0.119 26]);
nm='gc/gc2d7[aeh]/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_gc_fit_Y_X(c,'X','g');
ded_plot_PID('gc/gc2d7a/*',1);
ded_plot_PID('gc/gc2d7e/*',1);
ded_plot_PID('gc/gc2d7h/*',1);




nm='gc/f6/b/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);

nm='gc/gc2d7j/*';dc;figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);
ded_cmp_sims(nm);
ded_plot_PID('gc/gc2d7e/21');
nm='gc/gc2d7b/*';figure(1);c=ded_gc_fit_Y_X(nm,'g','Re','dpow');figure(2);ded_gc_fit_Y_X(c,'X','Re');figure(3);ded_plot_PID(nm);

ded_plot_PID('gc/qgc7a/*',1);

ded_print('gc/gc2d7a',1);

nms='pgc/qgc7a/*';
figure(1);c=ded_gc_fit_Y_X(nms,'g','Re','dpow');

figure(1);c=ded_gc_fit_Y_X('gc/gc2d*/*','g','Re','dpow');

ded_3d('gc/qgc2d'); 
ded_3d('gc/qgc2d'); 
a=ded_read_hdf('~/gc/qgc2c/01/noise/noise_s14.hdf5');
mesh(squeeze(a.noisey(:,1,:)))
subplot(2,1,1);
plot(a.z,sqrt(mean(mean(a.noisey.^2,3),2)))
subplot(2,1,2);
plot(a.z,sqrt(mean(mean(a.noisez.^2,3),2)))

ded_plot_PID('gc/qgc2d7a/*');
ded_print('gc/qgc2f/*',1);  

ded_plot_PID('gc/qgc2d*/*');

ded_gca_1('gc/qgc2a/08'); 
ded_gca_1('gc/qgc5a/08'); 
ded_gca_1('pgc/qgc7a/08'); 
ded_gca_1('pgc/qgc8a/08'); 
ded_gca_1('pgc/qgc10a/08'); 



ded_gc_X_U('5[23]*nms='gc/qgc*/*';
figure(1);c=ded_gc_fit_Y_X(nms,'g','Re');
','Re','X');


ded_gc_X_U('5[23]*','Re','X');

ded_gc_X_U('5[01]*','Re','X');  
dc;ded_gc_print('5[01]*',1);  

ded_gc_X_U('[78]*','Re','X');  



ded_cmp_sims({'[23478]*','5[01]*'});
ded_cmp_sims({'gc/[23478]*','gc/5[01]*'});

dc;
figure;ded_gc_fit_X_U('0*');
figure;ded_gc_fit_U_Re('9[0-2]*');
figure;ded_gc_fit_X_Re('5[01]*');
figure;ded_gc_fit_U_Re('5[23]*');
figure;ded_gc_fit_X_U('6*');
figure;ded_gc_fit_X_Re({'[23478]*','5[01]*'},'Re','X'); % Fixed U simulations forcing 5


figure;ded_gc_fit_U_Re('f6/*');

p.sz=[1024 512];
p.cva=29.8586;
p.aa=1;
ded_3d('pgc/f6/071','~/films/dedalus/pgc-f6-071',p);
p.aa=1;p.cva=29;p.sz=[1024 512];
ded_3d('gc/ccle/004','~/films/dedalus/gc-ccle-004',p);




Why does pm forcing 3 have non-zero divergence

x=squeeze(sum(sum((a.fw-permute(a.fv,[2 1 3])).^2,1),2));
x=squeeze(sum(sum((a.fu-permute(a.fu,[2 1 3])).^2,1),2));

for j=1:size(a.div,3)
  mesh(a.div(:,:,j));
  drawnow;
  pause;
end
a=ded_pm_plot_forcing('pm/pmf3');
div=a.fudx+a.fvdy+a.fwdz;
j=10;
subplot(2,2,1);mesh(a.fudx(:,:,j))
subplot(2,2,2);mesh(a.fvdy(:,:,j))
subplot(2,2,3);mesh(a.fwdz(:,:,j))
subplot(2,2,4);mesh(div(:,:,j))



[a p]=ded_gca_3('gc/gc2d7n/01');







