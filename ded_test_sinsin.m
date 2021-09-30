function ded_test_sinsin(num)

%A="-L 2 -H 1 --Nx 256 --Ny 1 -g None --Re 500 -T 10 -U 1"
%B="--dtjp 0.1 --dtju 0.1 --dtjw 0.1 -B None -S None " 
%rm -rf ~/gc/sinsin/010 ; mpiexec -n 4 ded_gc.py $A $B  --fuspecial1 4 --fuspecial2 2                             gc/sinsin/010
%rm -rf ~/gc/sinsin/011 ; mpiexec -n 4 ded_gc.py --reset --pfn sinsin/010 --rfn sinsin/010 -T 20 --fuspecial None gc/sinsin/011
%rm -rf ~/gc/sinsin/001 ; mpiexec -n 4 ded_gc.py $A $B                                                            gc/sinsin/001
%rm -rf ~/gc/sinsin/002 ; mpiexec -n 4 ded_gc.py --reset --pfn sinsin/001 --rfn sinsin/001 -T 50 --fuspecial None gc/sinsin/002
%rm -rf ~/gc/sinsin/003 ; mpiexec -n 4 ded_gc.py --reset --pfn sinsin/001 --rfn sinsin/001 -T 50 --fuspecial None --Scs 1 -S 1 --sinitial 1 --Skp 0.01 --dtjs 0.1 gc/inertial/003
%num='002';

nm=['gc/inertial/' num];
u=ded_read_g(nm,'u');
w=ded_read_g(nm,'w');
p=ded_read_g(nm,'p');
c=ded_coord(nm);
pp=ded_read_param(nm);


U=pp.U;
H=pp.H;
L=pp.L;

kx=2*pi/L;
kz=2*pi/H;
kk=sqrt((kx^2+kz^2)/2);
[xx zz]=ndgrid(kx*c.Jx,kz*c.Jz);
psi = U*sin(xx).*sin(zz)/kk;
fu = +U*sin(xx).*cos(zz)*kz/kk;
fw = -U*cos(xx).*sin(zz)*kx/kk;
fp = U^2*(kz^2*cos(2*xx) + kx^2*cos(2*zz))/kk/4;
fq = U^2*(     cos(2*xx) +      cos(2*zz))*kx^2*kz^2/kk^2;
fww= 2*U*sin(xx).*sin(zz)*kk;
fss= U^2/kk^2*(2*kx^2*kz^2*cos(xx).^2.*cos(zz).^2+sin(xx).^2.*sin(zz).^2*(kx^2-kz^2)^2/2);
[f g ww d ss rr q]=helmholtz_fft(fu,fw,c.dJx,c.dJz);
[fe fc fcx fcy]=errcorr(f,psi);
[qe qc qcx qcy]=errcorr(q,fq);
[we wc wcx wcy]=errcorr(ww,fww);
[re rc rcx rcy]=errcorr(rr,fww.^2/2);
[se sc scx scy]=errcorr(ss,fss);
disp([fe qe we re se]);
fu=fu';
fw=fw';
fp=fp';

ar=pp.H/pp.L;

figure(1);clf;
ah=jsubplot([1 2],[0.1 0.1],[0.05 0.05],[0.05 0.05]);
axes(ah(1));
imagesc(c.Jx,c.Jz,fp);
hold('on');
cnts=contourc(c.Jx,c.Jz,psi',20);
ch1=plot_contours(cnts);
set(ah(1),'dataaspectratio',[1 1 1],'ydir','normal');
gu_set_aspect_ratio(ah(1),ar);
axis([0 pp.L 0 pp.H]);
xlabel('x');
ylabel('z');
[hc hp hl]=jcolorbar(ah(1),[20 20],'p');



axes(ah(2));cla;
imagesc(c.Jx,c.Jz,fq);
hold('on');
cnts=contourc(c.Jx,c.Jz,psi',20);
ch1=plot_contours(cnts);
set(ah(2),'dataaspectratio',[1 1 1],'ydir','normal');
gu_set_aspect_ratio(ah(2),ar);
axis([0 pp.L 0 pp.H]);
xlabel('x');
ylabel('z');
[hc hp hl]=jcolorbar(ah(2),[20 20],'q');

fus=sqrt(mean(mean(fu.^2)));
fws=sqrt(mean(mean(fw.^2)));
fps=sqrt(mean(mean(fp.^2)));

u.t(1)=[];
u.u(:,:,1)=[];
w.w(:,:,1)=[];
p.p(:,:,1)=[];

t=u.t;
uc=squeeze(mean(mean(fu.*u.u,1),2));
wc=squeeze(mean(mean(fw.*w.w,1),2));
pc=squeeze(mean(mean(fp.*p.p,1),2));

us=sqrt(squeeze(mean(mean(u.u.^2,1),2)));
ws=sqrt(squeeze(mean(mean(w.w.^2,1),2)));
ps=sqrt(squeeze(mean(mean(p.p.^2,1),2)));

figure(2);clf;
tt=linspace(t(1),t(end),1e3);
e=exp(-(kx^2+kz^2)/pp.Re*tt);
ah=jsubplot([1 3],[0.1 0.1],[0.05 0.05],[0.05 0.05]);
axes(ah(1));plot(t,uc,tt,e*uc(1));ylabel('u');   axis([t(1) t(end) 0 max(uc)]);
axes(ah(2));plot(t,wc,tt,e*wc(1));ylabel('w');   axis([t(1) t(end) 0 max(wc)]);
axes(ah(3));plot(t,pc,tt,pc(1)*e.^2);ylabel('p');axis([t(1) t(end) 0 max(pc)]);
xlabel('t');

figure(3);clf;
e=exp(-(kx^2+kz^2)/pp.Re*t');
ah=jsubplot([1 3],[0.1 0.1],[0.05 0.05],[0.05 0.05]);
axes(ah(1));plot(t,uc./(e   .*uc(1))-1);ylabel('u');gu_setalim();
axes(ah(2));plot(t,wc./(e   .*wc(1))-1);ylabel('w');gu_setalim();
axes(ah(3));plot(t,pc./(e.^2.*pc(1))-1);ylabel('p');gu_setalim();
xlabel('t');



figure(4);clf;
ah=jsubplot([1 3],[0.1 0.1],[0.05 0.05],[0.05 0.05]);
axes(ah(1));plot(t,uc./us/fus-1);ylabel('u');axis('tight');
axes(ah(2));plot(t,wc./ws/fws-1);ylabel('w');axis('tight');
axes(ah(3));plot(t,pc./ps/fps-1);ylabel('p');axis('tight');
xlabel('t');



return;


ded_2d('gc/inertial/003',{'s'});
