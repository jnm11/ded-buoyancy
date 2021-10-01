nm='pm/f4/03';
DD=ded_dedalus_data_dir;
nmf=[DD '/' nm '/final/final_s1.hdf5'];

p=ded_read_param(nm);
a=ded_pm_plot_forcing(nm);
b=ded_read_g(nm,'yz');
[s ss]=ded_read_state(nmf,p);
[maxu f]=max(max(max(a.fu,[],1),[],2),[],3);
x=(p.x2+p.x3)/2;
f=findmin(abs(a.x-x));

[yy zz]=ndgrid(a.y,a.z);
rr=sqrt(yy.^2+zz.^2);
F=2*exp(-2*rr.^2/p.radius.^2);

figure;mesh(a.y,a.z,a.fu(:,:,f)-p.U*F);
figure;mesh(a.y,a.z,a.fb(:,:,f)-p.B*F);

dA=(a.y(2)-a.y(1))*(a.z(2)-a.z(1));
A=pi*p.radius.^2;
B=dA*squeeze(sum(sum(a.fb,1),2))/A;
V=dA*squeeze(sum(sum(a.fu,1),2))/A;
Q=dA*squeeze(sum(sum(a.fu.*a.fb,1),2))/A;
M=dA*squeeze(sum(sum(a.fu.*a.fu,1),2))/A;
TA=p.H*p.W;

subplot(4,1,1);
plot(a.x,V,b.x,b.u/A)
ylabel('V/A');

subplot(4,1,2);
plot(a.x,B,b.x,b.b/A)
line([1;1]*[p.x4 p.x5 p.x6 p.x7],[0 3],'color',0.7*[1 0 0]);
axis([0 p.L 0 3]);
ylabel('B/A');

subplot(4,1,3);
plot(a.x,Q,b.x,b.ub/A)
ylabel('Q/A');

subplot(4,1,4);
plot(a.x,M,b.x,b.uu/A)
line([1;1]*[p.x0 p.x1 p.x2 p.x3],[0 3],'color',0.7*[1 1 1]);
ylabel('M/A');
axis([0 p.L 0 3]);


% $$$ plot(b.x,b.u/A);
% $$$ plot(b.x,b.b/A);
% $$$ b=ded_read_hdf('yz/yz_s6.hdf5');
% $$$ 
% $$$ a=ded_read_hdf('u_s10.hdf5');
% $$$ plot(a.x,squeeze(sum(sum(a.u,1),2))*(a.y(2)-a.y(1))*(a.z(2)-a.z(1)))

ss.x=a.x;
ss.y=a.y;
ss.z=a.z;


u=conj(cat(3,ss.u,zeros(p.Nz-1,p.Ny-1,1),conj(ss.u(:,:,end:-1:2))));
u=ifftn(ss.u,'symmetric')*p.Nx*p.Ny*p.Nz;
for j=1:p.Nx;mesh(s.u(:,:,j));drawnow;end;




      u(k,j,i)=exp(i*a.x(I)*ss.kx+a.y(J)*ss.ky+
for kx=1:length(kx)