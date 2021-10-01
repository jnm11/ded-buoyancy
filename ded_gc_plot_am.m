function ded_gc_plot_am(nm)
%ded_gc_plot_am(nm) % look at angular momentum

nm='gc/qgcf6';
p=ded_read_param(nm);
a=ded_read_g(nm,'xyz');
b=ded_read_g(nm,'momb');
f=ded_read_g(nm,'force');


X = b.bx./b.b;
Z = b.bz./b.b;
U = a.bu./a.b;
W = a.bw./a.b;
XU= b.bux./a.bu;
ZU= b.buz./a.bu;
XW= b.bwx./a.bw; % These are poorly posed
ZW= b.bwz./a.bw;

w   = (b.buz-b.bwx)/b.b;
w0  = (a.bu.*Z-a.bw.*Z)/b.b;

subplot(4,1,1);
plot(b.t,X);
ylabel('bx');
axis('tight')
subplot(4,1,2);
plot(b.t,Z);
ylabel('bz');
axis('tight')
subplot(4,1,3);
plot(b.t,w0-w);
ylabel('angular momentum');
axis('tight')
subplot(4,1,4);
plot(a.t,a.Ex+a.Ey+a.Ez);
ylabel('Enstrophy');
axis('tight')


