
fn='~/Dropbox/Jim-Claudia-GC/mat/mat-ccle-024/024'
nm='-22-42';

fn='~/Dropbox/Jim-Claudia-GC/mat/mat-ccle-046/046'
nm='-20';

load([fn  nm '-steady-cheb.mat']);
load([fn '-param.mat']);
w=ichebintw(p.Nz,p.H);

bw=w*a.bw;
bu=w*a.bu;
b=w*a.b;
btol=0.01;
subplot(3,1,1);plot(a.x,b);
subplot(3,1,2);plot(a.x,bw./max(btol,b));
subplot(3,1,3);plot(a.x,bu./max(btol,b));
