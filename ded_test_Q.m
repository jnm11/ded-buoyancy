

w=0.2;
Q=0.3;
H=2;
h=1+0.5*rand(1);
U=1;
x=linspace(0,H,100);
dx=x(2)-x(1);
f=-0.9-1.1*tanh((x-h)/w);
f=-1-tanh((x-h)/w);
f=-U*f/sum(f*dx);
plot(x,f);


F0=-sum(f*dx);
F1=sum(max(f,0)*dx);


g=ded_Q_iter(dx,f,U,H,Q);
plot(x,f,x,g1,x,g2);
