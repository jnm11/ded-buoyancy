function f=ded_Q_iter(dx,f,U,H,Q)

F0=sum(f*dx);
F1=sum(max(f,0)*dx);
for j=1:10
  h=max(0.1*H,sum((f>0)*dx));
  A=(F1*H*U+F0*Q)/(F0*h-F1*H);
  B=(U*h+Q)*H/(-F0*h+F1*H);
  f=A+B*f;
  f=-H*U*f/sum(f*dx);
  F0=sum(f*dx);
  F1=sum(max(f,0)*dx);
  disp(sprintf('%u %5.3f %6.1e %6.1e %6.1e',j,F1,abs(F1-Q)/(H*U),H*U+F0,F1-Q));
  if abs(F1-Q)<H*U*1e-12
    break;
  end
end
