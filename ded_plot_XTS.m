function ded_plot_XTS(nm)


a=ded_read_g(nm,'yz');
aa=ded_read_g(nm,'xyz');
T=ded_convergence_T(nm);
b=ded_read_g(nm,'momb');
p=ded_read_param(nm);
t=a.t;
U=-aa.u/(p.L*p.H*p.W);
c=ded_gc_find_front(a,T);


mx=b.bx/b.b;
sx=sqrt(max(0,b.bxx/b.b-mx.^2));
X2=mx+sqrt(3)*sx;


clf;
trg=[min(min(a.t),min(c.t)) max(max(a.t),max(c.t))];

subplot(2,1,1);
mU=mean(U(floor(length(U)/2):end));
plot(aa.t,U,'r-',trg,mU*[1 1])
axis([trg -inf inf]);
title(sprintf('U=%8.6f',mU));
ylabel('U');
subplot(2,1,2);
plot(c.t,c.X,'b-',c.t([1 end]),p.PIDX*[1 1],'r--');
axis([trg -inf inf]);
ylabel('X');
xlabel('t');

return;

dt=median(diff(a.t));
s1 = iddata(c.X(:),u(:),dt);
s2 = detrend(s1);
plot(s2);

clf
m1 = impulseest(s2); % non-parametric (FIR) model
showConfidence(impulseplot(m1),3);

m1 = ssest(s2);

showConfidence(stepplot(m1,'b',3),3)


compare(s2,m1);  

m1 = nlarx(s2,[3 3 1],'sigmoidnet'); 
m2 = arx(s2,[2 2 1]); 
m3 = armax(s2,[2 2 1 2]); 
compare(s2,m1,'b',m2,'r',m3,'c');

m1 = nlarx(s2,[3 3 1],'sigmoidnet'); 
compare(s2,m1);





mi = impulseest(s1); % non-parametric (FIR) model
showConfidence(impulseplot(mi),3); %impulse response with 3 standard
                                   %deviations confidence region


m4 = ssest(s2);