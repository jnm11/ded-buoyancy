function e=ded_mix_entropy(x)
e=zeros(size(x));
f=find(x>0 & x<1);
y=x(f);
e(f)=-y.*log(y)-(1-y).*log(1-y);


