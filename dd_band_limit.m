


f1 = @(x,a)->1./(1+exp(-x/a));
f2 = @(x,a)->(1+erf(-x/a))/2;


a=1;
x=-linspace(-5,5);

plot(x,f1(x,a),x,f2(x,a));

