

xx=linspace(-pi,pi,101);
yy=linspace(-pi,pi,100);

[y x]=ndgrid(yy,xx);

f=147/256-171.*cos(x).*(1/512)+21.*cos(2.*x).*(1/256)-5.*cos(3.*x).*(1/512)-171.*cos(y).*(1/512)-225.*cos(x).*cos(y).*(1/1024)+51.*cos(2.*x).*cos(y).*(1/512)-15.*cos(3.*x).*cos(y).*(1/1024)+21.*cos(2.*y).*(1/256)+51.*cos(x).*cos(2.*y).*(1/512)+3.*cos(2.*x).*cos(2.*y).*(1/256)-3.*cos(3.*x).*cos(2.*y).*(1/512)-5.*cos(3.*y).*(1/512)-15.*cos(x).*cos(3.*y).*(1/1024)-3.*cos(2.*x).*cos(3.*y).*(1/512)-(1/1024).*cos(3.*x).*cos(3.*y);
mesh(xx,yy,f');
