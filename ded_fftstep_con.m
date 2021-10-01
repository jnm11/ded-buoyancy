function [a b]=ded_fftstep_con(f,n1,n2,s1,s2)

x=ifft(fft(f),n2,'symmetric')*n2/n1;
a=[-x x-1 s1.*diff(f([end 1:end 1]),2,2) s2.*diff(x([end 1:end 1]),2,2)];
b=[];
