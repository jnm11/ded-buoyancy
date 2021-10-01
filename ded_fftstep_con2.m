function [a b]=ded_fftstep_con2(f,n1,n2)

x=ifft(fft(f),n2,'symmetric')*n2/n1;
a=-x;%[-x x-1];
b=[];
