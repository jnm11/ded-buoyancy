function f=ded_resample(f,A,P)

f=ichebf2c(f,1);
f=f(1:size(f,1)*A,:);
f=ichebc2f(f,1);
f=fftresample(f,size(f).*([1 A]), [0 P]);

