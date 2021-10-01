function u=ded_direct_transform(v,s)
u=0;
for I=1:length(s.kz)
  disp(sprintf('%d/%d',I,length(s.kz)));
  for J=1:length(s.ky)
    for K=1:length(s.kx)
      [Z Y X]=ndgrid(exp(i*s.kz(I)*s.z),exp(i*s.ky(J)*s.y),exp(i*s.kx(K)*s.x));
      u=u+v(I,J,K)*X.*Y.*Z;
    end
  end
end
