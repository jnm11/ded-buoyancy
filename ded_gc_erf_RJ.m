function a=ded_gc_erf_RJ(a)
for j=1:length(a)
  a(j).R=a(j).uw./a(j).bw;
  a(j).du=abs(a(j).u1-a(j).u2);
  a(j).db=abs(a(j).b1-a(j).b2);
  a(j).dudz=a(j).du./a(j).uw/sqrt(pi);
  a(j).dbdz=a(j).db./a(j).bw/sqrt(pi);
  a(j).Rib=a(j).g.*a(j).db.*a(j).uh./a(j).du.^2;
  a(j).J=a(j).g.*a(j).dbdz./a(j).dudz.^2;
  a(j).Fr=a(j).du./sqrt(a(j).g.*a(j).bh);
end
