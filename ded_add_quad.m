function a=ded_add_quad(a,nn)
nn=setdiff(nn,fieldnames(a));
for j=1:length(nn)
  n=nn{j};
  n1=n(1);
  n2=n(2:end);
  a.(n)=a.(n1).*a.(n2);
end
