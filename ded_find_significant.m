function  [p r ss]=ded_find_significant(E,w)
% Find the best rms residual for every subset
m=size(E,1);
n=size(E,2);
o=size(E,3);

ss=sqrt(squeeze(mean(sum(w.*E.^2,1),2)));
s=findmax(ss);

oo=(1:o);
oo(s)=[];

for j=1:o
  f=nchoosek(oo,j-1);
  n=size(f,1);
  rr=zeros(n,1);
  for k=1:n
    rr(k)=sqrt(mean(sum(w.*sum(E(:,:,[s f(k,:)]),3).^2,1),2));      
  end
  [r(j) k]=min(rr);
  p{j}=sort([s f(k,:)]);
end
