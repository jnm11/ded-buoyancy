function ded_bg_Re_Nx
R = 100:5:1000;

N1=repmat(2,size(R));
N2=repmat(inf,size(R));
T1=repmat(inf,size(R));
T2=repmat(inf,size(R));

N=145;
n=length(R);
k=0;
j=1;
!rm -rf ~/data/dedalus/bg/*
while j<=n
  k=k+1;
  nm=sprintf('bg/%04u',k);
  cmd=sprintf('ded_gc.py -g 0 -L 10 --xa 1 --xb 2 --dtu 0.1 --Re %u --Wx 1 -T 5 --Nx %u %s &> /dev/null',R(j),2*N,nm);
  unix(cmd);
  a=ded_read_g(nm,'u');
  mtv=max(diff(sum(abs(diff([a.u;a.u(1,:)],1,1)),1)*(a.x(2)-a.x(1))));
  disp(sprintf('k=%4u, j=%3u R=%4u, N=[%4u,%4u,%4u] mtv=[%8.1e,%8.1e,%8.1e] %s',k,j,R(j),2*N1(j),2*N,2*N2(j),T1(j),mtv,T2(j),nm));
  if mtv>0
    N1(j)=N;
    T1(j)=mtv;
  else
    N2(j)=N;
    T2(j)=mtv;
  end
  if N1(j)+1==N2(j)
    if mtv>0
      N=N2(j);
    end
    M(j)=N;
    j=j+1;
    N1(j)=N1(j-1);
    N=round(N*R(j)/R(j-1));
    continue;
  end
  if isfinite(N2(j))
    N=round((T1(j)*N2(j)-T2(j)*N1(j))/(T1(j)-T2(j)));
    N=min(N2(j)-1,N);
  else
    N=round(1.1*N);
  end
  N=max(N1(j)+1,N);
end
R=R(1:length(M));
p=polyfit(R,M,1);

plot(R,M,R,polyval(p,R));

save('~/ded_bg_Re_Nx.mat',R,M,p);
