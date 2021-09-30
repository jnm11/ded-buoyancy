function b=ded_int(a,p,parity,c,nm)

if nargin<5
  nm=fieldnames(a);
end
T=flip(p.T);
X=flip(p.X);
dd=flip(c.dd);
cc=flip(c.c);
for j=1:length(nm)
  e=nm{j};
  for k=1:length(cc)
    f=[e 'I' cc(k)];
    if length(e)>3
      if strcmp(e([1 end-1:end]),['dd' cc(k)])
        f = e(2:end-2);
      end
    end
    switch(T{k})
      case 'Cheb'
        c=ichebintc(ichebf2c(a.(e),k),k,X(k));
        switch(k)
          case(1)
            c(end,:)=[];
         case(2)
            c(:,end,:)=[];
         case(3)
            c(:,:,end,:)=[];
         case(4)
            c(:,:,:,end,:)=[];
         case(5)
            c(:,:,:,end,:)=[];
        end
        b.(f)=ichebc2f(c,k,[],X(k));
      otherwise
        b.(f)=pr_int(a.(e),dd(k),k,[],parity.(nm{j})(4-k));
    end
  end
end
