function a=ded_read_avrg(nm,T)

if nargin<2
  T=[];
end

a=ded_read_g(nm,'avrg');
p=ded_read_param(nm);

nt=length(a.t);
nms=fieldnames(a);

if nt>1
  if isempty(T)
    T=a.t(end-1);
  end
  f=min(find(a.t>=T));
  a.t=a.t([f end]);
  a.dt=diff(a.t);
  for j=1:length(nms)
    c=nms{j};
    x=a.(c);
    if ndims(x)~=3 
      continue
    end
    a.(c) = diff(a.(c)(:,:,[f end]),1,3)/a.dt;
  end
else
  if max(abs(a.ab(:)))==0
    a=[];
    return;
  end
  keyboard
end
