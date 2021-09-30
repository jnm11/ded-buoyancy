function x=ded_parity(x,P,edge)

if nargin<3
  edge=[];
end
if isempty(edge)
  edge=true;
end
if edge
  for j=1:length(P)
    switch(P(j))
      case -1
        x=cat(j,x,flip(-x,j));
      case 0
      case 1
        x=cat(j,x,flip(x,j));
    end
  end
else
  for j=1:length(P)
    switch(P(j))
      case -1
        x=cat(j,flip(-x,j),x);
      case 0
      case 1
        x=cat(j,flip(x,j),x);
    end
  end
end
