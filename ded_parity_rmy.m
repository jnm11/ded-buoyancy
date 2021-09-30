function p=ded_parity_rmy(p)
fnm=fieldnames(p);
for j=1:length(fnm)
  a=fnm{j};
  if ndims(p.(a))==3
    p.(a)(2)==[];
  end
end
