function ded_ccle_tavrg
fnm={'Ex','Ey','Ez','S','b','p','u','uu','uv','uw','v','vv','vw','w','ww'};
nms={'gc-ccle-020','gc-ccle-021','gc-ccle-022','gc-ccle-024','gc-ccle-025','gc-ccle-046'};
kk= {4:9,           4:9,         4:9,           4:9,          4:9};          
for j=1:length(nms)
  nm=nms{j};
  fno=sprintf('~/data/dedalus/%s.mat',nms{j});
  if isfile(fno)
    continue;
  end
  T=0;
  for k=kk{j}
    nm=sprintf('~/data/dedalus/%s-%02u.mat',nms{j},k);
    b=load(nm);
    if k==4
      a=b.a;
    else
      dt=a.t2-a.t1;
      T=T+dt;
      for jj=1:length(fnm)
        a.(fnm{jj})=a.(fnm{jj})+dt*b.a.(fnm{jj});
      end
    end
    a.t1=min(a.t1,b.a.t1);
    a.t2=max(a.t2,b.a.t2);
  end
  for k=1:length(fnm)
    a.(fnm{k})=a.(fnm{k})/T;
  end
  save(fno,'a');
end
