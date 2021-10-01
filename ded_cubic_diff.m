% test function to demonstrate cubic diffusion

n=20;
f=randn(n,1).^2;
f([1:2 end-1:end])=0;
df=zeros(n+1,1);
rgf=[min(f) max(f)];
x=linspace(0,1,n);
while(1)
  %df(2:end-1)=diff(f);
  %df=diff([0 f([1 1:end end]));
  df=diff([0;f;0]);
  dfs=df.^2;
  maxdf=max(dfs);
  f=f+diff(dfs.*df/(4*maxdf));
  plot(x,f,'-s');
  axis([0 1 min(f) max(f)]);
  title(sprintf('%8.3f',mean(f)));
  drawnow;
  pause;
end

