function ded_gc3(n)
if nargin<1
  fns=sort(cellstr_ls('~/data/dedalus/gc/*/*.h5'));
  fn=fns{end};
else
  fn=sprintf('~/data/dedalus/gc/gc_s%i/gc_s%i_p0.h5',n,n);
end

if ~isfile(fn)
  disp(fn);
  return;
end
disp(fn);
h5disp(fn);
%a=h5info(fn);
a.x=h5read(fn,'/scales/x/2');
a.y=h5read(fn,'/scales/z/2');
a.dt=h5read(fn,'/scales/timestep');
a.t=h5read(fn,'/scales/sim_time');
a.i=h5read(fn,'/scales/iteration');
a.b=h5read(fn,'/tasks/b');
a.u=h5read(fn,'/tasks/u');
a.w=h5read(fn,'/tasks/w');

H=a.y(end)-a.y(1);
L=a.x(end)-a.x(1);
S=1024/L;

nx=length(a.x);
ny=length(a.y);
nt=length(a.t);
a.b=reshape(a.b,[ny nx nt]); 
a.u=reshape(a.u,[ny nx nt]); 
a.w=reshape(a.w,[ny nx nt]); 
figure(1);
clf;
y=linspace(a.y(1),a.y(end),round(H/L*length(a.x)));

for j=1:length(a.t)
  b=interp1(a.y,a.b(:,:,j),y,'linear');
  u=interp1(a.y,a.u(:,:,j),y,'linear');
  v=interp1(a.y,a.v(:,:,j),y,'linear');
  ts=sprintf('t=%5.4f, i=%4u',a.t(j),a.i(j));
  if j==1
    gcf;
    hi=imagesc(a.x,y,b);
    set(gca,'position',[0 0 1 1],'ydir','normal','dataaspectratio',[1 1 1]);
    set(gcf,'units','pixels');
    p=get(gcf,'position');
    set(gcf,'position',[p(1:2),round(S*L),round(S*H)]);
    ht=text(L/200,H,ts,'horizontalalignment','left','verticalalignment','top','color',.95*[1 1 1]);
  else
    set(hi,'CData',b);
    set(ht,'string',ts);
  end
  drawnow;
  pause(0.01);
end
