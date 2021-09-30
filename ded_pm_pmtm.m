function b=ded_pm_pmtm(nm,b,qq,dd)


fns=ded_get_fn(nm,'b');
for j=1:Length(fns)
  disp('
  b=ded_read_hdf(fn{j});
  if b.t<qq.rgt(1) | b.t>qq.rgt(2); continue; end
  pb(:,j)=squeeze(mean(abs(fft(b.b,[],1)).^2,2))+squeeze(mean(abs(fft(b.b,[],2)).^2,1));
end

b=rmfieldifexist(b,{'u','v','p','w'});
c=ded_coord(nm);
p=ded_read_param(nm);
% Do not flip in the 'x' direction
b.s=ded_parity(b.s,flip(b.s_parity).*[1 1 0],false);
b.b=ded_parity(b.b,flip(b.b_parity).*[1 1 0],false);
if strcmp(p.Ty,'SinCos') c.y=[c.y-p.W;c.y]; end
if strcmp(p.Tz,'SinCos') c.z=[c.z-p.H;c.z]; end

m=findmax(squeeze(max(max(b.s))));

fx=find(c.x>qq.rgx(1) & qq.rgx(2));
ps=squeeze(mean(mean(abs(fft(b.s(:,:,fx),[],1)).^2,2),3))+squeeze(mean(mean(abs(fft(b.s(:,:,fx),[],2)).^2,1),3))';
pb=squeeze(mean(mean(abs(fft(b.b(:,:,fx),[],1)).^2,2),3))+squeeze(mean(mean(abs(fft(b.b(:,:,fx),[],2)).^2,1),3))';
k=imag(fft_modes(length(pb),c.dy));
f=find(k>0);

figure(1);
clf;
w=6.88;
preprint([w 2.10],8);colour_lines;
ah=jsubplot([3 1],[0.1 0.2],[0.05 0.05],[0.02 0.05]);

axes(ah(1));
imagesc(c.y,c.z,max(0,2*b.s(:,:,m)));
set(gca,'clim',[0 6],'ydir','normal');
axis([-2 2 -2 2]);
text(-1.9,1.9,'a)','horizontalalignment','left','verticalalignment','top','fontweight','bold','color',[1 1 1]);
ylabel('$y$','interpreter','latex');
xlabel('$x$','interpreter','latex');
axes(ah(2));
imagesc(c.y,c.z,max(0,b.b(:,:,m)));
axis([-2 2 -2 2]);
set(gca,'clim',[0 1],'ydir','normal');
text(-1.9,1.9,'b)','horizontalalignment','left','verticalalignment','top','fontweight','bold','color',[1 1 1]);
xlabel('$x$','interpreter','latex');

%set(ah,'dataaspectratio',[1 1 1],'box','on');


axes(ah(3));
h=loglog(k(f)/pi,ps(f),k(f)/pi,pb(f));
xlabel('inverse length');
ylabel('power');
line([1 1],[2e-4 25],'color',[0 0 0]);
axis([0.1   20  2e-4 25]);
legend(h,{'s','b'},'location','NE');
a=axis;text(0.12,12,'c)','horizontalalignment','left','verticalalignment','top','fontweight','bold');
set(gca,'xtick',10.^(-1:2));
set(gca,'ytick',10.^(-4:2));

w0=0.1;
w1=1.7;
w2=0.1;
w2=0.1;

x1=0.26;
x2=x1 + w0 + w1 + w2 + 0.4;
x3=x2 + w0 + w1 + w2 + 0.7;
xx=x3+w1;

set(ah,'units','inches');
set(ah(1),'position',[x1 0.36 w1 w1],'clim',[0 6]);
set(ah(2),'position',[x2 0.36 w1 w1],'clim',[0 1]);
set(ah(3),'position',[x3 0.36 w1 w1],'clim',[0 1]);
[ch1 hp hl]=jcolorbar(ah(1));
[ch2 hp hl]=jcolorbar(ah(2));
set(ch1,'units','inches');
set(ch2,'units','inches');
set(ch1,'position',[x1+w1+w0 0.36 w2 w1],'clim',[0 6]);
set(ch2,'position',[x2+w1+w2 0.36 w2 w1],'clim',[0 1]);

if ~isempty(dd)
  print('-depsc2',[dd '/dns-cluster.eps']);
end


return;
nm='pm/f7/e/31';
fn=ded_get_fn(nm,'final',[],'state');
b=ded_read_hdf(fn{end});

q.sc=[0 4];
q.bc=[0 1];
q.uc=[-0.2 2];
q.rgx=[ 4 28];
q.rgy=[-2  2];

ded_pm_f7_e2('pm/f7/e/31',b,q,'~/Dropbox/Particle-Laden-Plumes/ofigs');



ded_pm_f7_e2('pm/f7/e/27',[],q,'~/Dropbox/Particle-Laden-Plumes/ofigs');

b=ded_pm_f7_e2(nm,b,q,'~/Dropbox/NSF-OCE-Sediments/2020NSF/ofigs');


ded_pm_f7_e2('pm/f7/e/28',[],q,'~/Dropbox/Particle-Laden-Plumes/ofigs');
