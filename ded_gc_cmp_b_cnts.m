function ded_gc_cmp_b_cnts(nms,c)

col=get(0,'defaultaxescolororder');
clf;hold('on');
h=zeros(length(nms),1);
for j=1:length(nms)
  nm=['~/Dropbox/Jim-Claudia-GC/mat/' nms{j} '-steady.mat'];
  load(nm,'a');
  for k=1:length(c)
    [X Y]=longest_contours(a.x,a.z,a.b,c(k),1);
    h(j)=plot(X,Y,'color',col(j,:));
  end
end
lh=legend(h,nms,'location','NE');
ah=gca;
aa=[-6 0.5 0 1];
axis(aa);
ar=(aa(4)-aa(3))/(aa(2)-aa(1));
pp=get(gcf,'paperposition');
par=pp(4)/pp(3);
set(gca,'dataaspect',[1 1 1],'box','on');
make_tight_axes;
set(lh,'position',[0.777    0.51    0.215    0.48]);
set(ah,'position',[0.049    0.18    0.95    ar/par*0.95]);

return

fd='~/Dropbox/Jim-Claudia-GC/ofigs';
preprint([5 1.3],6);colour_lines;
ded_gc_cmp_b_cnts({'emle/019/16','emle/017/16','ccle/024/18','ccle/046/16'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-16.eps']);

ded_gc_cmp_b_cnts({'emle/019/20','emle/017/20','ccle/024/22','ccle/046/20'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-20.eps']);
                
ded_gc_cmp_b_cnts({'emle/019/24','emle/017/24','ccle/024/26','ccle/046/24'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-24.eps']);

ded_gc_cmp_b_cnts({'emle/019/28','emle/017/28','ccle/024/30','ccle/046/28'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-28.eps']);

ded_gc_cmp_b_cnts({'emle/017/32'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-32.eps']);
