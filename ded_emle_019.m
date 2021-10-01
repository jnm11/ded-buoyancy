nm='gc/emle/019';
dd='~/Dropbox/Jim-Claudia-GC/mat/emle/019';
%ded_gc_ccle_process
p=ded_read_param(nm);
a=ded_gc_traj(nm,dd,[5 33.5],1);
ded_gc_steady(nm,dd,[16 33.5],a.fx,a.fu, 'y',[-20 5],'16');
ded_gc_steady(nm,dd,[20 33.5],a.fx,a.fu, 'y',[-20 5],'20');
ded_gc_steady(nm,dd,[24 33.5],a.fx,a.fu, 'y',[-20 5],'24');
ded_gc_steady(nm,dd,[28 33.5],a.fx,a.fu, 'y',[-20 5],'28');
ded_gc_steady(nm,dd,[32 33.5],a.fx,a.fu, 'y',[-20 5],'32');
ded_gc_yint(nm,dd);
ded_emlecmp
ded_gc_ccle_cplots({'emle/019','emle/017','ccle/024','ccle/046'},33.5,'017-33-')
ded_gc_slice('gc/emle/019');

!rsync -uavp ~/Dropbox/Jim-Claudia-GC asahi:Dropbox


fd='~/Dropbox/Jim-Claudia-GC/ofigs';
preprint([5 1.3],6);colour_lines;clf;
ded_gc_cmp_b_cnts({'emle/019/16','emle/017/16','ccle/024/18','ccle/046/16'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-16.eps']);

ded_gc_cmp_b_cnts({'emle/019/20','emle/017/20','ccle/024/22','ccle/046/20'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-20.eps']);
                
ded_gc_cmp_b_cnts({'emle/019/24','emle/017/24','ccle/024/26','ccle/046/24'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-24.eps']);

ded_gc_cmp_b_cnts({'emle/017/28','ccle/024/30','ccle/046/28'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-28.eps']);

ded_gc_cmp_b_cnts({'emle/017/32'},[0.02 0.5]);
print('-depsc2','-loose',[fd '/bcnt-32.eps']);




