function ded_plot_fx(nm)

c=ded_coord(nm);
p=ded_read_param(nm);
fx=ded_read_hdf(['~/' nm '/force/fx.hdf5']);

fn=intersect(fieldnames(fx),{'dudxrn','pmif','Bwxx','wmult','wbl','wbr','sbl','sbr','wdl','sdr','wdr','wul','wur','wux1','wux2','wux3','wux4','wux5','wux6','wux7','wux8','Iux1','Iux2','Iux3','Iux4','Iux5','sminf','bminf','wwl','wwr','ful','fur','fulv','furv','wdivx','wuli','wuri'});


fn=intersect(fieldnames(fx),{'wbl','sbl','wdl','sdl','wul'});
n=length(fn);