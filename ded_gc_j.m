figure;
[c f pp]=ded_gc_fit_g_Re('gc/f7/jb/*/*',8);

dc;

figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/jb/*/*','X','Re');
figure;[c f pp]=ded_gc_fit_Y_X('gc/f7/jb/*/*','g','Re');
figure;[c f pp]=ded_gc_fit_Z_Y_X('gc/f7/jb/*/*','Re','X','g');
h=legend;delete(h);
ded_plot_X('gc/f7/jb/*/*',2);
