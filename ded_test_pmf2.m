




rm -rf $DEDALUS_DATA/pm/qpmf2 ; mpiexec -n 4 ded_gc.py --preset qpmf2 pm/qpmf2


a=ded_pm_plot_forcing('pm/qpmf2');


a=ded_read_hdf([ded_dedalus_data_dir 'pm/qpmf2/force/force_s1.hdf5']);
