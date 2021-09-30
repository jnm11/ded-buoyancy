ded_pm_plot_fluxes('pm/f7/e/05',[140 inf]);
ded_pm_plot_fluxes('pm/f7/e/06',[140 inf]);
ded_pm_plot_fluxes('pm/f7/e/07',[80 inf]);


ded_pm_test_f7('pm/f7/e/05',[140 inf]);
ded_pm_test_f7('pm/f7/e/06',[140 inf]);
ded_pm_test_f7('pm/f7/e/07',[80 inf]);

ded_pm_MTT('pm/f7/e/05',[140 inf],[5 25]); %0.132
ded_pm_MTT('pm/f7/e/06',[140 inf],[5 25]); %0.131
ded_pm_MTT('pm/f7/e/07',[ 80 inf],[5 25]); %0.134

ded_pm_fit_fluxes('pm/f7/e/05',[140 inf],[5 25],1); %0.132
ded_pm_fit_fluxes('pm/f7/e/06',[140 inf],[5 25],1); %0.131
ded_pm_fit_fluxes('pm/f7/e/07',[ 80 inf],[5 25],1); %0.134

