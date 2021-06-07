% Quick script to read and plot all individual lcs on one lc (as structured in
% year1/shape_modelling/162p/162P_14_artificial_lc_data_plt.dat)

%% Read in lc data
file_dir = "~/Documents/year1/shape_modelling/162p/";
fid = fopen(file_dir+"162P_14_artificial_lc_data_plt.dat", 'r');

sizecols = [ 4 Inf ]
formatSpec = '%i %f %f %f';

lc_dat = fscanf(fid,formatSpec,sizecols)';

    lcno = lc_dat(:,1);
   phase = lc_dat(:,2);
 rel_mag = lc_dat(:,3);
    uncs = lc_dat(:,4);


 fontsize=18;
 %scatter(phase, rel_mag, 30, lcno, 'filled', 'd')
 errorbar(phase, rel_mag, uncs, 'o')
 colormap winter;
 line(xlim(), [0,0], 'LineWidth', 1, 'LineStyle', ':', 'Color', [17 17 17]/255);
 legend({'2007-2018'}, 'Interpreter', 'latex', 'Fontsize', 20, 'Location', 'southeast')