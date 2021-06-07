% Quick script to read orbital elements (as structured in
% year1/162p_orbital_els.py)
% Also calls function to calculate sub-earth latitude
% And make plot
%%
file_dir = "~/Documents/year1/";
fid = fopen(file_dir+"162p_orbital_els.txt", 'r');

sizeels = [ 7 Inf ]
formatSpec = '%f %f %f %f %f %f %f';
%%
els = fscanf(fid,formatSpec,sizeels)';

  epochs = els(:,1);
asc_node = els(:,2);
   omega = els(:,3);
     inc = els(:,4);
   phase = els(:,5);
      ra = els(:,6);
     dec = els(:,7);
     

%% Rotation pole (ecliptic coordinates):
lam_p = 300;
bet_p = 70;

lats = [ 1 size(epochs) ];

for i=1:size(epochs)
    lats(i) = AD_coords(lam_p, bet_p, asc_node(i), omega(i), inc(i), ra(i), dec(i));
end

% Data for the 1-1-2022: 
ra_jan = 192.7969;
dec_jan = 17.7512;
phase_jan = 16.2409;

inc_jan = 22.94927;
asc_node_jan = 4.273909;
omega_jan = 10.360803;

lat_jan = AD_coords(lam_p, bet_p, asc_node_jan, omega_jan, inc_jan, ra_jan, dec_jan)

phase(19) = phase_jan
lats(19) = lat_jan

%,'MarkerFaceColor', '[0.4 1 0.9]'
fontsize=18;
scatter(phase(1:18), lats(1:18), 300,'filled','x', 'MarkerEdgeColor', '[0.9 0 0.5]', 'Linewidth', 2.25)
hold on 
scatter(phase(19),lats(19),350,'d','MarkerFaceColor','[0.9 0 0.5]', 'MarkerEdgeColor', 'black', 'Linewidth', 0.75)
line(xlim(), [0,0], 'LineWidth', 1, 'LineStyle', ':', 'Color', [17 17 17]/255);
%grid on;
%title('Observing geometry: 162P/Siding-Spring','FontSize',18, 'FontWeight', 'Normal')
set(0,'defaulttextinterpreter','latex')
%set(gca,'XTickLabel', {'0', '2', '4', '6', '8', '10', '12', '14'}, 'TickLabelInterpreter', 'latex', 'Fontsize', fontsize-2);
%set(gca,'YTickLabel', {'-4', '-2', '0', '2', '4', '6', '8', '10', '12', '14'}, 'TickLabelInterpreter', 'latex', 'Fontsize', fontsize-2);
legend({'2007-2018', 'LT 2021B', '162P equator'}, 'Interpreter', 'latex', 'Fontsize', 20, 'Location', 'southeast')
%legend({'2007-2018 LCs', '162P equator'}, 'Interpreter', 'latex', 'Fontsize', fontsize-2, 'Location', 'northwest')
xlabel('Solar phase angle $\alpha$ (deg)', 'FontSize', fontsize) 
ylabel('sub-Earth latitude $\phi$ (deg)', 'Fontsize', fontsize)
hold off

%saveas(gcf, '162P_lat_phase', 'png')