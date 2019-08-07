function [] = persitence_cond2(concentration, ha_features, control_features)
% control plotted on top of hadase
hadase_pers = ha_features(:,5:6);
control_pers = control_features(:,5:6);

%% Scatter plot

persitence_fig_scatter = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(hadase_pers(:,1),hadase_pers(:,2),25,'LineWidth',1.5,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980]);
hold on;
scatter(control_pers(:,1),control_pers(:,2),25,'LineWidth',1.5,'MarkerEdgeColor',[0, 0.4470, 0.7410]);
hold off

ylim([0 5])
title(['PC3 Cells ', num2str(concentration), ' ug/ml Collagen'],'FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Av. Speed [um/min]','FontSize',20)
% ylabel('|Av. Velocity| [um/min]','FontSize',20)

legend('HAdase','Control')

figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\speed_vs_delta_theta\by_condition';
figure_file_name = [figure_directory, '\','scatter_',num2str(concentration),'.png'];
print(persitence_fig_scatter,'-dpng',figure_file_name)


%% Error plot


n_bins = 18;
max_angle = 180;
x_axis_offset = max_angle/n_bins/2;

% 1) for Hadase cells

[N_hadase, edges_hadase, bins_hadase] = histcounts(hadase_pers(:,1), n_bins);

for n = 1:n_bins
    bin_means(:,n) = nanmean(hadase_pers(bins_hadase == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N_hadase);

persitence_fig_bar = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;

% bPlotHA = bar(edges_hadase(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0.8500, 0.3250, 0.0980]);
% errorbar(edges_hadase(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')
errPlotHA = errorbar(edges_hadase(1:n_bins)+x_axis_offset,bin_means,bin_sem,'*','color',[0.8500, 0.3250, 0.0980]);

% Do a x^2 fit
x_array = edges_hadase(1:n_bins)+x_axis_offset;
y_array = bin_means;

ha_fit = fit(x_array', y_array', 'poly2');
f = feval(ha_fit,x_array');
ha_fit_rsq = 1 - sum((y_array'-f).^2)/sum((y_array'-mean(y_array')).^2);
ha_fit_equation_string = sprintf('y = %.3g x^2 + %.3g x + %.3g', ha_fit.p1, ha_fit.p2, ha_fit.p3);
ha_fit_rsq_string = sprintf('R^2 = %.3g', ha_fit_rsq);
ha_fit_string = [ha_fit_equation_string, ', ',ha_fit_rsq_string];

ha_fit_plot = plot(ha_fit);
set(ha_fit_plot, 'LineWidth',1.5,'color',[0.8500, 0.3250, 0.0980]) %'LineStyle','--'


% 2) for Control cells
[N_control, edges_control, bins] = histcounts(control_pers(:,1), n_bins);

for n = 1:n_bins
    bin_means(:,n) = nanmean(control_pers(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N_control);

% bPlotCN = bar(edges_control(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0, 0.4470, 0.7410]);
% errorbar(edges_control(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')
errPlotCN = errorbar(edges_control(1:n_bins)+x_axis_offset,bin_means,bin_sem,'*','color',[0, 0.4470, 0.7410]);

% Do a x^2 fit
x_array = edges_control(1:n_bins)+x_axis_offset;
y_array = bin_means;

control_fit = fit(x_array', y_array', 'poly2');
f = feval(control_fit,x_array');
control_fit_rsq = 1 - sum((y_array'-f).^2)/sum((y_array'-mean(y_array')).^2);
control_fit_equation_string = sprintf('y = %.3g x^2 + %.3g x + %.3g', control_fit.p1, control_fit.p2, control_fit.p3);
control_fit_rsq_string = sprintf('R^2 = %.3g', control_fit_rsq);
control_fit_string = [control_fit_equation_string, ', ',control_fit_rsq_string];

control_fit_plot = plot(control_fit);
% set(control_fit_plot(2), 'LineWidth',1.5,'color','k')
set(control_fit_plot, 'LineWidth',1.5,'color',[0, 0.4470, 0.7410])

% display fit equations
plot_ylim = ylim;
plot_miny = plot_ylim(1);
plot_maxy = plot_ylim(2);
text(10, plot_miny + 0.9*(plot_maxy-plot_miny), ha_fit_string,'color',[0.8500, 0.3250, 0.0980])
text(10, plot_miny + 0.8*(plot_maxy-plot_miny), control_fit_string,'color',[0, 0.4470, 0.7410])

xticks(edges_control)
title(['PC3 Cells ', num2str(concentration), ' ug/ml Collagen'],'FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Av. Speed [um/min]','FontSize',20)
% ylabel('|Av. Velocity| [um/min]','FontSize',20)
grid on;

legend([errPlotHA,ha_fit_plot, errPlotCN, control_fit_plot], 'Hadase','Hadase fit', 'Control', 'Control fit')
hold off;

figure_file_name = [figure_directory, '\','err_',num2str(concentration),'.png'];
print(persitence_fig_bar,'-dpng',figure_file_name)

%% Histograms of delta theta

persitence_fig_hist = figure('units','normalized','outerposition',[0 0 1 1]);

histogram('BinEdges', edges_control, 'BinCounts', N_control, 'Normalization','probability')
hold on;
histogram('BinEdges', edges_hadase, 'BinCounts', N_hadase, 'Normalization','probability')

xticks(edges_control)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
title(['PC3 Cells ', num2str(concentration), ' ug/ml Collagen'],'FontSize',20)

legend('Control','Hadase')

figure_file_name = [figure_directory, '\','hist_',num2str(concentration),'.png'];
print(persitence_fig_hist,'-dpng',figure_file_name)

