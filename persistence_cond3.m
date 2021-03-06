function [] = persistence_cond3(concentration, ha_features, control_features)
% hadase plotted on top of control
hadase_pers = ha_features(:,5:6);
control_pers = control_features(:,5:6);

%% Scatter plot

persitence_fig_scatter = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(hadase_pers(:,1),hadase_pers(:,2),25,[0.8500, 0.3250, 0.0980]);
hold on;
scatter(control_pers(:,1),control_pers(:,2),25,[0, 0.4470, 0.7410]);
hold off

title(['PC3 Cells ', num2str(concentration), ' ug/ml Collagen'],'FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
% ylabel('Av. Speed [um/min]','FontSize',20)
ylabel('|Av. Velocity| [um/min]','FontSize',20)

legend('HAdase','Control')

figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\speed_vs_delta_theta\by_condition';
figure_file_name = [figure_directory, '\','scatter_',num2str(concentration),'.png'];
print(persitence_fig_scatter,'-dpng',figure_file_name)


%% Bar plot


n_bins = 18;
max_angle = 180;
x_axis_offset = max_angle/n_bins/2;

% 1) for control cells

[N, edges, bins] = histcounts(control_pers(:,1), n_bins);

for n = 1:n_bins
    bin_means(:,n) = nanmean(control_pers(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N);

persitence_fig_bar = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;

bPlotCN = bar(edges(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0, 0.4470, 0.7410]);
errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')

% 2) for Hadase cells
[N, edges, bins] = histcounts(hadase_pers(:,1), n_bins);

for n = 1:n_bins
    bin_means(:,n) = nanmean(hadase_pers(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N);

bPlotHA = bar(edges(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')

xticks(edges)
title(['PC3 Cells ', num2str(concentration), ' ug/ml Collagen'],'FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Av. Speed [um/min]','FontSize',20)
% ylabel('|Av. Velocity| [um/min]','FontSize',20)
grid on;

legend([bPlotCN,bPlotHA], 'Control','Hadase')
hold off;

figure_file_name = [figure_directory, '\','bar_',num2str(concentration),'.png'];
print(persitence_fig_bar,'-dpng',figure_file_name)
