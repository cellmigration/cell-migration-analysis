function [] = persitence_cond(ctrl_or_ha,concentration,features)

persitence_fig = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(features(:,5),features(:,6))

title(['PC3 Cells ', ctrl_or_ha,', ',num2str(concentration), ' ug/ml Collagen'],'FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Speed [um/min]','FontSize',20)

figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\speed_vs_delta_theta\by_condition';
figure_file_name = [figure_directory, '\','scatter_',ctrl_or_ha,'_',num2str(concentration),'.png'];
print(persitence_fig,'-dpng',figure_file_name)

 
n_bins = 15;
x_axis_offset = 360/n_bins/2;

[N, edges, bins] = histcounts(features(:,5), n_bins);
 
for n = 1:n_bins
    bin_means(:,n) = nanmean(features(bins == n,6))'; % average cell speed in each bin/subinterval
end
 
bin_sem = bin_means./sqrt(N);

persitence_fig = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
bar(edges(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1);
errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')
hold off
 
xticks(edges)

title(['PC3 Cells ', ctrl_or_ha,', ',num2str(concentration), ' ug/ml Collagen'],'FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Speed [um/min]','FontSize',20)
grid on;

figure_file_name = [figure_directory, '\','bar_',ctrl_or_ha,'_',num2str(concentration),'.png'];
print(persitence_fig,'-dpng',figure_file_name)