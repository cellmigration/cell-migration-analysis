function [] = persistence_analysis(control_pers, hadase_pers)

figure
scatter(hadase_pers(:,1),hadase_pers(:,2),25,'LineWidth',1.5,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980]);
hold on;
scatter(control_pers(:,1),control_pers(:,2),25,'LineWidth',1.5,'MarkerEdgeColor',[0, 0.4470, 0.7410]);
hold off

title('PC3 Cells','FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Speed [um/min]','FontSize',20)
ylim([0 5])

legend('HAdase','Control')

n_bins = 18;
max_angle = 180;
x_axis_offset = max_angle/n_bins/2;

%% 1) for Hadase cells

[N, edges, bins] = histcounts(hadase_pers(:,1), n_bins);
% first plot a histogram of delta theta distribution
% histogram('BinEdges', edges, 'BinCounts', N)
% then average the speeds per bin
for n = 1:n_bins
    bin_means(:,n) = nanmean(hadase_pers(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N);

figure
hold on;

bPlotHA = bar(edges(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')

%% 2) for Control cells
[N, edges, bins] = histcounts(control_pers(:,1), n_bins);
% binlabels = categorical({
% first plot a histogram of delta theta distribution
% histogram('BinEdges', edges, 'BinCounts', N)
% then average the speeds per bin
for n = 1:n_bins
    bin_means(:,n) = nanmean(control_pers(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N);

bPlotCN = bar(edges(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0, 0.4470, 0.7410]);
errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')

xticks(edges)
title('PC3 Cells','FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Speed [um/min]','FontSize',20)
grid on;

legend([bPlotHA,bPlotCN], 'Hadase','Control')
hold off;
