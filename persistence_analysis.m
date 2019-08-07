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

[N, edges_hadase, bins] = histcounts(hadase_pers(:,1), n_bins);
% first plot a histogram of delta theta distribution
% histogram('BinEdges', edges, 'BinCounts', N)
% then average the speeds per bin
for n = 1:n_bins
    bin_means(:,n) = nanmean(hadase_pers(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N);

figure
hold on;

% bPlotHA = bar(edges(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0.8500, 0.3250, 0.0980]);
% errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')

errPlotHA = errorbar(edges_hadase(1:n_bins)+x_axis_offset,bin_means,bin_sem,'color',[0.8500, 0.3250, 0.0980]); %'*',

% Do a x^2 fit
% x_array = edges_hadase(1:n_bins)+x_axis_offset;
% y_array = bin_means;
% use the raw data instead of the binned data
x_array = rmmissing(hadase_pers(:,1)');
y_array = rmmissing(hadase_pers(:,2)');

ha_fit = fit(x_array', y_array', 'poly2');
f = feval(ha_fit,x_array');
ha_fit_rsq = 1 - sum((y_array'-f).^2)/sum((y_array'-mean(y_array')).^2);
ha_fit_equation_string = sprintf('y = %.3g x^2 + %.3g x + %.3g', ha_fit.p1, ha_fit.p2, ha_fit.p3);
ha_fit_rsq_string = sprintf('R^2 = %.3g', ha_fit_rsq);
ha_fit_string = [ha_fit_equation_string, ', ',ha_fit_rsq_string];

% Do a x^2 fit with the whole data

ha_fit_plot = plot(ha_fit);
set(ha_fit_plot, 'LineWidth',1.5,'color',[0.8500, 0.3250, 0.0980]) %'LineStyle','--'

%% 2) for Control cells
[N, edges_control, bins] = histcounts(control_pers(:,1), n_bins);
% binlabels = categorical({
% first plot a histogram of delta theta distribution
% histogram('BinEdges', edges, 'BinCounts', N)
% then average the speeds per bin
for n = 1:n_bins
    bin_means(:,n) = nanmean(control_pers(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N);

% bPlotCN = bar(edges(1:n_bins)+x_axis_offset,bin_means,'BarWidth',1,'FaceColor',[0, 0.4470, 0.7410]);
% errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'k.')

errPlotCN = errorbar(edges_control(1:n_bins)+x_axis_offset,bin_means,bin_sem,'color',[0, 0.4470, 0.7410]);%'*'

% Do a x^2 fit
% x_array = edges_control(1:n_bins)+x_axis_offset;
% y_array = bin_means;
% use the raw data instead of the binned data
x_array = rmmissing(control_pers(:,1)');
y_array = rmmissing(control_pers(:,2)');


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
text(10, plot_miny + 0.95*(plot_maxy-plot_miny), ha_fit_string,'color',[0.8500, 0.3250, 0.0980])
text(10, plot_miny + 0.85*(plot_maxy-plot_miny), control_fit_string,'color',[0, 0.4470, 0.7410])


xticks(edges_control)
title('PC3 Cells - All Conditions','FontSize',20)
xlabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
ylabel('Av. Speed [um/min]','FontSize',20)
grid on;

legend([errPlotHA,ha_fit_plot, errPlotCN, control_fit_plot], 'Hadase','Hadase fit', 'Control', 'Control fit')
hold off;
