function [] = persistence_cond_beta(treatment, features, fig_scatter, fig_bar, fig_hist)
% Grab delta theta and average speed columns
persistence_features = features(:,5:6);

%% Scatter plot

figure(fig_scatter)

legend_label = treatment;
scatter(persistence_features(:,1),persistence_features(:,2),25,'LineWidth',1.5, 'DisplayName',legend_label);


%% Error plot

figure(fig_bar)

n_bins = 18;
max_angle = 180;
x_axis_offset = max_angle/n_bins/2;

[N, edges, bins] = histcounts(persistence_features(:,1), n_bins);

for n = 1:n_bins
    bin_means(:,n) = nanmean(persistence_features(bins == n,2))'; % average cell speed in each bin/subinterval
end

bin_sem = bin_means./sqrt(N);

legend_label = treatment;
errPlot = errorbar(edges(1:n_bins)+x_axis_offset,bin_means,bin_sem,'*', 'DisplayName',legend_label);

% Do a x^2 fit
x_array = edges(1:n_bins)+x_axis_offset;
y_array = bin_means;

treat_fit = fit(x_array', y_array', 'poly2');
f = feval(treat_fit,x_array');
treat_fit_rsq = 1 - sum((y_array'-f).^2)/sum((y_array'-mean(y_array')).^2);
treat_fit_equation_string = sprintf('y = %.3g x^2 + %.3g x + %.3g', treat_fit.p1, treat_fit.p2, treat_fit.p3);
treat_fit_rsq_string = sprintf('R^2 = %.3g', treat_fit_rsq);
treat_fit_string = [treat_fit_equation_string, ', ',treat_fit_rsq_string];

legend_label = treat_fit_string;
fit_plot = plot(treat_fit);
set(fit_plot, 'color', errPlot.Color, 'LineWidth',1.5, 'DisplayName', legend_label)

xticks(edges) % is this the problem?

%% Histograms of delta theta

figure(fig_hist)

legend_label = treatment;
histogram('BinEdges', edges, 'BinCounts', N, 'Normalization', 'probability', 'DisplayName',legend_label)
xticks(edges)