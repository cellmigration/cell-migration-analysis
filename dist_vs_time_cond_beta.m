function [time, dist, equationString, rsqString, slope, SSresid, rsq] = dist_vs_time_cond_beta(cell_type, treatment, concentration,features)

% Extract tmin, dt, tmax from features
time_column = features(:,3);
dt = round(mode(diff(time_column)));
% diff_time_column = diff(time_column);
% dt = mean(diff_time_column(diff_time_column > 0));
tmin = round(min(time_column));
tmax = round(max(time_column));
time = (tmin:dt:tmax)';

n_timepts = numel(time);

dist = zeros(n_timepts,2); % col1 = averages, col1 = sem

for i = 1:n_timepts

    dist_array = features(time_column < (dt*i+dt/2) & time_column >= (dt*i- dt/2), 10);

    dist(i,1) = nanmean(dist_array);
    dist(i,2) = dist(i,1)/sqrt(length(dist_array)-sum(isnan(dist_array)));
end
% Compute linear fit, residual standard deviation (standard error of the regression) and rsq value
[linear_fit, s] = polyfitZero(sqrt(time),dist(:,1),1);
slope = linear_fit(1); 
f = polyval(linear_fit,sqrt(time));

SSresid = sqrt(s.normr/s.df);

rsq = 1 - sum((dist(:,1)-f).^2)/sum((dist(:,1)-mean(dist(:,1))).^2);
% coeff = [slope, SSresid, rsq];

equationString = sprintf('y = %.3g x^{1/2} + %.3g', linear_fit(1), linear_fit(2));
rsqString = sprintf('S = %.3g, R^2 = %.3g', SSresid, rsq);

% dist_vs_time_fig = figure;

legend_label = treatment;
err_bar = errorbar(time,dist(:,1),dist(:,2),'.','DisplayName',legend_label);

% ylim([0 60])

% Plot linear fit 

legend_label = [equationString,', ', rsqString];
plot(time, f, 'color', err_bar.Color, 'DisplayName', legend_label);

figure_directory = '\\phys34212\migrationdata\MigrationData\Migration1\figures\PC3 cells\dist_vs_time\by_condition';
figure_file_name = [figure_directory, '\',num2str(concentration),'.png'];
% Optional: save figures automatically
% print(dist_vs_time_fig,'-dpng',figure_file_name)