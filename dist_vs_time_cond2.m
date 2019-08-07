function [coeff_ha, coeff_control] = dist_vs_time_cond2(concentration,ha_features, control_features)

% 1) HA features
time = 0:5:36*5;

ha_dist = zeros(37,2); % col1 = averages, col1 = sem

for i = 1:37
%     speeds_array = features(find(floor(features(:,3))==5*i),4);
    ha_dist_array = ha_features(find((floor(ha_features(:,3))-mod(floor(ha_features(:,3)),5))==5*(i-1)),10);
    ha_dist(i,1) = nanmean(ha_dist_array);
    ha_dist(i,2) = ha_dist(i,1)/sqrt(length(ha_dist_array)-sum(isnan(ha_dist_array)));
end
% Compute linear fit, residual standard deviation (standard error of the regression) and rsq value
[linear_fit, s] = polyfitZero(sqrt(time'),ha_dist(:,1),1);
slope = linear_fit(1); 
f_ha = polyval(linear_fit,sqrt(time'));

SSresid = sqrt(s.normr/s.df);

rsq = 1 - sum((ha_dist(:,1)-f_ha).^2)/sum((ha_dist(:,1)-mean(ha_dist(:,1))).^2);
coeff_ha = [slope, SSresid, rsq];

equationString_ha = sprintf('y = %.3g x^{1/2} + %.3g', linear_fit(1), linear_fit(2));
rsqString_ha = sprintf('S = %.3g, R^2 = %.3g', SSresid, rsq);



% 2) control features
time = 0:5:36*5;

control_dist = zeros(37,2); % col1 = averages, col1 = sem

for i = 1:37
%     speeds_array = features(find(floor(features(:,3))==5*i),4);
    control_dist_array = control_features(find((floor(control_features(:,3))-mod(floor(control_features(:,3)),5))==5*(i-1)),10);
    control_dist(i,1) = nanmean(control_dist_array);
    control_dist(i,2) = control_dist(i,1)/sqrt(length(control_dist_array)-sum(isnan(control_dist_array)));
end
% Compute linear fit, residual standard deviation (standard error of the regression) and rsq value
[linear_fit, s] = polyfitZero(sqrt(time'),control_dist(:,1),1);
slope = linear_fit(1); 
f_control = polyval(linear_fit,sqrt(time'));

SSresid = sqrt(s.normr/s.df);

rsq = 1 - sum((control_dist(:,1)-f_control).^2)/sum((control_dist(:,1)-mean(control_dist(:,1))).^2);
coeff_control = [slope, SSresid, rsq];

equationString_control = sprintf('y = %.3g x^{1/2} + %.3g', linear_fit(1), linear_fit(2));
rsqString_control = sprintf('S = %.3g, R^2 = %.3g', SSresid, rsq);

dist_vs_time_fig = figure;

ha_err_bar = errorbar(time,ha_dist(:,1),ha_dist(:,2),'.','color',[0.8500, 0.3250, 0.0980]);
hold on
control_err_bar = errorbar(time,control_dist(:,1),control_dist(:,2),'.', 'color', [0, 0.4470, 0.7410]);
ylim([0 60])

% Plot linear fit 

plot(time, f_ha,'color',[0.8500, 0.3250, 0.0980]);plot(time, f_control,'color',[0, 0.4470, 0.7410]);
text(10, 55, equationString_ha,'color',[0.8500, 0.3250, 0.0980]) %max(speeds(:,1))+0.4
text(10, 50, rsqString_ha,'color',[0.8500, 0.3250, 0.0980]) %max(speeds(:,1))+0.2
text(10, 45, equationString_control,'color',[0, 0.4470, 0.7410]) %max(speeds(:,1))+0.4
text(10, 40, rsqString_control,'color',[0, 0.4470, 0.7410]) %max(speeds(:,1))+0.2
hold off

title(['PC3 Cells, ', num2str(concentration), ' ug/ml Collagen'])
xlabel('Time [min]','FontSize',20)
ylabel('Cell Distance [um]','FontSize',20)
legend([ha_err_bar,control_err_bar], 'Hadase','Control')

figure_directory = '\\phys34212\migrationdata\MigrationData\Migration1\figures\PC3 cells\dist_vs_time\by_condition';
figure_file_name = [figure_directory, '\',num2str(concentration),'.png'];
print(dist_vs_time_fig,'-dpng',figure_file_name)