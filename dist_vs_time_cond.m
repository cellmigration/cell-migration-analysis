function coeff = dist_vs_time_cond(ctrl_or_ha,concentration,features)

time = 0:5:36*5;

dist = zeros(37,2); % col1 = averages, col1 = sem

for i = 1:37
%     speeds_array = features(find(floor(features(:,3))==5*i),4);
    dist_array = features(find((floor(features(:,3))-mod(floor(features(:,3)),5))==5*(i-1)),10);
    dist(i,1) = nanmean(dist_array);
    dist(i,2) = dist(i,1)/sqrt(length(dist_array)-sum(isnan(dist_array)));
end
% Compute linear fit, residual standard deviation (standard error of the regression) and rsq value
[linear_fit, s] = polyfitZero(sqrt(time'),dist(:,1),1);
slope = linear_fit(1); 
f = polyval(linear_fit,sqrt(time'));

SSresid = sqrt(s.normr/s.df);

rsq = 1 - sum((dist(:,1)-f).^2)/sum((dist(:,1)-mean(dist(:,1))).^2);
coeff = [slope, SSresid, rsq];

equationString = sprintf('y = %.3g x^{1/2} + %.3g', linear_fit(1), linear_fit(2));
rsqString = sprintf('S = %.3g', SSresid');

dist_vs_time_fig = figure;

errorbar(time,dist(:,1),dist(:,2))
ylim([0 60])

% Plot linear fit 
hold on
plot(time, f)
text(10, 55, equationString) %max(speeds(:,1))+0.4
text(10, 45, rsqString) %max(speeds(:,1))+0.2
hold off

title([ctrl_or_ha,', ',num2str(concentration), ' ug/ml Fibronectin'])
xlabel('Time [min]','FontSize',20)
ylabel('Cell Distance [um]','FontSize',20)

figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\dist_vs_time\by_condition';
figure_file_name = [figure_directory, '\',ctrl_or_ha,'_',num2str(concentration),'.png'];
print(dist_vs_time_fig,'-dpng',figure_file_name)