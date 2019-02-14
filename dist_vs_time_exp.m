function slope = dist_vs_time_exp(experiment,features)

time = 0:5:36*5;

dist = zeros(37,2); % col1 = averages, col1 = sem

for i = 1:37
%     speeds_array = features(find(floor(features(:,3))==5*i),4);
    dist_array = features(find((floor(features(:,3))-mod(floor(features(:,3)),5))==5*(i-1)),10);
    dist(i,1) = nanmean(dist_array);
    dist(i,2) = dist(i,1)/sqrt(length(dist_array)-sum(isnan(dist_array)));
end

% Com;pute linear fit and rsq value
linear_fit = polyfit(time',dist(:,1),1);
slope = linear_fit(1); 
f = polyval(linear_fit,time');
rsq = 1 - sum((dist(:,1)-f).^2)/sum((dist(:,1)-mean(dist(:,1))).^2);

equationString = sprintf('y = %.3g x + %.3g', linear_fit(1), linear_fit(2));
rsqString = sprintf('R^2 = %.3g', rsq');

speed_vs_time_fig = figure;

errorbar(time,dist(:,1),dist(:,2))
% ylim([0 2])

% Plot linear fit 
hold on
plot(time, f)
% text(10, 1.9, equationString) %max(speeds(:,1))+0.4
% text(10, 1.7, rsqString) %max(speeds(:,1))+0.2
hold off

title(['Experiment ',experiment])
xlabel('Time [min]','FontSize',20)
ylabel('Cell distance [um]','FontSize',20)

% figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\speed_vs_time\by_experiment';
% figure_file_name = [figure_directory, '\',experiment,'.png'];
% print(speed_vs_time_fig,'-dpng',figure_file_name)