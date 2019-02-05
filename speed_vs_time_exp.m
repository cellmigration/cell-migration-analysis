function slope = speed_vs_time_exp(experiment,features)

time = 5:5:36*5;

speeds = zeros(36,2); % col1 = averages, col1 = sem

for i = 1:36
%     speeds_array = features(find(floor(features(:,3))==5*i),4);
    speeds_array = features(find((floor(features(:,3))-mod(floor(features(:,3)),5))==5*i),4);
    speeds(i,1) = nanmean(speeds_array);
    speeds(i,2) = speeds(i,1)/sqrt(length(speeds_array)-sum(isnan(speeds_array)));
end

% Com;pute linear fit and rsq value
linear_fit = polyfit(time',speeds(:,1),1);
slope = linear_fit(1); 
f = polyval(linear_fit,time');
rsq = 1 - sum((speeds(:,1)-f).^2)/sum((speeds(:,1)-mean(speeds(:,1))).^2);

equationString = sprintf('y = %.3g x + %.3g', linear_fit(1), linear_fit(2));
rsqString = sprintf('R^2 = %.3g', rsq');

speed_vs_time_fig = figure;

errorbar(time,speeds(:,1),speeds(:,2))
ylim([0 2])

% Plot linear fit 
hold on
plot(time, f)
text(10, 1.9, equationString) %max(speeds(:,1))+0.4
text(10, 1.7, rsqString) %max(speeds(:,1))+0.2
hold off

title(['Experiment ',experiment])
xlabel('Time (min)','FontSize',20)
ylabel('Cell Speed um/min','FontSize',20)

figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\speed_vs_time\by_experiment';
figure_file_name = [figure_directory, '\',experiment,'.png'];
% print(speed_vs_time_fig,'-dpng',figure_file_name)