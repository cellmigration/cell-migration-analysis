function coeff = angle_vs_time_cond(ctrl_or_ha,concentration,features)

time = 0:5:36*5;

angle = zeros(37,2); % col1 = averages, col1 = sem

for i = 1:37
    angle_array = features(find((floor(features(:,3))-mod(floor(features(:,3)),5))==5*(i-1)),5);
    angle(i,1) = nanmean(angle_array);
    angle(i,2) = angle(i,1)/sqrt(length(angle_array)-sum(isnan(angle_array)));
end

figure
errorbar(time,angle(:,1),angle(:,2))
title([ctrl_or_ha,', ',num2str(concentration), ' ug/ml Fibronectin'])
xlabel('Time [min]','FontSize',20)
ylabel(['Angle [', char(176) ,']'],'FontSize',20)
figure
polarhistogram(angle(:,1),12)
title([ctrl_or_ha,', ',num2str(concentration), ' ug/ml Fibronectin'])
% ylim([0 60])

% Plot linear fit 
% hold on
% plot(time, f)
% text(10, 55, equationString) %max(speeds(:,1))+0.4
% text(10, 45, rsqString) %max(speeds(:,1))+0.2
% hold off
% 
% title([ctrl_or_ha,', ',num2str(concentration), ' ug/ml Fibronectin'])
% xlabel('Time [min]','FontSize',20)
% ylabel('Cell Distance [um]','FontSize',20)
% 
% figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\dist_vs_time\by_condition';
% figure_file_name = [figure_directory, '\',ctrl_or_ha,'_',num2str(concentration),'.png'];
% print(dist_vs_time_fig,'-dpng',figure_file_name)