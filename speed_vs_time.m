function [] = speed_vs_time(ctrl_or_ha,concentration,features)

time = 5:5:36*5;

speeds = zeros(36,2); % col1 = averages, col1 = sem

for i = 1:36
%     speeds_array = features(find(floor(features(:,3))==5*i),4);
    speeds_array = features(find((floor(features(:,3))-mod(floor(features(:,3)),5))==5*i),4);
    speeds(i,1) = nanmean(speeds_array);
    speeds(i,2) = speeds(i,1)/sqrt(length(speeds_array)-sum(isnan(speeds_array)));
end

% figure(1)
% 
% plot(time,speeds(:,1))
% title([ctrl_or_ha,' ',num2str(concentration), ' ug/ml'])
% xlabel(

figure

errorbar(time,speeds(:,1),speeds(:,2))
title([ctrl_or_ha,', ',num2str(concentration), ' ug/ml Fibronectin'])
xlabel('Time (min)')
ylabel('Cell Speed um/min')