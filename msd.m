function [alpha,alpha_sem, rsq] = msd(features)

% initialize array for ensemble analysis
% tau_msd_array = [];

% loop over cells
cell_indices = find(isnan(features(:,4)));
cell_total = length(cell_indices);
cells_used = 0;

alphas = [];
rsqs = [];

for i = 1:length(cell_indices)-1 % fix this to include the last cell
    
   cell_features = features(cell_indices(i):cell_indices(i+1)-1,:);
   
   duration_cutoff = 15;
   if size(cell_features,1) > duration_cutoff
       
       cells_used = cells_used + 1;
     
       x = cell_features(:,1);
       y = cell_features(:,2);
       time = cell_features(:,3);

       timepts = size(cell_features,1);

    %    time = 0:5:5*timepts-1;
       msd = zeros(timepts-1,1);

       for tau = 1:timepts-1
           for tau2 = tau:timepts-1
               msd(tau) = msd(tau)+ (x(tau2+1)-x(tau2+1-tau))^2 + (y(tau2+1)-y(tau2+1-tau))^2;
           end
           msd(tau) = msd(tau)/(timepts-tau+1);
       end

%        close all
%        figure(51)
%        plot(x,y)
%        xlabel('x-coordinate (um)')
%        ylabel('y-coordinate (um)')
%        xlim([-50 50])
%        ylim([-50 50])
%        figure(52)
%        plot(time(2:end),msd)
% 
%        ylabel('MSD')
%        xlabel('Tau')
%        figure(53)
% 
%        plot(log(time(2:end)),log(msd))
%        ylabel('Ln(MSD)')
%        xlabel('Ln(Tau)')

       linear_fit = polyfit(log(time(2:end)),log(msd),1);
       alpha_cell = linear_fit(1);
       f = polyval(linear_fit,log(time(2:end)));
       rsq_cell = 1 - sum((log(msd)-f).^2)/sum((log(msd)-mean(log(msd))).^2);

       alphas = [alphas; alpha_cell];
       rsqs = [rsqs; rsq_cell];
       
       % update array for ensemble analysis
%        tau_msd_array = [tau_msd_array;time(2:end),msd];

%        equationString = sprintf('y = %.3g x + %.3g', linear_fit(1), linear_fit(2));
%        rsqString = sprintf('R^2 = %.3g', rsq_cell);
% 
%        % Plot linear fit 
%        hold on
%        plot(log(time(2:end)), f)
%        text(1, max(log(msd))-1, equationString)
%        text(1, max(log(msd))-1.5, rsqString)
%        hold off   

%        display('check')
   end
       
end

alpha = nanmean(alphas);
alpha_sem = nanstd(alphas)/sqrt(cells_used);
rsq = mean(rsqs);

% % perform ensemble statistics
% tau_length = 72;
% msd_ensemble = zeros(tau_length,2); % col1 = averages, col1 = sem
% 
% for i = 1:tau_length
% %     speeds_array = features(find(floor(features(:,3))==5*i),4);
%     msd_array = tau_msd_array(find((floor(tau_msd_array(:,1))-mod(floor(tau_msd_array(:,1)),2.5))==2.5*(i-1)),2);
%     msd_ensemble(i,1) = nanmean(msd_array);
%     msd_ensemble(i,2) = msd_ensemble(i,1)/sqrt(length(msd_array)-sum(isnan(msd_array)));
% end
% 
% dtau = 2.5; 
% tau_array = dtau:dtau:tau_length*dtau;
% figure
% errorbar(tau_array,msd_ensemble(:,1),msd_ensemble(:,2));
% ylabel('MSD')
% xlabel('Tau')
% 
% figure
% plot(log(tau_array),log(msd_ensemble(:,1)))
% ylabel('Ln(MSD)')
% xlabel('Ln(Tau)')
% 
% linear_fit = polyfit(log(tau_array)',log(msd_ensemble(:,1)),1);
% alpha_condition = linear_fit(1);
% f = polyval(linear_fit,log(tau_array)');
% rsq_condition = 1 - sum((log(msd_ensemble(:,1))-f).^2)/sum((log(msd_ensemble(:,1))-mean(log(msd_ensemble(:,1)))).^2);
% 
% equationString = sprintf('y = %.3g x + %.3g', linear_fit(1), linear_fit(2));
% rsqString = sprintf('R^2 = %.3g', rsq_condition);
% 
% % Plot linear fit 
% hold on
% plot(log(tau_array), f)
% text(1, max(log(msd_ensemble(:,1)))-1, equationString)
% text(1, max(log(msd_ensemble(:,1)))-1.5, rsqString)
% hold off   
% 
% disp('check')