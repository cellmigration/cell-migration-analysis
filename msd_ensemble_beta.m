function [alpha_condition, SSresid_condition, rsq_condition] = msd_ensemble_beta(cell_type, substrate_type, treat, conc, features, time_threshold)

% initialize array for ensemble analysis
tau_msd_array = [];

% loop over cells
cell_indices = find(isnan(features(:,4)));
cell_total = length(cell_indices);
cells_used = 0;

for i = 1:length(cell_indices)-1 % fix this to include the last cell
    
   cell_features = features(cell_indices(i):cell_indices(i+1)-1,:);
   % Implement time threshold
   cell_features = cell_features(cell_features(:,3)<= time_threshold, :);

   duration_cutoff = 0;
   if size(cell_features,1) > duration_cutoff
       
       cells_used = cells_used + 1;
     
       x = cell_features(:,1);
       y = cell_features(:,2);
       time = cell_features(:,3);

       timepts = size(cell_features,1);

       msd = zeros(timepts-1,1);

       for tau = 1:timepts-1
           for tau2 = tau:timepts-1
               msd(tau) = msd(tau)+ (x(tau2+1)-x(tau2+1-tau))^2 + (y(tau2+1)-y(tau2+1-tau))^2;
           end
           msd(tau) = msd(tau)/(timepts-tau+1);
       end
       
       % update array for ensemble analysis
       tau_indices = 1:timepts-1;
       tau_msd_array = [tau_msd_array;tau_indices',msd];
       
   end
       
end

% perform ensemble statistics
tau_length = 72;
msd_ensemble = zeros(tau_length,3); % col1 = averages, col1 = sem

for i = 1:tau_length
%     speeds_array = features(find(floor(features(:,3))==5*i),4);
%     msd_array = tau_msd_array(find((floor(tau_msd_array(:,1))-mod(floor(tau_msd_array(:,1)),2.5))==2.5*(i-1)),2);
    msd_array = tau_msd_array(find(tau_msd_array(:,1)== i),2);
    msd_ensemble(i,1) = nanmean(msd_array);
    msd_ensemble(i,2) = msd_ensemble(i,1)/sqrt(size(msd_array,1)-sum(isnan(msd_array)));
    msd_ensemble(i,3) = size(msd_array,1)-sum(isnan(msd_array));
end

dtau = 2.5; 
tau_array = dtau:dtau:tau_length*dtau;

fraction_to_plot = 0.5;
index_f = floor(tau_length*fraction_to_plot);
tau_array = tau_array(1:index_f);
msd_ensemble = msd_ensemble(1:index_f,:);
data_points_min = msd_ensemble(index_f,3);

msd_figure = figure;
errorbar(tau_array,msd_ensemble(:,1),msd_ensemble(:,2));
ylabel('MSD')
xlabel('Tau')
title([cell_type,' Cells, ', treat, ' ', num2str(conc), ' ug/ml', substrate_type, ' (> ',num2str(data_points_min), ' points)'])

figure
plot(tau_array, msd_ensemble(:,3))
xlabel('Tau')
ylabel('Data points used')

alpha_figure = figure;
hold on
plot(log(tau_array),log(msd_ensemble(:,1)))
ylabel('Ln(MSD)')
xlabel('Ln(Tau)')
title([cell_type,' Cells, ', treat, ' ', num2str(conc), ' ug/ml', substrate_type])

[linear_fit, s] = polyfit(log(tau_array)',log(msd_ensemble(:,1)),1);
alpha_condition = linear_fit(1);
f = polyval(linear_fit,log(tau_array)');
rsq_condition = 1 - sum((log(msd_ensemble(:,1))-f).^2)/sum((log(msd_ensemble(:,1))-mean(log(msd_ensemble(:,1)))).^2);
SSresid_condition = sqrt(s.normr/s.df);

equationString = sprintf('y = %.3g x + %.3g', linear_fit(1), linear_fit(2));
rsqString = sprintf('R^2 = %.3g', rsq_condition);

% Plot linear fit 
legend_label = [equationString, ', ', rsqString];
plot(log(tau_array), f)
legend(legend_label)
hold off

% 
% figure_directory = '\\phys34212\migrationdata\MigrationData\Migration1\figures\MSD';
% msd_figure_file_name = [figure_directory, '\','msd-',treat,'-', num2str(conc),'.png'];
% print(msd_figure,'-dpng',msd_figure_file_name)
% alpha_figure_file_name = [figure_directory, '\','alpha-',treat,'-', num2str(conc),'.png'];
% print(alpha_figure,'-dpng',alpha_figure_file_name)