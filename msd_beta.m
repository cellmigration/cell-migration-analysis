function [alpha,alpha_sem, rsq] = msd_beta(features, time_threshold)

% loop over cells
cell_indices = find(isnan(features(:,4)));
cell_total = length(cell_indices);
cells_used = 0;

alphas = [];
rsqs = [];

for i = 1:length(cell_indices)-1 % fix this to include the last cell
    
   cell_features = features(cell_indices(i):cell_indices(i+1)-1,:);
   % implement time threshold
   cell_features = cell_features(cell_features(:,3)<= time_threshold, :);
   
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
       
   end
       
end

alpha = nanmean(alphas);
alpha_sem = nanstd(alphas)/sqrt(cells_used);
rsq = mean(rsqs);
