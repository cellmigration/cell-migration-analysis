%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: this script plots speed vs. concentration for a given
% experiment. Concentrations inferred from the data rather than hard-coded.
% 
% Date: 9/23/18 PYA
%
% Update: refactoring this code to make it more user-friendly and
% incorporate in the main code
%
% Date: 1/22/2019 PYA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = visualizeFeatures(output, analysisplan)

cd(output)

[num,~,raw] = xlsread(analysisplan,'experiments');
filenums = find(num(:,7)==0)+1;

%initialize features variable
all_features = [];

for curfilenum = 1:length(filenums)
    curfile = [num2str(raw{filenums(curfilenum),1}) '_features'];
    load(curfile);
    all_features = [all_features;features];
end

features = all_features;

% Calculate total number of cells

unique_cell_ind = find(isnan(features(:,4)));
cells_total = length(unique_cell_ind);

% Migration ratio calculation

threshold = 0.2;

% Extract concentrations from data
concentrations = unique(features(:,7));
n_conc = length(concentrations);

ctr_speed = zeros(n_conc,4); % col1 = av, col2 = sem, col3 = total cells, col4 = migration ratio
ha_speed = zeros(n_conc,4); % col1 = av, col2 = sem, col3 = total cells, col4 = migration ratio

for i = 1:n_conc
    
    % Calculations for control cells
    ctrl_tmp_ind = features(:,7)==concentrations(i) & features(:,8) == 1;
    ctrl_tmp_val = features(ctrl_tmp_ind,4);
    ctr_speed(i,1) = nanmean(ctrl_tmp_val);
    sample_size = length(ctrl_tmp_val);
    ctr_speed(i,2)= nanstd(ctrl_tmp_val)/sqrt(sample_size);
    
    % compute total number of cells
    unique_cell_ind = find(isnan(ctrl_tmp_val));
    cond_cells_total = length(unique_cell_ind);
    ctr_speed(i,3) = cond_cells_total;
    
    % migration calculations
    initial_condition_rows = length(find(ctrl_tmp_ind > 0))-cond_cells_total;
    migrating_rows = find(ctrl_tmp_val > threshold);
    migration_ratio = length(migrating_rows)/initial_condition_rows*100;
    ctr_speed(i,4) = migration_ratio;
    
    % Calculations for HAS
    ha_tmp_ind = features(:,7) == concentrations(i) & features(:,8) == 2;
    ha_tmp_val = features(ha_tmp_ind,4);
    ha_speed(i,1) = nanmean(ha_tmp_val);
    sample_size = length(ha_tmp_val);
    ha_speed(i,2)= nanstd(ha_tmp_val)/sqrt(sample_size);
    
    % compute total number of cells
    unique_cell_ind = find(isnan(ha_tmp_val));
    cond_cells_total = length(unique_cell_ind);
    ha_speed(i,3) = cond_cells_total;
    
    % migration calculations
    initial_condition_rows = length(find(ha_tmp_ind>0))-cond_cells_total;
    migrating_rows = find(ha_tmp_val > threshold);
    migration_ratio = length(migrating_rows)/initial_condition_rows*100;
    ha_speed(i,4) = migration_ratio;
    
end

%% Plot 
close all
clf

% Set the "0" concentration to be 1 for logarithmic plot
concentrations(1) = 1;

figure(1)

% Uncomment for plot without error bars
% semilogx(concentrations,ha_speed(:,1),'r-o','Linewidth',3)
% hold on
% semilogx(concentrations,ctr_speed(:,1),'b-o','Linewidth',3)

% Uncomment to plot all concentrations with error bars
% errorbar(concentrations,ha_speed(:,1),ha_speed(:,2),'r-o','Linewidth',3)
% hold on
% errorbar(concentrations,ctr_speed(:,1),ctr_speed(:,2),'b-o','Linewidth',3)

% Uncomment to plot concentrations up to 50 ug/ml with error bars
errorbar(concentrations(1:6),ha_speed(1:6,1),ha_speed(1:6,2),'r-o','Linewidth',3)
text(concentrations,ha_speed(:,1)-0.02,num2str(ha_speed(:,3)),'FontWeight','bold')
hold on
errorbar(concentrations(1:6),ctr_speed(1:6,1),ctr_speed(1:6,2),'b-o','Linewidth',3)
text(concentrations,ctr_speed(:,1)+0.02,num2str(ctr_speed(:,3)),'FontWeight','bold')
xlim([0 100])

set(gca,'XScale','log');
xlabel('Fibronectin Concentration [ug/mL]','FontSize',20);
ylabel('Cell Speed um/min','FontSize',20);
title(['MEF Cells - Combined Experiments - ',num2str(cells_total),' cells'])
grid on; 

% legend('Hase','Control');
legend(['Hase',' (',num2str(sum(ha_speed(:,3))),' cells)'],['Control ','(',num2str(sum(ctr_speed(:,3))),' cells)'],'Location','Southwest');
%%  Plot percentage of migrating cells per condition

figure(2)

semilogx(concentrations,ha_speed(:,4),'r-o','Linewidth',3)
hold on
semilogx(concentrations,ctr_speed(:,4),'b-o','Linewidth',3)

xlabel('Fibronectin Concentration [ug/mL]','FontSize',20);
ylabel('Cell Migration Ratio (%)','FontSize',20);
% title(['PC3 Cells; Threshold = ',num2str(threshold),'; Migration ratio: ',num2str(100*cells_mig/cells_total),'%' ])
title(['MEF Cells - Combined Experiments - ; Threshold = ',num2str(threshold),' um/min'])

grid on; legend(['Hase',' (',num2str(sum(ha_speed(:,3))),' cells)'],['Control ','(',num2str(sum(ctr_speed(:,3))),' cells)'],'Location','Southeast');