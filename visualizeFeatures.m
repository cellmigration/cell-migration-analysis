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

cond_vs_exp = input('Plot 1)speed vs time by condition 2) speed vs time by experiment? 3) dist vs. time by condition 4) dist vs. time by experiment 5) angle vs. time by condition 6) delta theta vs. distance');
conc_range = input('Plot 1) all concentrations or 2) up to 75 ug/mL? ');

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

slope = 0;
persis_coeff = [];

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
    
    % plot speed vs. time for each condition
    if cond_vs_exp == 1
        ctrl_cond_features = features(ctrl_tmp_ind,:);
        slope_ctrl = speed_vs_time_cond('Control',concentrations(i),ctrl_cond_features);
        slope = slope + slope_ctrl*cond_cells_total/cells_total;
    end
    
    % plot dist vs. time for each condition
    
    if cond_vs_exp == 3
        ctrl_cond_features = features(ctrl_tmp_ind,:);
        coeff_rsq = dist_vs_time_cond('Control',concentrations(i),ctrl_cond_features);
        persis_coeff = [persis_coeff; coeff_rsq];
    end
    
    % plot msd persistence for each condition
    % note; only useful for global behavior
    
%     if cond_vs_exp == 5
%         ctrl_cond_features = features(ctrl_tmp_ind,:);
%         msd_cond = msd('Control',concentrations(i),ctrl_cond_features);
%     end

    if cond_vs_exp == 5
        ctrl_cond_features = features(ctrl_tmp_ind,:);
        angle_vs_time_cond('Control',concentrations(i),ctrl_cond_features);
    end
    
    
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
    
    % plot speed vs. time for each condition
    if cond_vs_exp ==1
        ha_cond_features = features(ha_tmp_ind,:);
        slope_has = speed_vs_time_cond('Hase',concentrations(i),ha_cond_features);
        slope = slope + slope_has*cond_cells_total/cells_total;
    end
    
%     if cond_vs_exp == 3
%         ha_cond_features = features(ha_tmp_ind,:);
%         slope_has = dist_vs_time_cond('Hase',concentrations(i),ha_cond_features);
%         slope = slope + slope_has*cond_cells_total/cells_total;
%     end
    
    % migration calculations
    initial_condition_rows = length(find(ha_tmp_ind>0))-cond_cells_total;
    migrating_rows = find(ha_tmp_val > threshold);
    migration_ratio = length(migrating_rows)/initial_condition_rows*100;
    ha_speed(i,4) = migration_ratio;
    
end

%% Plot speed vs concentration for each condition
% close all
% clf

% Set the "0" concentration to be 1 for logarithmic plot
concentrations(1) = 1;

figure


% Uncomment for plot without error bars
% semilogx(concentrations,ha_speed(:,1),'r-o','Linewidth',3)
% hold on
% semilogx(concentrations,ctr_speed(:,1),'b-o','Linewidth',3)

% Uncomment to plot all concentrations with error bars ~ text near points
% errorbar(concentrations,ha_speed(:,1),ha_speed(:,2),'r-o','Linewidth',3)
% hold on
% errorbar(concentrations,ctr_speed(:,1),ctr_speed(:,2),'b-o','Linewidth',3)
% text(concentrations,ctr_speed(:,1)+0.02,num2str(ctr_speed(:,3)),'FontWeight','bold')
% ylim([0.65 1.07])

% Uncomment to plot all concentrations with error bars ~ text aligned, automatic plot limits
if conc_range == 1
    
    [ymin, ymax, ytext_ctr, ytext_ha] = plotParam(ctr_speed, ha_speed);

    errorbar(concentrations,ha_speed(:,1),ha_speed(:,2),'r','Linewidth',1.5)
    text(concentrations,ytext_ha*ones(1,n_conc),num2str(ha_speed(:,3)),'Color','red','FontWeight','bold')

    hold on
    errorbar(concentrations,ctr_speed(:,1),ctr_speed(:,2),'b','Linewidth',1.5)
    text(concentrations,ytext_ctr*ones(1,n_conc),num2str(ctr_speed(:,3)),'Color','blue','FontWeight','bold')
    ylim([ymin ymax])
    xlim([0.65 650])
end

% Uncomment to plot concentrations up to 50 ug/ml with error bars ~ text near points
% errorbar(concentrations(1:6),ha_speed(1:6,1),ha_speed(1:6,2),'r-o','Linewidth',3)
% text(concentrations(1:6),ha_speed(1:6,1)-0.02,num2str(ha_speed(1:6,3)),'FontWeight','bold')
% hold on
% errorbar(concentrations(1:6),ctr_speed(1:6,1),ctr_speed(1:6,2),'b-o','Linewidth',3)
% text(concentrations(1:6),ctr_speed(1:6,1)+0.02,num2str(ctr_speed(1:6,3)),'FontWeight','bold')
% text(concentrations(1:6),1.05*ones(1,6),num2str(ctr_speed(1:6,3)),'FontWeight','bold')
% xlim([0 100])
% ylim([0.7 1.07])

% Uncomment to plot concentrations up to 50 ug/ml with error bars ~ text aligned, automatic plot limits
if conc_range == 2
    
    [ymin, ymax, ytext_ctr, ytext_ha] = plotParam(ctr_speed(1:6,:), ha_speed(1:6,:));

    errorbar(concentrations(1:6),ha_speed(1:6,1),ha_speed(1:6,2),'r','Linewidth',1.5)
    text(concentrations(1:6),ytext_ha*ones(1,6),num2str(ha_speed(1:6,3)),'Color','red','FontWeight','bold')

    hold on

    errorbar(concentrations(1:6),ctr_speed(1:6,1),ctr_speed(1:6,2),'b','Linewidth',1.5)
    text(concentrations(1:6),ytext_ctr*ones(1,6),num2str(ctr_speed(1:6,3)),'Color','blue','FontWeight','bold')

    xlim([0.65 100])
    ylim([ymin ymax])
end

set(gca,'XScale','log');
xlabel('Fibronectin Concentration [ug/mL]','FontSize',20);
ylabel('Cell Speed [um/min]','FontSize',20);
title(['MEF Cells - Combined Experiments - ',num2str(cells_total),' cells'])
grid on; 

% legend('Hase','Control');
legend(['Hase',' (',num2str(sum(ha_speed(:,3))),' cells)'],['Control ','(',num2str(sum(ctr_speed(:,3))),' cells)'],'Location','Southwest');

%% Plot persistence coefficient vs. concentration

% persis_coeff(:,2) = (1-persis_coeff(:,2)).*persis_coeff(:,1);
if cond_vs_exp == 3
    figure
    % plot(concentrations, persis_coeff(:,1),'Linewidth',1.5)
    errorbar(concentrations(1:6,:), persis_coeff(1:6,1),persis_coeff(1:6,2),'Linewidth',1.5)

    set(gca,'XScale','log');
    xlabel('Fibronectin Concentration [ug/mL]','FontSize',20);
    ylabel('Persist. coeff. [um/min^{1/2}]','FontSize',20);
    title(['MEF Cells - Combined Experiments - ',num2str(cells_total),' cells'])
    xlim([0.65 100])
    grid on; 
end


%%  Plot percentage of migrating cells per condition

figure

semilogx(concentrations,ha_speed(:,4),'r-o','Linewidth',3)
hold on
semilogx(concentrations,ctr_speed(:,4),'b-o','Linewidth',3)

xlabel('Fibronectin Concentration [ug/mL]','FontSize',20);
ylabel('Cell Migration Ratio (%)','FontSize',20);
% title(['PC3 Cells; Threshold = ',num2str(threshold),'; Migration ratio: ',num2str(100*cells_mig/cells_total),'%' ])
title(['MEF Cells - Combined Experiments - ; Threshold = ',num2str(threshold),' um/min'])

grid on; legend(['Hase',' (',num2str(sum(ha_speed(:,3))),' cells)'],['Control ','(',num2str(sum(ctr_speed(:,3))),' cells)'],'Location','Southeast');


%% Plot speed vs time per experiment

if cond_vs_exp == 2
    
    slope = 0;
    
    for curfilenum = 1:length(filenums)
        experiment = num2str(raw{filenums(curfilenum),1});
        curfile = [experiment '_features'];
        load(curfile);

        % Calculate total number of cells per experiment

        unique_cell_ind = find(isnan(features(:,4)));
        cells_experiment = length(unique_cell_ind);
        
        % update slope

        experiment_slope = speed_vs_time_exp(experiment, features);
        slope = slope + experiment_slope*cells_experiment;
    end

    slope = slope/cells_total;
    
    display(['Weighted slope (speed vs. time) is: ', num2str(slope)])
    
end



%% Plot distance vs. time per experiment

if cond_vs_exp == 4
    
    slope = 0;
    
    for curfilenum = 1:length(filenums)
        experiment = num2str(raw{filenums(curfilenum),1});
        curfile = [experiment '_features'];
        load(curfile);

        % Calculate total number of cells per experiment

        unique_cell_ind = find(isnan(features(:,4)));
        cells_experiment = length(unique_cell_ind);
        
        % update slope

        experiment_slope = dist_vs_time_exp(experiment, features);
        slope = slope + experiment_slope*cells_experiment;
    end

    slope = slope/cells_total;
    
    display(['Weighted slope (dist vs. time) is: ', num2str(slope)])
    
end

%% Plot distribution of delta thetas for all data

if cond_vs_exp == 5
    non_nan_rows = not(isnan(features(:,5)));
%     migrating_rows = non_nan_rows(find(features(:,4) > 0.5));
    all_angles = features(non_nan_rows,5); %angles excluding spurious 0's corresponding to the first two data points for each cell
    figure
    polarhistogram(all_angles,12)
    title('MEF cells \Delta \theta distribution')
    
end

%% Plot delta theta vs. distance 

if cond_vs_exp == 6
    figure
    scatter(features(:,6),features(:,5))
    xlabel('Distance [um]','FontSize',20)
    ylabel(['Delta theta [', char(176), ']'],'FontSize',20)
%     xlim([0 4])
    figure
    scatter(features(:,6)/10,features(:,5)) % the time step between points 1 and 3, 2 and 4 etc.. is 10  minutes
    title('MEF Cells','FontSize',20)
    xlabel('Speed [um/min]','FontSize',20)
    ylabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
    xlim([0 4])
end


    