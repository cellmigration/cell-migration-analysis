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

cell_type = input('Input cell type: \n 1) PC3 \n 2) MEF \n 3) MCF-7 \n');
substrate_type = input('Input substrate: \n 1) Collagen \n 2) Fibronectin \n');

switch cell_type
    case 1
        label_cell_type = 'PC3';
    case 2
        label_cell_type = 'MEF';
    case 3
        label_cell_type = 'MCF-7';
    otherwise
        label_cell_type = 'Unknown';
end

switch substrate_type
    case 1
        label_substrate_type = 'Collagen';
    case 2
        label_substrate_type = 'Fibronectin';
    otherwise
        label_substrate_type = 'Unknown';
end

cond_vs_exp = input('Plot \n 0) speed vs. concentration \n 1) speed vs time by condition \n 2) speed vs time by experiment? \n 3) dist vs. time by condition \n 4) dist vs. time by experiment \n 5) angle vs. time by condition \n 6) delta theta vs. distance \n 7) speed vs. delta theta by condition \n 8) MSD \n 9) speed vs. delta theta all conditions\n 10) MCF7 cells\n');

cd(output)

[num,~,raw] = xlsread(analysisplan,'experiments');
filenums = find(num(:,7)==0)+1;
treatment = num(filenums,2);

%initialize features variable
all_features = [];
control_features = [];
hadase_features = [];


for curfilenum = 1:length(filenums)
    curfile = [num2str(raw{filenums(curfilenum),1}) '_features'];
    load(curfile);
    all_features = [all_features;features];
    if treatment(curfilenum) == 1
        control_features = [control_features; features];
    else % treatment(curfilenum) ==2
        hadase_features = [hadase_features; features];
    end
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
doxy_0p5_speed = zeros(n_conc,4); % col1 = av, col2 = sem, col3 = total cells, col4 = migration ratio
doxy_0p1_speed = zeros(n_conc,4);
doxy_0_speed = zeros(n_conc,4);

slope = 0;
persis_coeff = [];

% initialize msd variables
msd_control = [];
msd_ha = [];
msd_ensemble_control = [];
msd_ensemble_ha = [];

for i = 1:n_conc
    
    % Calculations for PC3 or MEF cells (eperiments using HAdase)
    if cell_type == 1 || cell_type ==2
    
        % Calculations for control cells
        ctrl_tmp_ind = features(:,7)==concentrations(i) & features(:,8) == 1;
        ctrl_tmp_val = features(ctrl_tmp_ind,4);
%         ctrl_tmp_val = ctrl_tmp_val(ctrl_tmp_val > threshold);
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



        % plot msd persistence for each condition
        % note; only useful for global behavior

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
%         ha_tmp_val = ha_tmp_val(ha_tmp_val > threshold);
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

        if cond_vs_exp == 3
            ctrl_cond_features = features(ctrl_tmp_ind,:);
            ha_cond_features = features(ha_tmp_ind,:);
            [coeff_rsq_ha, coeff_rsq_control] = dist_vs_time_cond2(concentrations(i),ha_cond_features, ctrl_cond_features);
            persis_coeff = [persis_coeff; coeff_rsq_ha, coeff_rsq_control];
        end

        if cond_vs_exp == 7
            ha_cond_features = features(ha_tmp_ind,:);
            ctrl_cond_features = features(ctrl_tmp_ind,:);
            persistence_cond2(concentrations(i), ha_cond_features, ctrl_cond_features)
        end

        if cond_vs_exp == 8
            ctrl_cond_features = features(ctrl_tmp_ind,:);
            ha_cond_features = features(ha_tmp_ind,:);
            % process per cell
            [alpha_cond_control, alpha_sem_cond_control, rsq_cond_control] = msd(ctrl_cond_features);
            [alpha_cond_ha, alpha_sem_cond_ha, rsq_cond_ha] = msd(ha_cond_features);
            msd_control = [msd_control;concentrations(i),alpha_cond_control,alpha_sem_cond_control,rsq_cond_control];
            msd_ha = [msd_ha;concentrations(i),alpha_cond_ha,alpha_sem_cond_ha, rsq_cond_ha];
            % process with ensemble statistics
            [alpha_ens_cond_control, alpha_ens_ssresid_cond_control, rsq_ens_cond_control] =  msd_ensemble('Control', concentrations(i), ctrl_cond_features);
            [alpha_ens_cond_ha, alpha_ens_ssresid_cond_ha, rsq_ens_cond_ha] =  msd_ensemble('Hadase', concentrations(i), ha_cond_features);
            msd_ensemble_control = [msd_ensemble_control; concentrations(i),alpha_ens_cond_control, alpha_ens_ssresid_cond_control, rsq_ens_cond_control];
            msd_ensemble_ha = [msd_ensemble_ha; concentrations(i), alpha_ens_cond_ha, alpha_ens_ssresid_cond_ha, rsq_ens_cond_ha];
        end
    
    elseif cell_type == 3
        
        % calculations for doxy = 0.5
        doxy_0p5_tmp_ind = features(:,7)==concentrations(i) & features(:,8) == 0.5;
        doxy_0p5_tmp_val = features(doxy_0p5_tmp_ind,4);
        doxy_0p5_speed(i,1) = nanmean(doxy_0p5_tmp_val);
        sample_size = length(doxy_0p5_tmp_val);
        doxy_0p5_speed(i,2)= nanstd(doxy_0p5_tmp_val)/sqrt(sample_size);
        
        % compute total number of cells
        unique_cell_ind = find(isnan(doxy_0p5_tmp_val));
        cond_cells_total = length(unique_cell_ind);
        doxy_0p5_speed(i,3) = cond_cells_total;
        
        %calculations for doxy = 0.1
        doxy_0p1_tmp_ind = features(:,7)==concentrations(i) & features(:,8) == 0.1;
        doxy_0p1_tmp_val = features(doxy_0p1_tmp_ind,4);
        doxy_0p1_speed(i,1) = nanmean(doxy_0p1_tmp_val);
        sample_size = length(doxy_0p1_tmp_val);
        doxy_0p1_speed(i,2)= nanstd(doxy_0p1_tmp_val)/sqrt(sample_size);
        
        % compute total number of cells
        unique_cell_ind = find(isnan(doxy_0p1_tmp_val));
        cond_cells_total = length(unique_cell_ind);
        doxy_0p1_speed(i,3) = cond_cells_total;
        
        %calculations for doxy = 0
        doxy_0_tmp_ind = features(:,7)==concentrations(i) & features(:,8) == 0;
        doxy_0_tmp_val = features(doxy_0_tmp_ind,4);
        doxy_0_speed(i,1) = nanmean(doxy_0_tmp_val);
        sample_size = length(doxy_0_tmp_val);
        doxy_0_speed(i,2)= nanstd(doxy_0_tmp_val)/sqrt(sample_size);
        
        % compute total number of cells
        unique_cell_ind = find(isnan(doxy_0_tmp_val));
        cond_cells_total = length(unique_cell_ind);
        doxy_0_speed(i,3) = cond_cells_total;
        
    end    
end

% Plot results for MCF7 cells
if cond_vs_exp == 10
    
    concentrations(1) = 1;
    
    figure(107)
        
%     errorbar(concentrations,doxy_0p5_speed(:,1),doxy_0p5_speed(:,2),'-o','Linewidth',1.5)
    hold on;
%     errorbar(concentrations,doxy_0p1_speed(:,1),doxy_0p1_speed(:,2),'-o','Linewidth',1.5)
    errorbar(concentrations,doxy_0_speed(:,1),doxy_0_speed(:,2),'r','Linewidth',1.5)
    xlabel([label_substrate_type, ' Concentration (ug/ml)'], 'FontSize',20)
    ylabel('Cell Speed (um/min)','FontSize',20);
        
    xlim([0.95 1050])
    set(gca,'XScale','log');

%     title([label_cell_type, ' cells; ',num2str(cells_total),' cells'], 'FontSize',16)
    title([label_cell_type, ' cells; ',num2str(cells_total),' cells'])
    grid on;
%     legend(['0.5 Doxy',' (', num2str(sum(doxy_0p5_speed(:,3))), ' cells)'],['0.1 Doxy',' (', num2str(sum(doxy_0p1_speed(:,3))), ' cells)'],['0 Doxy',' (', num2str(sum(doxy_0_speed(:,3))), ' cells)']', 'FontSize', 20, 'Location','Northeast');
    legend(['0 ug/ml Doxy',' (', num2str(sum(doxy_0_speed(:,3))), ' cells)']','Location','Southwest');

    set(gca,'FontSize',20)
    
end

%% Plot speed vs concentration for each condition
% close all
% clf

% Set the "0" concentration to be 1 for logarithmic plot

if cond_vs_exp == 0
    conc_range = input('Plot 1) all concentrations or 2) up to 75 ug/mL? ');

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
%         ylim([0 1])
%         xlim([0.65 150])
        ylim([ymin ymax])
        xlim([0.65 1500])
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

        num_concentrations = 6;

        [ymin, ymax, ytext_ctr, ytext_ha] = plotParam(ctr_speed(1:num_concentrations,:), ha_speed(1:num_concentrations,:));

        errorbar(concentrations(1:num_concentrations),ha_speed(1:num_concentrations,1),ha_speed(1:num_concentrations,2),'r','Linewidth',1.5)
        text(concentrations(1:num_concentrations),ytext_ha*ones(1,num_concentrations),num2str(ha_speed(1:num_concentrations,3)),'Color','red','FontWeight','bold')

        hold on

        errorbar(concentrations(1:num_concentrations),ctr_speed(1:num_concentrations,1),ctr_speed(1:num_concentrations,2),'b','Linewidth',1.5)
        text(concentrations(1:num_concentrations),ytext_ctr*ones(1,num_concentrations),num2str(ctr_speed(1:num_concentrations,3)),'Color','blue','FontWeight','bold')
%         ylim([0 1])
%         xlim([0.65 150])
        xlim([0.65 100])
        ylim([ymin ymax])
    end

    set(gca,'XScale','log');
    xlabel([label_substrate_type,' Concentration [ug/mL]'],'FontSize',20);
    ylabel('Cell Speed [um/min]','FontSize',20);
    title([label_cell_type,' Cells - Combined Experiments - ',num2str(cells_total),' cells'])
    grid on; 

    % legend('Hase','Control');
    legend(['w/o HA',' (',num2str(sum(ha_speed(:,3))),' cells)'],['w/ HA',' (',num2str(sum(ctr_speed(:,3))),' cells)'],'Location','Southwest');
    
end

%% Plot persistence coefficient vs. concentration

% persis_coeff(:,2) = (1-persis_coeff(:,2)).*persis_coeff(:,1);
if cond_vs_exp == 3
    figure
    % plot(concentrations, persis_coeff(:,1),'Linewidth',1.5)
    concentrations(1) = 1;
    
    num_concentrations = 5;
    
    errorbar(concentrations(1:num_concentrations,:), persis_coeff(1:num_concentrations,1),persis_coeff(1:num_concentrations,2),'Linewidth',1.5, 'color',[0.8500, 0.3250, 0.0980])
    hold on
    errorbar(concentrations(1:num_concentrations,:), persis_coeff(1:num_concentrations,4),persis_coeff(1:num_concentrations,5),'Linewidth',1.5, 'color', [0, 0.4470, 0.7410])

    set(gca,'XScale','log');
    xlabel([label_substrate_type,' Concentration [ug/mL]'],'FontSize',20);
    ylabel('Persist. coeff. [um/min^{1/2}]','FontSize',20);
    title([label_cell_type,' Cells - Combined Experiments - ',num2str(cells_total),' cells'])
    xlim([0.65 1500])
    grid on; 
end


%%  Plot percentage of migrating cells per condition

figure

semilogx(concentrations,ha_speed(:,4),'r-o','Linewidth',1.5)
hold on
semilogx(concentrations,ctr_speed(:,4),'b-o','Linewidth',1.5)

xlim([0.65 100])
ylim([0 100])

xlabel([label_substrate_type,' Concentration [ug/mL]'],'FontSize',20);
ylabel('Cell Migration Ratio [%]','FontSize',20);
title([label_cell_type,' Cells; Threshold = ',num2str(threshold),' um/min']);%'; Migration ratio: ',num2str(100*cells_mig/cells_total),'%' 

grid on; legend(['w/o HA',' (',num2str(sum(ha_speed(:,3))),' cells)'],['w/ HA ','(',num2str(sum(ctr_speed(:,3))),' cells)'],'Location','Northeast');


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

    figure(233)
    scatter(features(:,5),features(:,6),'ro') % the time step between points 1 and 3, 2 and 4 etc.. is 10  minutes
    title('MEF Cells','FontSize',20)
    xlabel('Speed [um/min]','FontSize',20)
    ylabel(['\Delta \theta [', char(176), ']'],'FontSize',20)
    hold on
    hadase_pers = [features(:,5) features(:,6)];
end

%% Display MSD analysis

if cond_vs_exp == 8
    table(msd_control)
    table(msd_ha)
    
    figure
    errorbar(msd_control(:,1), msd_control(:,2), msd_control(:,3), 'Linewidth',1.5)
    hold on
    errorbar(msd_ha(:,1), msd_ha(:,2), msd_ha(:,3), 'Linewidth',1.5)
    hold off
    set(gca,'XScale','log');
    xlim([0.3 1500])
    
%     semilogx(msd_control(:,1),msd_control(:,2),msd_ha(:,1),msd_ha(:,2),'Linewidth',1.5)
    
    title('PC3 Cells - MSD analysis per cell')
    xlabel('Collagen Concentration [ug/mL]','FontSize',20);
    ylabel('\alpha value','FontSize',20);
    legend('Control', 'Hase')
    
    figure
    errorbar(msd_ensemble_control(:,1), msd_ensemble_control(:,2), msd_ensemble_control(:,3), 'Linewidth',1.5)
    hold on
    errorbar(msd_ensemble_ha(:,1), msd_ensemble_ha(:,2), msd_ensemble_ha(:,3), 'Linewidth',1.5)
    hold off
    set(gca,'XScale','log');
    xlim([0.3 1500])
    
%     semilogx(msd_control(:,1),msd_control(:,2),msd_ha(:,1),msd_ha(:,2),'Linewidth',1.5)
    
    title('PC3 Cells - MSD ensemble analysis')
    xlabel('Collagen Concentration [ug/mL]','FontSize',20);
    ylabel('\alpha value','FontSize',20);
    legend('Control', 'Hase')
    
end

%% Plot speed vs. delta theta for all conditions

if cond_vs_exp == 9
    control_pers = [control_features(:,5) control_features(:,6)];
    hadase_pers = [hadase_features(:,5) hadase_features(:,6)];
    persistence_analysis(control_pers, hadase_pers)
end