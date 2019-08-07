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
function [] = visualizeFeatures_beta(output, analysisplan)

% Import labels (cell and substrate type) from analysis plan
[~,experiment_labels,~] = xlsread(analysisplan,'labels');
label_cell_type = experiment_labels{1,2};
label_substrate_type = experiment_labels{2,2};

% cell_type = input('Input cell type: \n 1) PC3 \n 2) MEF \n 3) MCF-7 \n');

cond_vs_exp = input(['Plot \n 0) speed vs. concentration \n '...
    '1) speed vs. time by condition \n '...
    '2) speed vs. time by experiment \n '...
    '3) dist vs. time by condition \n '...
    '4) dist vs. time by experiment \n '...
    '5) delta theta distribution \n '...
    '6) speed vs. delta theta (all data)\n '...
    '7) speed vs. delta theta by condition \n '...
    '8) MSD \n '...
    '9) speed vs. delta theta by treatment\n '...
    '10) percent migration\n']);

% Set data analysis threshold (speed threshold and upper time limit)
speed_threshold = input('Enter desired speed threshold in um/min or d for default value: ','s');
if strcmp(speed_threshold, 'd')
    speed_threshold = 0.0;
else
    speed_threshold = str2double(speed_threshold);
end

time_threshold = input('Enter upper time limit for experiment in min or d for default value: ','s');
if strcmp(time_threshold, 'd')
    time_threshold = 180;
else
    time_threshold = str2double(time_threshold);
end

cd(output)

[num,~,raw] = xlsread(analysisplan,'experiments');
filenums = find(num(:,7)==0)+1;
treatment = num(filenums,2);

%initialize features variable
all_features = [];
control_features = [];
hadase_features = [];

% create map for decoding treatment
treatmentMap = containers.Map({1, 2, 0, 0.1, 0.5},{'w/ HA', 'w/o HA', '0 ug/mL Doxy', '0.1 ug/mL Doxy', '0.5 ug/mL Doxy'});

for curfilenum = 1:length(filenums)
    curfile = [num2str(raw{filenums(curfilenum),1}) '_features'];
    load(curfile);
    all_features = [all_features;features];
end

features = all_features;

% Calculate total number of cells

unique_cell_ind = find(isnan(features(:,4)));
cells_total = length(unique_cell_ind);

% Extract concentrations from data
concentrations = unique(features(:,7));
n_conc = length(concentrations);

% Extract treatments from data
treatments = unique(features(:,8));
n_treat = numel(treatments);

speed_matrix = zeros(n_conc*n_treat,6); % col1 = treatment, col2 = av, col3 = sem, col4 = total cells, col5 = subst. concentration, col6 = migration ratio
msd_matrix = zeros(n_conc*n_treat,5); % col1 = treatment, col2 = alpha, col3 = alpha_sem, col4 = rsq, col5 = subst. concentration
msd_ensemble_matrix = zeros(n_conc*n_treat,5);

% these variables will be deprecated
ctr_speed = zeros(n_conc,4); % col1 = av, col2 = sem, col3 = total cells, col4 = migration ratio
ha_speed = zeros(n_conc,4); % col1 = av, col2 = sem, col3 = total cells, col4 = migration ratio

slope = 0;
persis_coeff = [];

% initialize msd variables
msd_control = [];
msd_ha = [];
msd_ensemble_control = [];
msd_ensemble_ha = [];

row_ind=1;

for i = 1:n_conc
    
    if cond_vs_exp == 3
        figure
        hold on
    end
    
    if cond_vs_exp == 7
        persistence_fig_scatter = figure('units','normalized','outerposition',[0 0 1 1]);
        hold on
        persistence_fig_bar = figure('units','normalized','outerposition',[0 0 1 1]);
        hold on
        persistence_fig_hist = figure('units','normalized','outerposition',[0 0 1 1]);
        hold on        
    end
    
    conc_cell_total = 0; % total cell number per concentration (many treatments)
    
    % loop over treatments i.e. for j = 1:n_treat; use treatment_map
    % for labels!
    for j = 1:n_treat
        
        % Extract condition speeds (specific concentration and treatment with applied time and speed thresholds)
        tmp_treat = treatments(j);
        tmp_conc = concentrations(i);
        tmp_ind = features(:,7)== tmp_conc & features(:,8) == tmp_treat;
        tmp_val = features(tmp_ind,4);
        
        % proceed if the condition has data
%         if not(isempty(tmp_val))
        
            % Compute num of cells for condition
            unique_cell_ind = find(isnan(tmp_val));
            cond_cells_total = numel(unique_cell_ind);
            conc_cell_total = conc_cell_total + cond_cells_total;
            speed_matrix(row_ind,4) = cond_cells_total;

            % Implement thresholds
            tmp_ind = features(:,7)== tmp_conc & features(:,8) == tmp_treat & features(:,3) <= time_threshold & features(:,4) >= speed_threshold;
            tmp_val = features(tmp_ind,4);

            % Populate speed matrix
            speed_matrix(row_ind,1) = tmp_treat;
            speed_matrix(row_ind,5) = tmp_conc;
            speed_matrix(row_ind,2) = nanmean(tmp_val); % Compute average speed
            sample_size = length(tmp_val);
            speed_matrix(row_ind,3)= nanstd(tmp_val)/sqrt(sample_size); % Compute SEM

            % Migration ratio calculations
            if cond_vs_exp == 10
                tmp_ind = features(:,7)== tmp_conc & features(:,8) == tmp_treat;
                tmp_val = features(tmp_ind,4);
                initial_condition_rows = length(find(tmp_ind > 0))-cond_cells_total;
                migrating_rows = find(tmp_val >= speed_threshold);
                migration_ratio = numel(migrating_rows)/initial_condition_rows*100;
                speed_matrix(row_ind,6) = migration_ratio;
            end

            % plot speed vs. time for each condition (Treatment/Concentration)
            if cond_vs_exp == 1
                cond_features = features(tmp_ind,:);
                cond_slope = speed_vs_time_cond_beta(label_cell_type, treatmentMap(tmp_treat),cond_cells_total, concentrations(i),label_substrate_type,cond_features);
                slope = slope + cond_slope*cond_cells_total/cells_total;
            end

            % distance vs. time by condition
            if cond_vs_exp == 3
                cond_features = features(tmp_ind,:);
                [time, dist, equationString, rsqString, slope, SSresid, rsq] = dist_vs_time_cond_beta(label_cell_type, treatmentMap(tmp_treat), concentrations(i),cond_features);
                persis_coeff = [persis_coeff;tmp_treat, slope, SSresid];
            end

            % speed vs. delta theta by condition
            if cond_vs_exp == 7
                cond_features = features(tmp_ind,:);
                persistence_cond_beta(treatmentMap(tmp_treat), cond_features, persistence_fig_scatter, persistence_fig_bar, persistence_fig_hist)
            end

            % MSD analysis
            if cond_vs_exp == 8
                % we have to redo this to recover the NaN values (cell markers)
                tmp_ind = features(:,7)== tmp_conc & features(:,8) == tmp_treat;
                cond_features = features(tmp_ind,:);

                % process per cell
                [alpha_cond, alpha_sem_cond, rsq_cond] = msd_beta(cond_features, time_threshold);
                msd_matrix(row_ind,:) = [tmp_treat, alpha_cond, alpha_sem_cond,rsq_cond, tmp_conc];
    %             msd = [msd; tmp_treat, alpha_cond, alpha_sem_cond,rsq_cond];

                % process with ensemble statistics
                [alpha_ens_cond, alpha_ens_ssresid_cond, rsq_ens_cond] =  msd_ensemble_beta(label_cell_type, label_substrate_type, treatmentMap(tmp_treat), tmp_conc, cond_features, time_threshold);
                msd_ensemble_matrix(row_ind,:) = [tmp_treat, alpha_ens_cond, alpha_ens_ssresid_cond,rsq_ens_cond, tmp_conc];

            end

            row_ind = row_ind + 1;
%         end
        
    end
    
    if cond_vs_exp == 3
        hold off
        legend show
        legend('Location','southeast')
        grid on
        box on
        
        xlabel('Time (min)','FontSize',20)
        ylabel('Cell Distance (um)','FontSize',20);
        title([label_cell_type, ' cells; ', num2str(tmp_conc), ' ug/mL ', label_substrate_type,'; ' num2str(conc_cell_total),' cells'], 'FontSize',20)
    end
    
    if cond_vs_exp == 7
        
        figure(persistence_fig_scatter)
        
        hold off
        legend show
        box on
        
        ylim([0 5])
        xlabel(['\Delta \theta (', char(176), ')'],'FontSize',20)
        ylabel('Av. Speed (um/min)','FontSize',20)
        title([label_cell_type, ' cells; ', num2str(tmp_conc), ' ug/mL ', label_substrate_type,'; ' num2str(conc_cell_total),' cells'], 'FontSize',20)

        figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\speed_vs_delta_theta\by_condition';
        figure_file_name = [figure_directory, '\',label_cell_type, '_','scatter_',num2str(tmp_conc),'_',label_substrate_type,'.png'];
        print(persistence_fig_scatter,'-dpng',figure_file_name)
        
        figure(persistence_fig_bar)
        
        hold off
        grid on
        box on
        legend show

        xlabel(['\Delta \theta (', char(176), ')'],'FontSize',20)
        ylabel('Av. Speed (um/min)','FontSize',20)
        title([label_cell_type, ' cells; ', num2str(tmp_conc), ' ug/mL ', label_substrate_type,'; ' num2str(conc_cell_total),' cells'], 'FontSize',20)

        figure_file_name = [figure_directory, '\',label_cell_type, '_','bar_',num2str(tmp_conc),'_',label_substrate_type,'.png'];
        print(persistence_fig_bar,'-dpng',figure_file_name)

        figure(persistence_fig_hist)
        
        hold off
        box on
        legend show
        
        xlabel(['\Delta \theta (', char(176), ')'],'FontSize',20)
        title([label_cell_type, ' cells; ', num2str(tmp_conc), ' ug/mL ', label_substrate_type,'; ' num2str(conc_cell_total),' cells'], 'FontSize',20)
        
        figure_file_name = [figure_directory, '\',label_cell_type, '_','hist_',num2str(tmp_conc),'_',label_substrate_type,'.png'];
        print(persistence_fig_hist,'-dpng',figure_file_name)

    end
end

%% Plot speed vs. concentration
if cond_vs_exp == 0
    
    % offset the 0 ug/ml concentration for logarithmic plot (also because there
    % is some substrate in media)
    concentrations(concentrations == 0) = 1; 
    
    [xmin, xmax, ymin, ymax, ytext] = plotParam_beta(speed_matrix, n_treat);
    
%     figure
    figure('units','normalized','outerposition',[0 0 1 1])
    
    hold on;
    for i = 1:n_treat
        
        tmp_treat = treatments(i);
        treat_features = speed_matrix(speed_matrix(:,1)==tmp_treat,:);
        legend_label = [treatmentMap(tmp_treat),' (', num2str(sum(treat_features(:,4))), ' cells)'];
        
        erb = errorbar(concentrations, treat_features(:,2),treat_features(:,3),'Marker','.','MarkerSize', 18,'Linewidth',1.5, 'DisplayName',legend_label);
        text(concentrations,ytext(i)*ones(1,n_conc),num2str(treat_features(:,4)),'Color', erb.Color,'FontWeight','bold')
        
    end
    hold off
    grid on
    legend show
    legend('Location','southeast')
    box on
    
    xlim([xmin xmax])
    ylim([ymin ymax])

    set(gca,'XScale','log','FontSize',20);
    xlabel([label_substrate_type, ' Concentration (ug/ml)'], 'FontSize',20)
    ylabel('Cell Speed (um/min)','FontSize',20);
    title([label_cell_type, ' cells; ',num2str(cells_total),' cells'], 'FontSize',20)
    
end

%% Output weighted slope for 

if cond_vs_exp == 1
   display(['Weighted slope (speed vs. time) is: ', num2str(slope)]) 
end

%% Plot persistence coefficient vs. concentration

% persis_coeff(:,2) = (1-persis_coeff(:,2)).*persis_coeff(:,1);
if cond_vs_exp == 3
    
    % offset the 0 ug/ml concentration for logarithmic plot (also because there
    % is some substrate in media)
    concentrations(concentrations == 0) = 1; 
    
%     [xmin, xmax, ymin, ymax, ytext] = plotParam_beta(persis_coeff, n_treat);
    
%     figure(figure('units','normalized','outerposition',[0 0 1 1]))
    figure
    
    hold on;
    for i = 1:n_treat
        
        tmp_treat = treatments(i);
        treat_features = speed_matrix(speed_matrix(:,1)==tmp_treat,:);
        treat_pers_coeff = persis_coeff(persis_coeff(:,1) == tmp_treat,:);
        legend_label = [treatmentMap(tmp_treat),' (', num2str(sum(treat_features(:,4))), ' cells)'];
        
        erb = errorbar(concentrations, treat_pers_coeff(:,2),treat_pers_coeff(:,3),'Marker','.','MarkerSize', 18,'Linewidth',1.5, 'DisplayName',legend_label);
        
    end
    hold off
    grid on
    legend show
    legend('Location','southeast')
    box on
    
%     xlim([xmin xmax])
%     ylim([ymin ymax])

    set(gca,'XScale','log','FontSize',20);
    xlabel([label_substrate_type, ' Concentration (ug/ml)'], 'FontSize',20)
    ylabel('Persist. coeff. (um/min^{1/2})','FontSize',20);
    title([label_cell_type, ' cells; ',num2str(cells_total),' cells'], 'FontSize',20)
    
end


%%  Plot percentage of migrating cells per condition
if cond_vs_exp == 10
      
    % offset the 0 ug/ml concentration for logarithmic plot (also because there
    % is some substrate in media)
    concentrations(concentrations == 0) = 1; 
    
    [xmin, xmax, ymin, ymax, ytext] = plotParam_beta(speed_matrix, n_treat);
    
%     figure(figure('units','normalized','outerposition',[0 0 1 1]))
    figure
    
    hold on;
    for i = 1:n_treat
        
        tmp_treat = treatments(i);
        treat_features = speed_matrix(speed_matrix(:,1)==tmp_treat,:);
        legend_label = [treatmentMap(tmp_treat),' (', num2str(sum(treat_features(:,4))), ' cells)'];
        
        plot(concentrations,treat_features(:,6),'-o','Linewidth',1.5, 'DisplayName',legend_label)
        
    end
    hold off
    grid on
    box on
    legend show
    legend('Location','southeast')
    
    xlim([xmin xmax])
    ylim([0 100])
    set(gca,'XScale','log');

    xlabel([label_substrate_type, ' Concentration (ug/ml)'], 'FontSize',20)
    ylabel('Cell Migration Ratio (%)','FontSize',20);
    title([label_cell_type, ' cells; ',num2str(cells_total),' cells; Threshold = ',num2str(speed_threshold),' um/min'], 'FontSize',20)
    
end

%% Plot speed vs time per experiment

if cond_vs_exp == 2
    
    slope = 0;
    
    for curfilenum = 1:length(filenums)
        experiment = num2str(raw{filenums(curfilenum),1});
        curfile = [experiment '_features'];
        load(curfile);

        % Calculate total number of cells per experiment

        unique_cell_ind = find(isnan(features(:,4)));
        cells_experiment = numel(unique_cell_ind);
        
        % update slope

        experiment_slope = speed_vs_time_exp_beta(experiment, features, label_cell_type, label_substrate_type, cells_experiment);
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

        experiment_slope = dist_vs_time_exp(label_cell_type, label_substrate_type, experiment, features, cells_experiment);
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
    title([label_cell_type, ' \Delta \theta distribution'])
    
end

%% Plot speed vs. delta theta

if cond_vs_exp == 6

    figure 
    
    scatter(features(:,5),features(:,6),'o')
    title(label_cell_type,'FontSize',20)
    xlabel(['\Delta \theta (', char(176), ')'],'FontSize',20)
    ylabel('Cell Speed (um/min)','FontSize',20)
    
    box on

end

%% Display MSD analysis

if cond_vs_exp == 8
    
    % offset the 0 ug/ml concentration for logarithmic plot (also because there
    % is some substrate in media)
    concentrations(concentrations == 0) = 1; 
    
    [xmin, xmax, ymin, ymax, ytext] = plotParam_beta(msd_matrix, n_treat);
    
    figure
    
    hold on;
    for i = 1:n_treat
        
        tmp_treat = treatments(i);
        treat_msd = msd_matrix(msd_matrix(:,1)==tmp_treat,:);
        
        legend_label = treatmentMap(tmp_treat);
        errorbar(concentrations, treat_msd(:,2),treat_msd(:,3),'Linewidth',1.5, 'DisplayName',legend_label);
      
        
    end
    hold off
    grid on
    legend show
    box on
    
    xlim([xmin xmax])
    ylim([ymin ymax])
    
    title([label_cell_type, ' Cells - MSD analysis per cell'])
    xlabel([label_substrate_type, ' Concentration (ug/ml)'])
    ylabel('\alpha value','FontSize',20);
    
    set(gca,'XScale','log','FontSize',20);
    
end

%% Plot speed vs. delta theta for all treatments

if cond_vs_exp == 9
    
    persistence_fig_scatter = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    persistence_fig_bar = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    persistence_fig_hist = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    
    for i = 1:n_treat
        
        tmp_treat = treatments(i);
        tmp_ind = features(:,8) == tmp_treat;
        treat_pers_features = features(tmp_ind,:);
        % Compute num of cells for treatment
        tmp_val = features(tmp_ind,4); % need this for the NaN cell indicators
        unique_cell_ind = find(isnan(tmp_val));
        treat_cells_total = numel(unique_cell_ind);
        
        persistence_cond_beta(treatmentMap(tmp_treat), treat_pers_features, persistence_fig_scatter, persistence_fig_bar, persistence_fig_hist)
        
    end
    
    figure(persistence_fig_scatter)

    hold off
    legend show
    box on

    ylim([0 5])
    xlabel(['\Delta \theta (', char(176), ')'],'FontSize',20)
    ylabel('Av. Speed (um/min)','FontSize',20)
    title([label_cell_type, ' Cells on ', label_substrate_type], 'FontSize',20)

    figure_directory = '\\PHYS34212\MigrationData\MigrationData\Migration1\figures\speed_vs_delta_theta\by_condition';
    figure_file_name = [figure_directory, '\',label_cell_type, '_','scatter_',label_substrate_type,'.png'];
    print(persistence_fig_scatter,'-dpng',figure_file_name)

    figure(persistence_fig_bar)

    hold off
    grid on
    box on
    legend show

    xlabel(['\Delta \theta (', char(176), ')'],'FontSize',20)
    ylabel('Av. Speed (um/min)','FontSize',20)
    title([label_cell_type, ' Cells on ', label_substrate_type], 'FontSize',20)

    figure_file_name = [figure_directory, '\',label_cell_type, '_','bar_',label_substrate_type,'.png'];
    print(persistence_fig_bar,'-dpng',figure_file_name)

    figure(persistence_fig_hist)

    hold off
    box on
    legend show

    xlabel(['\Delta \theta (', char(176), ')'],'FontSize',20)
    title([label_cell_type, ' Cells on ', label_substrate_type], 'FontSize',20)

    figure_file_name = [figure_directory, '\',label_cell_type, '_','hist_',label_substrate_type,'.png'];
    print(persistence_fig_hist,'-dpng',figure_file_name)
end