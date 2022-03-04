%% plot figure 4C of the manuscript

%% Important: run code while being in folder 'grasp_and_speech_decoding'

clc
clear all 
close all

%% changeable parameters

SavedData = [pwd '\Data\ClassificationComparisonMotorImagerySpeech\'];

number_sessions = 5;
number_phases = 4;

unit_regions = {'SMG', 'PMV', 'S1X'};
number_regions = length(unit_regions);
phase_names = {'ITI', 'Cue', 'Delay', 'Action'};
speech_cues = {'MotorImagery','Grasps', 'Colors'};
number_tasks = length(speech_cues);

keep_means = cell(number_regions,number_tasks);
keep_shuffled_means = cell(number_regions,number_tasks);
keep_sig_matrix = cell(number_regions,number_tasks);
prc_val_test = [95, 99, 99.9]; %percentile values to calculate

%loop through all brain regions
for n_reg = 1:number_regions
    %loop through three different tasks
    for n__task = 1:number_tasks
        %load data
        data = load([SavedData [unit_regions{n_reg} '_' speech_cues{n__task} '_Error']]);
        Shuffle = load([SavedData [unit_regions{n_reg} '_' speech_cues{n__task} '_ErrorShuffle']]);
        
        %preallocate variables
        keep_mean_acc = zeros(number_phases,number_sessions);
        keep_mean_acc_shuffle = zeros(number_phases,number_sessions);
        sig_matrix = zeros(number_phases, length(prc_val_test));
        %calculate percentile significance for each phase
        for n_phase = 1:number_phases

            err_test = squeeze(data.errTest(:,:,n_phase,:));
            acc = (1 - err_test)*100;
            testing_err_shuffle = squeeze(Shuffle.errTest(:,:,n_phase,:));
            mean_acc = mean(acc);

            acc_shuffle = (1 - squeeze(mean(mean(testing_err_shuffle)))')*100;
            keep_mean_acc(n_phase,:) = mean_acc;
            keep_mean_acc_shuffle(n_phase,:) = acc_shuffle;
            %calculate if actual data is significantly different from shuffled data
            prc_shuffle = (1 - reshape(mean(testing_err_shuffle,2), [], 1))*100;
            for prc = 1:length(prc_val_test)
                prc_test_val = prc_val_test(prc);
                prc_val = prctile(prc_shuffle,prc_test_val);
                
                if mean(mean_acc) < prc_val
                    sig_matrix(n_phase,prc) = 0;
                else
                    sig_matrix(n_phase,prc) = 1;
                    disp([unit_regions{n_reg} ' is significant at ' num2str(prc_test_val)  ' p during '  phase_names{n_phase} ' phase ' speech_cues{n__task}])
                end 
            end                    
        end
        keep_means{n_reg,n__task} = keep_mean_acc;
        keep_shuffled_means{n_reg,n__task} = keep_mean_acc_shuffle;
        keep_sig_matrix{n_reg,n__task} = sig_matrix;
    end
end 
 %% plot figure 4C
phase_to_compute = 4; 
%extract data for selected phase
data_mean_phase = cell2mat(cellfun(@(x) mean(x(phase_to_compute,:)), keep_means, 'UniformOutput', false));
data_phase = cellfun(@(x) x(phase_to_compute,:), keep_means, 'UniformOutput', false);
data_shuffle_phase = cellfun(@(x) x(phase_to_compute,:), keep_shuffled_means, 'UniformOutput', false);
err_ci_phase = cellfun(@(x) calculate_err_ci(x(phase_to_compute,:)), keep_means, 'UniformOutput', false);
% extract significant percentile values for phase
keep_sig_mat_phase = cellfun(@(x) find(x(phase_to_compute,:),1,'last'), keep_sig_matrix, 'UniformOutput', false);
keep_sig_mat_phase( cellfun(@isempty, keep_sig_mat_phase) ) = {0};
keep_sig_mat_phase = cell2mat(keep_sig_mat_phase);

color_brain_regions = util.get_color_rgb_codes(unit_regions);
colors_all = {'blue', 'green', 'magenta'};
title_names = {'Motor Imagery', 'Spoken Grasp', 'Spoken Colors'};

% plot figure
fig = figure('units','normalized','outerposition',[0 0 0.3 0.4]);
distances = [1,2,3];
 
 for n_reg_2 = 1:number_regions
    subplot(1,3,n_reg_2)
    action_data_tmp = cell2mat(data_phase(:,n_reg_2));
    action_data_shuffle = cell2mat(data_shuffle_phase(:,n_reg_2));

    b = bar(data_mean_phase(:,n_reg_2),'FaceColor','flat');
    hold on
    b.CData = cell2mat(color_brain_regions');
    hold on 
    
    for n_phase_2 = 1:number_regions
        plot(n_phase_2*ones(1,5),action_data_tmp(n_phase_2,:), 'ko','markerfacecolor','k', 'MarkerSize', 2)
        hold on  
        plot(n_phase_2*ones(1,5),action_data_shuffle(n_phase_2,:), 'ro','markerfacecolor','k', 'MarkerSize', 2)
        hold on  
    end 
    
    hold on
    ErrCi = cell2mat(err_ci_phase(:,n_reg_2)');
    er = errorbar(1:3, data_mean_phase(:,n_reg_2),ErrCi(1,:),ErrCi(2,:));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    ylim([0 110]);
    xticklabels(unit_regions)
    xtickangle(-45)
    if n_reg_2 == 1; ylabel('Classification accuracy (%)', 'FontSize', 15); end
    title(title_names{n_reg_2})
        
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',13) 
    if n_reg_2 ~= 1; set(gca,'YTickLabel',[]); end

    
   for n_sig = 1:length(keep_sig_mat_phase(:,n_reg_2))
       
       if keep_sig_mat_phase(n_sig,n_reg_2) == 3
          text(distances(n_sig) - 0.4,max(action_data_tmp(n_sig,:)) + 5 , '***', 'Color', colors_all{n_sig}, 'LineWidth', 2,'FontSize', 17); %, 'Interpretation', 'latex')

       elseif keep_sig_mat_phase(n_sig,n_reg_2) == 2
          text(distances(n_sig) - 0.3,max(action_data_tmp(n_sig,:)) + 5 , '**', 'Color', colors_all{n_sig}, 'LineWidth', 2,'FontSize',17); %, 'Interpretation', 'latex')

       elseif keep_sig_mat_phase(n_sig,n_reg_2) == 1
          text(distances(n_sig) - 0.15,max(action_data_tmp(n_sig,:)) + 5 , '*', 'Color', colors_all{n_sig}, 'LineWidth', 2,'FontSize', 17); %, 'Interpretation', 'latex')
       end 
   end 

end

sgtitle(['Decoding of different tasks during ' phase_names{phase_to_compute} ' phase']);

%% functions
%calculate confidence interval
function err_ci = calculate_err_ci(Data)
    NA = size(Data,2);
    SEM_A = std(Data') / sqrt(NA); % Standard Error Of The Mean
    CI95 = tinv([0.025 0.975], NA-1);   
    err_ci =  abs(bsxfun(@times, SEM_A, CI95(:)));
end 
 


