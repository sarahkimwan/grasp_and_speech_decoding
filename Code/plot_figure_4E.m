%% plot cross-task classification figure 4E of manuscript

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all
close all

SavedData = [pwd '\Data\CrossModalityClassification\'];
%load Leave-One-Out results from Cross Modality Classification
Data = load([SavedData 'RealData']);
%load Leave-One-Out shuffled results (1000x) from Cross Modality Classification
Shuffle = load([SavedData 'ShuffledData']);

%% plot figure

task_names =  {'Motor Imagery','Spoken Grasps','Spoken Colors'};
number_tasks = 3;
number_phases = 4;
number_sessions = 5;
%percentages employed to compute significance 
prc_tests = [95,99,99.9];

%predeclare variables
data_mean_all = cell(1,number_tasks);
shuffle_mean_all = cell(1,number_tasks);
prc_out_all = cell(1,number_tasks);

%loop through tasks - the classification model was trained one on each
%task, and each time evalauted on all of them
for n_mod = 1:length(task_names)
    
    %predeclare temporary variable
    prc_outcome = nan(number_tasks, number_tasks, number_phases);
    task_names_per_classification = Data.DataNamesAll{n_mod}{1};
    
    %real data
    data_err_tmp = Data.errorsAll{n_mod};
    %shuffled data
    shuffle_err_tmp = Shuffle.errorsAll(:,n_mod);
   
    %calculate error over leave-one-out folds 
    data_mean_per_session = cellfun(@(x) 1 -squeeze(mean(x,2)), data_err_tmp, 'UniformOutput', false);
    %calculate mean error over all sessions
    data_mean_all_sessions = cell2mat(cellfun(@(x) mean(x,1), data_mean_per_session,'UniformOutput', false)');   
    shuffle_err = cell(1,number_tasks);
    
    %compute percentile value for each task and each percentile value
    for n_tasks = 1:length(task_names_per_classification)        
        shuffle_err{n_tasks} = cell2mat(cellfun(@(x) 1 -squeeze(mean(x{n_tasks},2)), shuffle_err_tmp, 'UniformOutput', false));             
        for n_prc = 1:length(prc_tests)
            prcV = prctile(shuffle_err{n_tasks},prc_tests(n_prc));
            prc_outcome(n_tasks,n_prc,:) = data_mean_all_sessions(n_tasks,:) >  prcV;
        end        
    end 
   
    %save classification and percentile results  
    data_mean_all{n_mod} = data_mean_per_session;
    shuffle_mean_all{n_mod} = cell2mat(cellfun(@(x) mean(x,1), shuffle_err, 'UniformOutput', false)');
    prc_out_all{n_mod} = prc_outcome;
end 

%% plot figure

figure();
phase_names = {'ITI', 'Cue', 'Delay', 'Action'};

for n_modality = 1:length(task_names)
   
x = 1:length(phase_names);

labels_per_task = Data.DataNamesAll{n_modality}{1};

colors_per_task = util.get_color_rgb_codes(labels_per_task);

err_test{1} = data_mean_all{n_modality}{1}*100;
err_test{2} = data_mean_all{n_modality}{2}*100;
err_test{3} = data_mean_all{n_modality}{3}*100;

subplot(1,3,n_modality)
offset = [-0.04, -0.02,0, 0.02, 0.04];
a = cell(1,number_tasks);
    for err = 1:length(err_test)
        err_test1 = err_test{err};
        for i = 1:number_sessions
            plot((1:number_phases) + offset(i),err_test1(i,:), '.', 'MarkerSize',13, 'Color', colors_per_task{err})     
            hold on
        end
        %95 Percent Confidence interval 
        NA = size(err_test1',2);
        SEM_A = std(err_test1) / sqrt(NA); % Standard Error Of The Mean
        CI95 = tinv([0.025 0.975], NA-1);   
        err_ci =  abs(bsxfun(@times, SEM_A, CI95(:)));

        a{err}= errorbar(1:number_phases, mean(err_test1,1), err_ci(1,:), err_ci(2,:),'-s','Marker', 's', 'LineStyle', '-','MarkerSize',10,'Color', colors_per_task{err},'LineWidth',1);
        hold on
        if err == 1
            y_data = max(err_test1);
        end 
        %plot significance 
        sig_values = squeeze(prc_out_all{n_modality}(err,:,:));
        for n_p = 1:length(phase_names)
           sig_val_per_phase = sig_values(:,n_p);
            if sig_val_per_phase(3)
                text(n_p - 0.3,y_data(n_p) + 10 - 2*err , '***', 'Color', colors_per_task{err}, 'FontSize', 15); 
            elseif sig_val_per_phase(2)
                text(n_p - 0.2,y_data(n_p) + 10 - 2*err , '**', 'Color', colors_per_task{err}, 'FontSize', 15); 
            elseif sig_val_per_phase(1)
                text(n_p - 0.1,y_data(n_p) + 10 - 2*err , '*', 'Color', colors_per_task{err}, 'FontSize', 15); 
            end 

        end 
    end
    
    %plot mean shuffled results
    hold on;
    for kk = 1:length(err_test)
        plot(shuffle_mean_all{n_mod}(kk,:)*100,'r', 'Marker', '.', 'MarkerSize', 15,'LineStyle', '--')
    end 
    
    title(['Train: ' labels_per_task{1} ])

    set(gca,'FontSize', 12)
    ylim([0 110])
    xlim([0.5, 4.5])
    xticks(1:length(phase_names));
    xticklabels(phase_names);
    xtickangle(-45)
    if n_modality == 1
        ylabel('Classification accuracy')
    end
    
    legend([a{1},a{2},a{3}],labels_per_task)
end
sgtitle('Cross Modality decoding')