%% plot figure 2 of the manuscript

%% Important: run code while being in folder 'grasp_and_speech_decoding'

clc
clear all 
close all

%% changeable parameters

flagGoData = true; 
%if true  - computes figure for Go trials 
%if false - computes figure for NoGo trials 

SavedData = [pwd '\Data\ClassificationMotorImagery\'];

%% load data

if flagGoData
    %Go data
    SMG = load([SavedData 'SMG_Error.mat']);
    PMV = load([SavedData 'PMV_Error.mat']);
    S1 = load([SavedData 'S1_Error.mat']);
    %Go shuffle data
    SMG_shuffle = load([SavedData  'SMG_ErrorShuffle.mat']);
    PMV_shuffle = load([SavedData  'PMV_ErrorShuffle.mat']);
    S1_shuffle = load([SavedData  'S1_ErrorShuffle.mat']);
    title_add = 'Go trials';
else
    %NoGo data
    SMG = load([SavedData  'SMG_Error_NoGo.mat']);
    PMV = load([SavedData  'PMV_Error_NoGo.mat']);
    S1 = load([SavedData  'S1_Error_NoGo.mat']);
    %NoGo shuffle Data
    SMG_shuffle = load([SavedData  'SMG_ErrorShuffle_NoGo.mat']);
    PMV_shuffle = load([SavedData  'PMV_ErrorShuffle_NoGo.mat']);
    S1_shuffle = load([SavedData  'S1_ErrorShuffle_NoGo.mat']);
    title_add = 'No-Go trials';
end 

%% figure 2A/B
%variables
number_phases = 4;
marker_size = 5;

%Mean classification accuracy values 
SMG.keepMeanAcc =  (1-squeeze(mean(squeeze(SMG.errTest),1)))*100;
%remove PMV and S1 sessions that have no data 
PMV.keepSessionIdxToRemove = [5,6,7];
PMVToKeep = setdiff(1:SMG.number_sessions,PMV.keepSessionIdxToRemove);
PMV.keepMeanAcc = (1-squeeze(mean(squeeze(PMV.errTest(:,:,:,PMVToKeep)),1)))*100;
S1.keepSessionIdxToRemove = [1,3];
S1ToKeep = setdiff(1:SMG.number_sessions,S1.keepSessionIdxToRemove);
S1.keepMeanAcc = (1-squeeze(mean(squeeze(S1.errTest(:,:,:,S1ToKeep)),1)))*100;
    
%Mean shuffled labels classification accuracy values
SMG_mean_acc_shuffle = (1-squeeze(mean(squeeze(mean(SMG_shuffle.errTest(:,:,:,:),1)),1)))*100;
PMV_mean_acc_shuffle = (1-squeeze(mean(squeeze(mean(PMV_shuffle.errTest(:,:,:,PMVToKeep),1)),1)))*100;
S1_mean_acc_shuffle = (1-squeeze(mean(squeeze(mean(S1_shuffle.errTest(:,:,:,S1ToKeep),1)),1)))*100;

%figure parameters
fig = figure('units','normalized','outerposition',[0 0 0.4 0.5]);
font_size = 13; 
marker_style = 's';
marker_size = 5;

%SMG data 

%plot individual session days classification and average classification
subplot(1,3,1)
plot(repmat([1; 2; 3; 4;], [1, length(SMG.keepMeanAcc)]), squeeze(SMG.keepMeanAcc), 'Marker', '.', 'MarkerSize', marker_size, 'LineStyle','none', 'Color', 'black');
hold on
plot(repmat([1; 2; 3; 4;], [1, length(SMG_mean_acc_shuffle)]), squeeze(SMG_mean_acc_shuffle), 'Marker', '.', 'MarkerSize', marker_size, 'LineStyle','none', 'Color', 'red');
hold on
plot([1; 2; 3; 4;], mean(SMG_mean_acc_shuffle,2), 'Marker','.', 'MarkerSize', marker_size, 'LineStyle','--', 'Color', 'red');
%add number of sessions on figure
n_sessions = size(SMG_mean_acc_shuffle,2); 
text(1 ,98 , ['N = ' num2str(n_sessions)], 'Color', 'black', 'LineWidth', 4,'FontSize', font_size); %, 'Interpretation', 'latex')
  
%Calculate and plot 95% Confidence interval    
err_ci = calculate_err_ci(SMG.keepMeanAcc); 
errorbar(1:number_phases, mean(SMG.keepMeanAcc,2)', err_ci(1,:), err_ci(2,:),'-s','Marker', marker_style, 'LineStyle', '-','MarkerSize',marker_size,'Color', 'blue',...
    'Color', 'blue','LineWidth',1);
hold on
d1 = plot(mean(SMG.keepMeanAcc,2)','Marker', marker_style, 'MarkerSize',marker_size, 'LineStyle','none', 'Color', 'blue','LineWidth',2);
hold on 
title('SMG', 'FontSize', font_size)

%Calculate percentile of shuffled dataset and compare them to averaged classification accuracy 
shuffle_data_SMG = (1 -reshape((mean(SMG_shuffle.errTest(:,:,:,:),2)), [],number_phases))*100;
sig_plot(SMG.keepMeanAcc, shuffle_data_SMG);


%plot parameters
figInfo()
ylabel('Classification accuracy [%]', 'FontSize', 13)
yticks([0:20:100])

%PMV data

%plot individual session days classification and average classification
subplot(1,3,2)
plot(repmat([1; 2; 3; 4;], [1, length(PMV.keepMeanAcc)]), squeeze(PMV.keepMeanAcc), 'Marker','.', 'MarkerSize', marker_size, 'LineStyle','none', 'Color', 'black');
hold on
plot(repmat([1; 2; 3; 4;], [1, length(PMV_mean_acc_shuffle)]), squeeze(PMV_mean_acc_shuffle), 'Marker', '.', 'MarkerSize', marker_size, 'LineStyle','none', 'Color', 'red');
hold on
plot([1; 2; 3; 4;], mean(PMV_mean_acc_shuffle,2), 'Marker','.', 'MarkerSize', marker_size, 'LineStyle','--', 'Color', 'red');
title('PMV','FontSize', font_size)
set(gca,'YTickLabel',[]);
%add number of sessions on figure
n_sessions = size(PMV_mean_acc_shuffle,2); 
text(1 ,98 , ['N = ' num2str(n_sessions)], 'Color', 'black', 'LineWidth', 4,'FontSize', font_size); %, 'Interpretation', 'latex')

%Calculate and plot 95% Confidence interval 
err_ci = calculate_err_ci(PMV.keepMeanAcc); 
errorbar(1:number_phases, mean(PMV.keepMeanAcc,2)', err_ci(1,:), err_ci(2,:),'-s','Marker', marker_style, 'LineStyle', '-','MarkerSize',marker_size,'Color', 'blue',...
'Color', 'green','LineWidth',1);
hold on
%plot average marker
d2 = plot(mean(PMV.keepMeanAcc,2)','Marker', marker_style, 'MarkerSize', marker_size, 'LineStyle','none', 'Color', 'green','LineWidth',2);
hold on 

%Calculate percentile of shuffled dataset and compare them to averaged classification accuracy 
shuffle_data_PMV = (1 -reshape((mean(PMV_shuffle.errTest(:,:,:,PMVToKeep),2)), [],number_phases))*100;
sig_plot(PMV.keepMeanAcc, shuffle_data_PMV);

figInfo();

%S1 data%
subplot(1,3,3)
d = plot(repmat([1; 2; 3; 4;], [1, length(S1.keepMeanAcc)]), squeeze(S1.keepMeanAcc), 'Marker', '.', 'MarkerSize', marker_size, 'LineStyle','none', 'Color', 'black');
hold on
s = plot(repmat([1; 2; 3; 4;], [1, length(S1_mean_acc_shuffle)]), squeeze(S1_mean_acc_shuffle), 'Marker', '.', 'MarkerSize', marker_size, 'LineStyle','none', 'Color', 'red');
hold on
s = plot([1; 2; 3; 4;], mean(S1_mean_acc_shuffle,2), 'Marker', '.', 'MarkerSize', marker_size, 'LineStyle','--', 'Color', 'red');
title('S1','FontSize', font_size)
set(gca,'YTickLabel',[]);
%add number of sessions on figure
n_sessions = size(S1_mean_acc_shuffle,2); 
text(1 ,98 , ['N = ' num2str(n_sessions)], 'Color', 'black', 'LineWidth', 4,'FontSize', font_size); 

%Calculate and plot 95% Confidence interval 
err_ci = calculate_err_ci(S1.keepMeanAcc); 
errorbar(1:number_phases, mean(S1.keepMeanAcc,2)', err_ci(1,:), err_ci(2,:),'-s','Marker', 's', 'LineStyle', '--','MarkerSize',5,'Color', 'blue',...
'Color', 'magenta','LineWidth',1);
hold on
% plot average marker
d3 = plot(mean(S1.keepMeanAcc,2)','Marker', 's', 'MarkerSize', 5, 'LineStyle','none', 'Color', 'magenta','LineWidth',2);
hold on 

%Calculate percentile of shuffled dataset and compare them to averaged classification accuracy 
shuffle_data_S1 = (1 -reshape((mean(S1_shuffle.errTest(:,:,:,S1ToKeep),2)), [],number_phases))*100;
sig_plot(S1.keepMeanAcc, shuffle_data_S1);

figInfo()

%title of general plot 
sgtitle(title_add, 'FontSize', 15);

%% figure 2C/D

%SMG Error matrix 
SMG_test_labels  = squeeze(SMG.keepTestingLabels);
SMG_predicted_labels  = squeeze(SMG.keepPredictedLabels);
plot_error_matrix(SMG_test_labels, SMG_predicted_labels, 'SMG', title_add);

%PMV error matrix       
PMV_test_labels  = squeeze(PMV.keepTestingLabels);
PMV_predicted_labels  = squeeze(PMV.keepPredictedLabels);
plot_error_matrix(PMV_test_labels, PMV_predicted_labels, 'PMV', title_add);

%% 
%there is no significant difference in classification accuracy
%between Go trials from only Go blocks, or from Go trials from GoNoGo
%blocks.

GoNoGo_sessions = SMG.Sessions(1:9);
Go_sessions = SMG.Sessions(10:end);

GoNoGo_means = SMG.keepMeanAcc(:,1:9);
Go_means = SMG.keepMeanAcc(:,10:end);

[SMG_hC,SMG_pC] = ttest2(GoNoGo_means(2,:),Go_means(2,:));
[SMG_hD,SMG_pD] = ttest2(GoNoGo_means(3,:),Go_means(3,:));
[SMG_hA,SMG_pA]= ttest2(GoNoGo_means(4,:),Go_means(4,:));

PMV_NoGo_idx = ismember(PMV.Sessions(PMVToKeep),GoNoGo_sessions);
PMV_Go_idx = ismember(PMV.Sessions(PMVToKeep),Go_sessions);

PMV_NoGo = PMV.keepMeanAcc(:,PMV_NoGo_idx);
PMV_Go = PMV.keepMeanAcc(:,PMV_Go_idx);

[PMV_hC,PMV_pC] = ttest(PMV_NoGo(2,:),PMV_Go(2,:));
[PMV_hD,PMV_pD] = ttest(PMV_NoGo(3,:),PMV_Go(3,:));
[PMV_hA,PMV_pA] = ttest(PMV_NoGo(4,:),PMV_Go(4,:));

S1_NoGo_idx = ismember(S1.Sessions(S1ToKeep), GoNoGo_sessions);
S1_Go_idx = ismember(S1.Sessions(S1ToKeep), Go_sessions);

S1_NoGo = S1.keepMeanAcc(:,S1_NoGo_idx);
S1_Go = S1.keepMeanAcc(:,S1_Go_idx);

[S1_hC,S1_pC] = ttest2(S1_NoGo(2,:),S1_Go(2,:));
[S1_hD,S1_pD] = ttest2(S1_NoGo(3,:),S1_Go(3,:));
[S1_hA,S1_pA] = ttest2(S1_NoGo(4,:),S1_Go(4,:));

%% functions 

%calculate and plot significant classification results 
function sig_plot(mean_classification, shuffled_classification)
    number_phases = 4;
    font_size = 16;
    %values for percentile tests
    prc_tests = [95,99,99.9];    significant_values = zeros(nnz(prc_tests), number_phases); %predeclare variables
    for n_perc = 1:length(prc_tests)  
        significant_values(n_perc,:) = prctile(shuffled_classification,prc_tests(n_perc)) < mean(mean_classification,2)';
    end 
    %plot significant values 
    max_values = max(mean_classification');
    
    for n_phases = 1:number_phases
       sig_values_per_phase = significant_values(:,n_phases);
        if sig_values_per_phase(3)
            text(n_phases - 0.5,max_values(n_phases) + 5 , '***', 'Color', 'black', 'LineWidth', 4, 'FontSize', font_size); 
        elseif sig_values_per_phase(2)
            text(n_phases - 0.4,max_values(n_phases) + 5 , '**', 'Color', 'black', 'LineWidth', 4, 'FontSize', font_size); 
        elseif sig_values_per_phase(1)
            text(n_phases - 0.3,max_values(n_phases) + 5 , '*', 'Color', 'black','LineWidth', 4, 'FontSize', font_size); 
        end 
    end 
end 

function err_ci = calculate_err_ci(Data)
    NA = size(Data,2);
    SEM_A = std(Data') / sqrt(NA); % Standard Error Of The Mean
    CI95 = tinv([0.025 0.975], NA-1);   
    err_ci =  abs(bsxfun(@times, SEM_A, CI95(:)));
end 

%plot figure labels, fixed xticks and yticks
function figInfo()

    xticks(1:4)
    xticklabels({'ITI', 'Cue', 'Delay', 'Action'})
    xtickangle(-45)
    ylim([0 105])
    yticks(0:10:100)
    xlim([0 5])
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',15) 
end 

% plot error matrix 
function plot_error_matrix(test_labels, predicted_labels, fig_title, title_add)

    number_phases = 4;
    test_labels = reshape(cell2mat(test_labels(number_phases,:)), [],1);
    PredictLabels = reshape(cell2mat(predicted_labels(number_phases,:)), [],1);
    grasps_to_test = arrayfun(@(x) util.image2class_simple(x), unique(test_labels), 'UniformOutput', false);
    number_grasps = length(grasps_to_test); 
    c_mat = confusionmat(test_labels, PredictLabels);
    number_repetition_per_grasp = unique(histc(test_labels, unique(test_labels)));
    c_mat = c_mat/number_repetition_per_grasp*100;

    figure()
    colormap(jet)

    imagesc(c_mat)
    xticks(1:number_grasps)
    xticklabels(grasps_to_test)
    xtickangle(45)
    yticks(1:number_grasps)
    yticklabels(grasps_to_test)
    xlabel('Predicted Class')
    ylabel('True Class')
    caxis([0 100])
    colorbar

    title([fig_title ' ' title_add ' - Action phase'])

end 