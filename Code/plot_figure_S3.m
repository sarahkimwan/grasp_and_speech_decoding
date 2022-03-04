%% plot figure S3 of the manuscript

%% Important: run code while being in folder 'grasp_and_speech_decoding'
clc
clear all
close all

%% changeable parameters

flagGraspSpeaking = false; 
%if true  - computes figure for Spoken Grasps task (S3A)
%if false - computes figure for Spoken Colors task (S3B)

SavedData = [pwd '\Data\CrossPhaseNeuronDropping\Speech\'];

%%
training_cues = {'Cue', 'Action'};
unit_region = 'SMG';
%load data
if flagGraspSpeaking
    load(fullfile(SavedData, 's2_SMG_PC_20_GraspSpeaking.mat'))
else
    load(fullfile(SavedData, 's2_SMG_PC_20_ColorSpeaking.mat'))
end 

figure()
%plot figure
for n_training_cue = 1:length(training_cues)
    subplot(1,2,n_training_cue)
    result = resultsPC{n_training_cue};

    data1name = CueNamesPC{n_training_cue}{1};
    data2name = CueNamesPC{n_training_cue}{2};

    acc1tmp = cellfun(@(x) squeeze(x(:,1,:)),result', 'UniformOutput', false);
    acc2tmp = cellfun(@(x) squeeze(x(:,2,:)),result', 'UniformOutput', false);

    acc1 = (1-cell2mat(cellfun(@(x) mean(x,1),acc1tmp, 'UniformOutput', false)))*100;
    acc2 = (1-cell2mat(cellfun(@(x) mean(x,1),acc2tmp, 'UniformOutput', false)))*100;

    acc1A = (1-cell2mat(acc1tmp))*100;
    acc2A = (1-cell2mat(acc2tmp))*100;

    %calculate bootstrapped confidence intervals
    CI_acc1 = bootci(1e3,{@nanmean,acc1});
    CI_acc2 = bootci(1e3,{@nanmean,acc2});
    CI1 = [CI_acc1(2,:) - mean(acc1); mean(acc1) - CI_acc1(1,:)];
    CI2 = [CI_acc2(2,:) - mean(acc2); mean(acc2) - CI_acc2(1,:)];

    color_phases = util.get_color_rgb_codes({data1name, data2name});
    %plot error bars    
    b1 = errorbar(1:size(acc1,2), mean(acc1), CI1(1,:), CI1(2,:),'-s', 'LineStyle', '-','MarkerSize',5,...
        'MarkerEdgeColor',color_phases{1},'MarkerFaceColor',color_phases{1} , 'Color', color_phases{1}) ;
    hold on 
    b2 = errorbar(1:size(acc2,2), mean(acc2), CI2(1,:), CI2(2,:),'-s', 'LineStyle', '-','MarkerSize',5,...
        'MarkerEdgeColor',color_phases{2},'MarkerFaceColor',color_phases{2} , 'Color', color_phases{2}) ;
    hold on 

    %figure properties
    xticks([0:200:size(acc1,2)])
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',15) 
    hold on
    chance_level = 1/5*100;
    l = line([0 size(acc1,2) ],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
    yticks([0:20:100])
    acc80 = find(mean(acc1)>80,1);
    disp(['Reach 80% accuracy training with ' data1name ' with ' num2str(acc80) ' unit region ' unit_region])
    ylim([0 100])
    xlim([0 size(acc1,2)])
    title(['Training: ' data1name ' Phase'])
    ylabel('Classification accuracy [%]')
end 

plot_title_sup = { 'Spoken Colors', 'Spoken Grasps'};
sgtitle([plot_title_sup{flagGraspSpeaking+1}, ' cross phase decoding'])
legend([b1, b2], CueNamesPC{2})