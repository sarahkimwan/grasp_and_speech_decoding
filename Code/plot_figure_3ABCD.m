%% plot figure 3 of the manuscript

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all 
close all

%% variables
%saved data folder
SavedData = [pwd '\Data\CrossPhaseNeuronDropping\MotorImagery\'];

%brain regions involved
unit_regions = {'SMG', 'PMV', 'S1'};

accuracy1_all = cell(3,2);
CI_all = cell(3,2);
training_cues = {'Cue', 'Action'};

%% compute figures 3A,B,C
for n_unit_region = 1:length(unit_regions)
    
    unit_region = unit_regions{n_unit_region};   
    %load precomputed data
    load([SavedData,'s2_' unit_region '_PC_20' ])
    fig = figure();

    for n_training_cue = 1:length(training_cues)
        subplot(1,2,n_training_cue)
        
        result = resultsPC{n_training_cue};
        data1name = CueNamesPC{n_training_cue}{1};
        data2name = CueNamesPC{n_training_cue}{2};

        acc1tmp = cellfun(@(x) squeeze(x(:,1,:)),result', 'UniformOutput', false);
        acc2tmp = cellfun(@(x) squeeze(x(:,2,:)),result', 'UniformOutput', false);

        acc1 = (1-cell2mat(cellfun(@(x) mean(x,1),acc1tmp, 'UniformOutput', false)))*100;
        acc2 = (1-cell2mat(cellfun(@(x) mean(x,1),acc2tmp, 'UniformOutput', false)))*100;
        %compute 95% bootstrapped confidence interval of the mean
        CI_acc1 = bootci(1e3,{@nanmean,acc1});
        CI_acc2 = bootci(1e3,{@nanmean,acc2});

        CI1 = [CI_acc1(2,:) - mean(acc1); mean(acc1) - CI_acc1(1,:)];
        CI2 = [CI_acc2(2,:) - mean(acc2); mean(acc2) - CI_acc2(1,:)];
        %save classification accuracy and confidence intervals for plot 3D
        accuracy1_all{n_unit_region,n_training_cue} = acc1;
        CI_all{n_unit_region,n_training_cue} = CI1;
        %colors for data in plot
        color_phases = util.get_color_rgb_codes({data1name, data2name});
        %plot error bars
        b1 = errorbar(1:size(acc1,2), mean(acc1), CI1(1,:), CI1(2,:),'-s', 'LineStyle', '-','MarkerSize',5,...
            'MarkerEdgeColor',color_phases{1},'MarkerFaceColor',color_phases{1} , 'Color', color_phases{1});
        hold on 

        b2 = errorbar(1:size(acc2,2), mean(acc2), CI2(1,:), CI2(2,:),'-s', 'LineStyle', '-','MarkerSize',5,...
           'MarkerEdgeColor',color_phases{2},'MarkerFaceColor',color_phases{2} , 'Color', color_phases{2});
        hold on 
        
        %figure parameters
        if strcmp(unit_region, 'SMG') 
            xticks([0:200:size(acc1,2)])
        end   
        if strcmp(unit_region, 'S1X')
            xticks([0:300:size(acc1,2)])   
        end 
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',15) 
        hold on
         chance_level = 1/5*100;
        l = line([0 size(acc1,2) ],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
        yticks([0:20:100])
        ylim([0 100])
        xlim([0 size(acc1,2)])
        title(['Training: ' data1name ' Phase'])

        %calculate average number of units needed to obtain > 80% classification accuracy
        Acc80 = find(mean(acc1)>80,1);
        disp(['Reach 80% accuracy training with ' data1name ' with ' num2str(Acc80) ' unit region ' unit_region])        
    end 
    sgtitle(unit_region)
end

legend([b1, b2], {'Cue phase', 'Action phase'})

%% figure 3D

brain_region_colors = util.get_color_rgb_codes({'SMG', 'PMV', 'S1X'});
fig_titles = {'Cue phase', 'Action phase'};

direction = {{'right', 'left'}, {'left','right'}}; %needed for plotting

figure();
for n_phase = 1:length(fig_titles)
    subplot(1,2,n_phase)

    %plot results for the first 146 units (total number of PMV units)
    acc1 = accuracy1_all{1,n_phase}(:,1:146);
    acc2 = accuracy1_all{2,n_phase};
    acc3 = accuracy1_all{3,n_phase}(:,1:146);

    CI1 = CI_all{1, n_phase}(:,1:146);
    C12 = CI_all{2,n_phase};
    C13 = CI_all{3,n_phase}(:,1:146);

    %plot error bars 
    b1= errorbar(1:size(acc1,2), mean(acc1), CI1(1,:), CI1(2,:),'-s', 'LineStyle', '-','MarkerSize',5,...
    'MarkerEdgeColor',brain_region_colors{1},'MarkerFaceColor',brain_region_colors{1} , 'Color', brain_region_colors{1});
    hold on 

    b2 =errorbar(1:size(acc2,2), mean(acc2), C12(1,:), C12(2,:),'-s', 'LineStyle', '-','MarkerSize',5,...
    'MarkerEdgeColor',brain_region_colors{2},'MarkerFaceColor',brain_region_colors{2} , 'Color', brain_region_colors{2});

    hold on 
    b3 = errorbar(1:size(acc3,2), mean(acc3), C13(1,:), C13(2,:),'-s', 'LineStyle', '-','MarkerSize',5,...
    'MarkerEdgeColor',brain_region_colors{3},'MarkerFaceColor',brain_region_colors{3} , 'Color', brain_region_colors{3}); 

    %figure parameters
    ylim([0 100])
    xlim([0 size(acc1,2)])
    hold on
    chance_level = 1/5*100;
    l = line([0 size(acc1,2) ],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
    yticks([0:20:100])
    xlabel('Number units')
    ylabel('Classification accuracy')

    %plot number of units needed to obtain 80% classification accuracy
    acc80_1 = find(mean(acc1)>80,1);
    acc80_2 = find(mean(acc2)>80,1);
    %figure parameters
    y1 = yline(80, 'k--', {'80 %'});
    y1.FontSize = 14;
    y1.LabelVerticalAlignment = 'top';
    y1.LabelHorizontalAlignment = 'left';
    hold on
    x1=xline(acc80_1, '--', {['N = ' num2str(acc80_1)]},'Color',brain_region_colors{1},'LineWidth',1.5);
    x1.LabelVerticalAlignment = 'middle';
    x1.LabelHorizontalAlignment = direction{1, n_phase}{1, 1};
    x1.FontSize =14;
    hold on
    x2=xline(acc80_2, '--', {['N = ' num2str(acc80_2)]},'Color',brain_region_colors{2},'LineWidth',1.5);
    x2.LabelVerticalAlignment = 'middle';
    x2.LabelHorizontalAlignment = direction{1, n_phase}{1, 2};
    x2.FontSize =14;
    title(fig_titles{n_phase})
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',15) 
    
end 

sgtitle('SMG, PMV and S1 decoding comparison')
legend([b1,b2, b3],unit_regions)
