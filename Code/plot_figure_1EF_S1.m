%% plot figure 1 E,F, supplemetary figure 1 A,B

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all
close all

%% changeable parameters

flagGoData = true; 
%if true  - computes figures for Go trials 
%if false - computes figures for NoGo trials 

phases_to_compute = [1,2,4];
% determines which task phases are plotted
% 1 = ITI, 2 = Cue, 3 = Delay, 4 = Action. 

saved_data = [pwd '\Data\GraspTuning'];
%% variables

if flagGoData
    %Go data
    data_Go = load([saved_data '\GoData.mat']);
    shuffle_all = load([saved_data '\GoDataShuffle.mat']);
else    
    %NoGo data
    data_Go = load([saved_data '\NoGoData.mat']);
    shuffle_all = load([saved_data '\NoGoDataShuffle.mat']);
end 

data_names = {'SMG', 'PMV', 'S1X'};
grasp_names = data_Go.GraspNames;
phase_names = {'ITI', 'Cue', 'Delay', 'Action'};

number_phases = length(phase_names);
number_grasps = length(grasp_names);
number_regions = length(data_names); 
number_sessions = size(data_Go.sumNbrTunedChannelsPerGrasp,2);
number_shuffle_repetitions = size(shuffle_all.sumNbrTunedChannelsPerGrasp{1},2);
number_phases_to_compute = length(phases_to_compute); 

colors_grasps = util.get_color_rgb_codes(grasp_names);
channel_count = data_Go.ChannelCount;

%predeclare variables
prc_values = cell(1,3); 
data_in_percentage = cell(1,3);
shuffle_accuracy = cell(1,3);
data_accuracy = cell(1,3); 
sum_tuning = zeros(number_regions,number_phases_to_compute);

%% plot figure

%prepare data to plot
%prepare shuffled data
for n_region = 1:number_regions 

    data_per_cue = shuffle_all.sumNbrTunedChannelsPerGrasp(n_region,:);
    summed_data_per_rep = zeros(number_phases_to_compute, number_grasps, number_shuffle_repetitions);
    summed_data_per_rep_all_grasps = zeros(number_phases_to_compute, number_shuffle_repetitions);
        
    for n_rep = 1:number_shuffle_repetitions 
        summed_data = zeros(length(phases_to_compute), number_grasps);
        for n_session = 1:number_sessions
            try % if a session is empty, ignore it 
                summed_data = data_per_cue{n_session}{n_rep}(phases_to_compute,:) + summed_data;
            catch
                continue
            end 
        end 
        
        summed_data_per_rep(:,:,n_rep) = summed_data;
        summed_data_per_rep_all_grasps(:,n_rep) = sum(summed_data,2);
    end 

    summed_data = mean(summed_data_per_rep,3);
    data_percentage_tmp =  zeros(number_phases_to_compute, number_shuffle_repetitions);

    for n_phase = 1:length(phases_to_compute )
        data_percentage_tmp(n_phase,:) = summed_data_per_rep_all_grasps(n_phase,:)/sum(channel_count(n_region,:))*100';
    end   

    data_in_percentage{n_region} = data_percentage_tmp;

    %compute percentile based on sum of tuned units of shuffled data for 95%, 99% and 99.9%
    prc_values{n_region}(1,:) = prctile(data_percentage_tmp',95); 
    prc_values{n_region}(2,:) = prctile(data_percentage_tmp',99); 
    prc_values{n_region}(3,:) = prctile(data_percentage_tmp',99.9); 

    shuffle_accuracy{n_region} = summed_data/sum(channel_count(n_region,:))*100;

end

%prepare actual data
for n_region = 1:number_regions
    data_per_cue = data_Go.sumNbrTunedChannelsPerGrasp(n_region,:);
    summed_data = zeros(number_phases_to_compute, number_grasps);

    for n_session = 1:number_sessions
        try % if a session is empty, ignore it 
            summed_data = data_per_cue{n_session}(phases_to_compute,:) + summed_data;
        catch
            continue
        end 
    end 
   
    sum_tuning(n_region,:) = sum(summed_data/sum(channel_count(n_region,:)),2)*100;
    data_accuracy{n_region} = summed_data/sum(channel_count(n_region,:))*100;
end

%% plot figure 2C/D

figure('units','normalized','outerposition',[0 0 1*0.4 0.35]);
for n_region = 1:number_regions
    subplot(1,number_regions,n_region)
    data_combined = zeros(number_regions,2,5);
    data_combined(:,1,:) = shuffle_accuracy{n_region};
    data_combined(:,2,:) = data_accuracy{n_region};
    util.plotBarStackGroups(data_combined, [1 2 3], colors_grasps)

    xticklabels(phase_names(phases_to_compute))
    xtickangle(-45)
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',13) 
    ylim([0 110])
    title(data_names{n_region})
    ylabel('Percentage of units (stacked)', 'FontSize', 15)

     
    %compute if sum of tuned data is significantly different from sum of
    %shuffled data
    significant_values = sum_tuning(n_region,:) > prc_values{n_region};
    hold on 
    y_data = sum_tuning(n_region,:);

    for n_phase_to_compute = 1:length(phases_to_compute)
       sig_values_per_phase = significant_values(:,n_phase_to_compute);
        if sig_values_per_phase(3)
            text(n_phase_to_compute - 0.15,y_data(n_phase_to_compute) + 7 , '***', 'Color', 'black', 'LineWidth', 2); 
            util.sigline([n_phase_to_compute-0.1625,n_phase_to_compute+0.1625],'', [], y_data(n_phase_to_compute));
        elseif sig_values_per_phase(2)
            text(n_phase_to_compute - 0.15,y_data(n_phase_to_compute) + 7 , '**', 'Color', 'black', 'LineWidth', 2); 
            util.sigline([n_phase_to_compute-0.1625,n_phase_to_compute+0.1625],'', [], y_data(n_phase_to_compute));
        elseif sig_values_per_phase(1)
            text(n_phase_to_compute - 0.05,y_data(n_phase_to_compute) + 7 , '*', 'Color', 'black', 'LineWidth', 2); 
            util.sigline([n_phase_to_compute-0.1625,n_phase_to_compute+0.1625],'', [], y_data(n_phase_to_compute));
        end 
    end 
end

legend_name = {'Lateral', 'Shuffle Lateral', 'WritingTripod', 'Shuffle WritingTripod', 'MediumWrap', 'Shuffle MediumWrap', 'PalmarPinch', 'Shuffle PalmarPinch', 'Sphere3Finger', 'Shuffle Sphere3Finger' };
legend(legend_name);

%% prepare figure S1

% prepare shuffle data
for n_region = 1:number_regions    

    data_per_cue = shuffle_all.sumNbrTunedToSeveralGrasps(n_region,:);
    summed_data_per_rep = zeros(number_phases_to_compute, number_grasps, number_shuffle_repetitions);
    summed_data_per_rep_all_grasps = zeros(number_phases_to_compute, number_shuffle_repetitions);

    for n_rep = 1:number_shuffle_repetitions
        summed_data = zeros(length(phases_to_compute), number_grasps+1);
        for n_session = 1:number_sessions

            try
                summed_data = data_per_cue{n_session}{n_rep}(phases_to_compute,:) + summed_data;
            catch
                %When a session is empty, skip it
            end 
        summed_data_per_rep(:,:,n_rep) = summed_data(:,2:end);
        summed_data_per_rep_all_grasps(:,n_rep) = sum(summed_data_per_rep(:,:,n_rep),2);
        end 
    end 
    
    summed_data = mean(summed_data_per_rep,3);
    data_percentage_tmp = zeros(number_phases_to_compute, number_shuffle_repetitions);
    
    for n_phase = 1:length(phases_to_compute) 
        data_percentage_tmp(n_phase,:) = summed_data_per_rep_all_grasps(n_phase,:)/sum(channel_count(n_region,:))*100';
    end   

    data_in_percentage{n_region} = data_percentage_tmp;
    %compute percentile based on sum of tuned units of shuffled data for 95%, 99% and 99.9%
    prc_values{n_region}(1,:) = prctile(data_percentage_tmp',95); 
    prc_values{n_region}(2,:) = prctile(data_percentage_tmp',99); 
    prc_values{n_region}(3,:) = prctile(data_percentage_tmp',99.9); 
   
    shuffle_accuracy{n_region} = summed_data/sum(channel_count(n_region,:))*100;
end

% prepares actual data
for n_region = 1:number_regions    
    data_per_cue = data_Go.sumNbrTunedToSeveralGrasps(n_region,:);
    summed_data = zeros(length(phases_to_compute), number_grasps+1);
    for n_session = 1:number_sessions
       try
        summed_data = data_per_cue{n_session}(phases_to_compute,:) + summed_data;
        catch
            %When a session is empty, skip it
            continue
        end 
    end 

    sum_tuning(n_region,:) = sum(summed_data(:,2:end),2)/sum(channel_count(n_region,:))*100;
    data_accuracy{n_region} = summed_data(:,2:end)/sum(channel_count(n_region,:))*100;
end

%% plot figure S1
color_number_units_tuned_to = {'#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499'};
figure('units','normalized','outerposition',[0 0 1*0.4 0.35]);

for n_region = 1:number_regions
    subplot(1,number_regions,n_region)
    data_combined = zeros(number_regions,2,5);
    data_combined(:,1,:) = shuffle_accuracy{n_region};
    data_combined(:,2,:) = data_accuracy{n_region};
    util.plotBarStackGroups(data_combined, [1 2 3],color_number_units_tuned_to)

    xticklabels(phase_names(phases_to_compute))
    xtickangle(-45)
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',13) 
    ylim([0 70])
        title(data_names{n_region})
    
    ylabel('Percentage of units (stacked)', 'FontSize', 15)    
     
    %Significance 
    significant_values = sum_tuning(n_region,:) > prc_values{n_region};
    hold on 
    y_data = sum_tuning(n_region,:);

    for n_phase_to_compute = 1:length(phases_to_compute)
       sig_values_per_phase = significant_values(:,n_phase_to_compute);
        if sig_values_per_phase(3)
            text(n_phase_to_compute - 0.15,y_data(n_phase_to_compute) + 7 , '***', 'Color', 'black', 'LineWidth', 2); %, 'Interpretation', 'latex')
            util.sigline([n_phase_to_compute-0.1625,n_phase_to_compute+0.1625],'', [], y_data(n_phase_to_compute));
        elseif sig_values_per_phase(2)
            text(n_phase_to_compute - 0.15,y_data(n_phase_to_compute) + 7 , '**', 'Color', 'black', 'LineWidth', 2); %, 'Interpretation', 'latex')
            util.sigline([n_phase_to_compute-0.1625,n_phase_to_compute+0.1625],'', [], y_data(n_phase_to_compute));
        elseif sig_values_per_phase(1)
            text(n_phase_to_compute - 0.05,y_data(n_phase_to_compute) + 7 , '*', 'Color', 'black', 'LineWidth', 2); %, 'Interpretation', 'latex')
            util.sigline([n_phase_to_compute-0.1625,n_phase_to_compute+0.1625],'', [], y_data(n_phase_to_compute));
        end 
    end 
end

legend_names = {'1 Grasp','1 Grasp Shuffle', '2 Grasps', '2 Grasp Shuffle','3 Grasps', '3 Grasp Shuffle','4 Grasps', '4 Grasp Shuffle','5 Grasps','5 Grasp Shuffle'};
legend(legend_names);

