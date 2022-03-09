%% compute tuning in 50ms time bins for motor imagery (figure 1C,D)

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all 
close all

%% changeable parameters

flagGoData = true; 
%if true  - computes tuning for Go trials 
%if false - computes tuning for NoGo trials - only exists for motor imagery data, not for spoken data
flag50msTuning = true; 
%if true - computes tuning in 50ms time bins
%if false - only computes tuning for average firing rate per task phase

flagSaveData = false; 
%if true - saves the data in save folder 
data_save_folder = [pwd '\Data\UnitTuning50msBins\']; %folder where computed data gets saved
%%   
if flagGoData
    %sessions with Go trials
    sessions_all = { {'20171003','20190404','20190409',    '20190417',    '20190423',...
    '20190509',    '20190510',    '20190522',    '20190528', '20190618', '20190911', '20190923', '20190924', '20190930','20191016'},...
    {'20171003','20190404',  '20190409',  '20190417', '20190522', '20190528', '20190618', ...
        '20190911', '20190923', '20190924', '20190930','20191016'}, ...
        {'20190404', '20190417',    '20190423',...
    '20190509',    '20190510',    '20190522',    '20190528', '20190618', '20190911', '20190923', '20190924', '20190930','20191016'}};
else
    %sessions with NoGo trials
    sessions_all = {{'20171003','20190404',   '20190409', '20190417',    '20190423',...
    '20190509',    '20190510',    '20190522',    '20190528'}, ...
    {'20171003','20190404','20190409',    '20190417','20190522','20190528'},...
    {'20190404',  '20190417',    '20190423',...
    '20190509',    '20190510',    '20190522',    '20190528'}};
end
%load data
data = load(fullfile(pwd, 'Data/NeuralData/MotorImagery', 'Table_sorting_aligned_thr_-4.5_ImageCue_perTrials.mat'));
task_names = {'Motor Imagery', 'Motor Imagery', 'Motor Imagery'};
unit_regions = {'SMG', 'PMV', 'S1'};


%variables
phase_names = {'ITI', 'Cue', 'Delay','Action'};

number_tasks = length(task_names);
number_phases = length(phase_names);
number_sessions = length(sessions_all{1});
    
tuned_channels_per_phase_all =cell(number_tasks,number_sessions);
sum_bin_all = cell(number_tasks, number_sessions);

%loop through number of brain regions 
for n_task = 1:length(unit_regions) 
    
    unit_region = unit_regions{n_task};
    sessions = sessions_all{n_task};
    number_sessions = length(sessions);

    %use Go or NoGo Data (only for motor imagery task)
    if flagGoData
        Go_data = data.Go_data;
    else
        Go_data = data.NoGo_data;
    end 

    labels = Go_data.GoLabels;
    session_idx = find(ismember(sessions_all{1}, sessions));

    for n_session = [session_idx]
        
        disp(['Tuning session ' sessions_all{1}{n_session} ]);       
        idxThisSession = ismember(Go_data.session_date, sessions_all{1}(n_session));
        
        if strcmp('SMG', unit_region)
            SessionData = Go_data.SMG_Go(idxThisSession,:);
        elseif strcmp('PMV', unit_region) 
            SessionData = Go_data.PMV_Go(idxThisSession,:);
        elseif strcmp('S1', unit_region)
            SessionData = Go_data.S1X_Go(idxThisSession,:);
        else
            error([unit_region ' does not exist '])
        end

        labels_per_session = labels(idxThisSession);
        time_phase_labels = Go_data.time_phase_labels(idxThisSession);
        %compute tuned units
        [~, tuned_units_per_phase,~,~,sumBin] = analysis.get_tuned_units(SessionData,labels_per_session,time_phase_labels, flag50msTuning);         
        sum_bin_all{n_task, n_session} = sumBin';
        tuned_channels_per_phase_all{n_task,n_session} = tuned_units_per_phase;        
    end 
end

% calculate number of total units 
num_channels = cell2mat(cellfun(@(x) size(x,1),tuned_channels_per_phase_all, 'UniformOutput', false));
num_channels_total = sum(num_channels');
% calculate percentage of tuned units
data1 = sum(cell2mat(sum_bin_all(1,:)'))/sum(num_channels(1,:))*100;
data2 = sum(cell2mat(sum_bin_all(2,:)'))/sum(num_channels(2,:))*100;
data3 = sum(cell2mat(sum_bin_all(3,:)'))/sum(num_channels(3,:))*100;


if flagSaveData
    data_save_name = {'NoGoData', 'GoData'};
    save([data_save_folder data_save_name{flagGoData + 1} '_MotorImagery'], 'data1', 'data2', 'data3', 'sum_bin_all', 'num_channels', 'num_channels_total'); 
end 

