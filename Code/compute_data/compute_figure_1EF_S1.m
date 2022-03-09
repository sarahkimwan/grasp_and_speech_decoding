%% compute tuning per analysis cue and actiion phase window (figure 1EF)

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
flag50msTuning = false; 
%if true - computes tuning in 50ms time bins
%if false - only computes tuning for average firing rate per phase

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
sumNbrTunedChannelsPerGrasp = cell(number_tasks,number_sessions);
sumNbrTunedToSeveralGrasps = cell(number_tasks,number_sessions);
ChannelCount = zeros(number_tasks, number_sessions);
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
    Grasp_Names = util.label_number_to_grasp_name(unique(labels)); 
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
        flagSelectTimeBins = true; 
        if flagSelectTimeBins 
            time_phase_labels_old = time_phase_labels{1};
            start_phase_time = 0.05;
            start_phase_time(2:4) = (find(diff(time_phase_labels_old))+1)*0.05;
            time_phase_labels_adapted = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;...
                            2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;...
                            0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;...
                            0;0;0;0;0;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;0;0;0;0;0;0;0;0;0;0;0;0;...
                            0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

            tmp = cell(size(time_phase_labels));
            [tmp{:}] = deal(time_phase_labels_adapted);
            time_phase_labels = tmp;
        end
    
        %compute tuned units
        [~, tuned_units_per_phase,~,sumPhase, ~, nbr_tuned_units_per_grasp,nbr_tuned_units_per_grasp_bin] = analysis.get_tuned_units(SessionData,labels_per_session,time_phase_labels, flag50msTuning);         

        tuned_channels_per_phase_all{n_task,n_session} = tuned_units_per_phase;   
       
        sumNbrTunedChannelsPerGrasp{n_task, n_session} = nbr_tuned_units_per_grasp;
        sumNbrTunedToSeveralGrasps{n_task, n_session} =  nbr_tuned_units_per_grasp_bin;
        ChannelCountTmp = size(SessionData{1},2);
        ChannelCount(n_task,n_session) = ChannelCountTmp;

    end 
end

if flagSaveData
    if flagGoData
        save([data_save_folder 'GoData2'],'ChannelCount', 'GraspNames', 'Sessions', 'sumNbrTunedChannelsPerGrasp','sumNbrTunedToSeveralGrasps' )
    else
        save([data_save_folder 'NoGoData2'],'ChannelCount', 'GraspNames', 'Sessions', 'sumNbrTunedChannelsPerGrasp','sumNbrTunedToSeveralGrasps' )
    end  
end 

