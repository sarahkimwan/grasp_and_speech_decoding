%% compute tuning in 50ms time bins for speaking data (figure 4D)

%% Important: run code while being in folder 'grasp_and_speech_decoding'

clc
clear all 
close all

%% changeable parameters

flag50msTuning = true; 
%if true - computes tuning in 50ms time bins
%if false - only computes tuning for average firing rate per phase

flagSaveData = false; 
%if true - saves the data in save folder 
DataSaveFolder = [pwd '\Data\UnitTuning50msBins\']; %folder where computed data gets saved
%%
   
%sessions for Spoken Grasps and Spoken Colors
sessions_all = {{'20190911','20190923','20190924','20190930','20191016'},...
    {'20190911','20190923','20190924','20190930','20191016'}, ...
    {'20190911','20190923','20190924','20190930','20191016'}}; %Grasp Speaking All;
task_names = {'Grasps', 'Colors'};
unit_regions = {'SMG', 'SMG'};

%variables
phase_names = {'ITI', 'Cue', 'Delay','Action'};

number_tasks = length(task_names);
number_phases = length(phase_names);
number_sessions = length(sessions_all{1});
    
tuned_channels_per_phase_all =cell(number_tasks,number_sessions);
sum_bin_all = cell(number_tasks, number_sessions);

for n_task = 1:number_tasks 
    
    unit_region = unit_regions{n_task};
    data_name = task_names{n_task};     
    sessions = sessions_all{n_task};
    number_sessions = length(sessions);

    %load speaking SMG data 
    data = load(fullfile(pwd, 'Data/NeuralData/SpeechDataset2', ['Table_sorting_aligned_thr_-4.5_Speaking_' data_name 'PerModality.mat']));
    %extract go trial data
    Go_data = data.Go_data;
    %extract labels
    labels = Go_data.GoLabels;
    grasp_names = util.label_number_to_grasp_name(unique(labels)); 
    %find idx of current session
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
%calculate percentage of tuned units
data1 = sum(cell2mat(sum_bin_all(1,:)'))/sum(num_channels(1,:))*100;
data2 = sum(cell2mat(sum_bin_all(2,:)'))/sum(num_channels(2,:))*100;


if flagSaveData
    save([DataSaveFolder 'Speech'], 'data1', 'data2', 'sum_bin_all', 'num_channels', 'num_channels_total');
end 

