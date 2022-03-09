%% compute and plot data for Venn diagram - figure S2

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all 
close all

%% parameters

saved_data = [pwd '\Data\NeuralData\SpeechDataset2'];
unit_region = 'SMG'; 
task_names = {'MotorImagery', 'Grasps', 'Colors'};
sessions = {'20190911','20190923','20190924','20190930','20191016',}; %Grasp Speaking All;
phase_names = {'ITI', 'Cue', 'Delay','Action'};

%predeclare variables 
number_tasks = length(task_names);
number_sessions = length(sessions);
number_phases = length(phase_names);
tuned_channels_per_phase_all =cell(number_tasks,number_sessions);
data_per_phase = cell(number_tasks,number_phases, number_sessions);

for n_task = 1:number_tasks 
   
    data_name = task_names{n_task};
    %load dataset
    dataset_name = fullfile(saved_data,['Table_sorting_aligned_thr_-4.5_Speaking_' data_name 'PerModality.mat']);
    disp(['Load ' dataset_name]);
    data = load(dataset_name);  
    Go_data = data.Go_data;

    %predeclare variables
    labels = Go_data.GoLabels;
    
    %loop through sessions and calculate unit tuning
    for n_session = 1:number_sessions
        if number_sessions ~=1
            disp(['Tuning session ' sessions{n_session} ]);
        end
        idx_current_session = ismember(Go_data.session_date, sessions(n_session));       

        if strcmp('SMG', unit_region)
            SessionData = Go_data.SMG_Go(idx_current_session,:);
        elseif strcmp('PMV', unit_region) 
            SessionData = Go_data.PMV_Go(idx_current_session,:);
        elseif strcmp('S1X', unit_region)
            SessionData = Go_data.S1X_Go(idx_current_session,:);
        else
            error([unit_region ' does not exist '])
        end

        labels_per_session = labels(idx_current_session);
        time_phase_labels = Go_data.time_phase_labels(idx_current_session);
           
        %Select which time bins to keep in the analysis 
        flagSelectTimeBins = true; 
        
        if flagSelectTimeBins 
            time_phase_labels_adapted = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;...
                            2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;...
                            0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;...
                            0;0;0;0;0;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;0;0;0;0;0;0;0;0;0;0;0;0;...
                            0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
                        
            tmp = cell(size(time_phase_labels));
            [tmp{:}] = deal(time_phase_labels_adapted);
            time_phase_labels = tmp;
        end
        
        %computed tuned channels per phase 
        [~, tuned_units_per_phase] = analysis.get_tuned_units(SessionData,labels_per_session,time_phase_labels,false);                    
          tuned_channels_per_phase_all{n_task,n_session} = tuned_units_per_phase;      
    end
end


%% plot figures S2 (Venn diagrams)

%calculate percentage of overlapping tuning units
phases_to_evaluate_all = [2, 4];

for n_phase = 1:length(phases_to_evaluate_all)
    %compute for each phase separately
    phase_to_evaluate = phases_to_evaluate_all(n_phase);
    phase_tuning = cellfun(@(x) x(:,phase_to_evaluate),tuned_channels_per_phase_all, 'UniformOutput', false);
    %predeclare variables 
    common_ = zeros(1,number_sessions);
    MI_ = zeros(1,number_sessions);
    GS_ = zeros(1,number_sessions);
    CS_ = zeros(1,number_sessions);
    I_MI_GS = zeros(1,number_sessions);
    I_MI_CS = zeros(1,number_sessions);
    I_GS_CS = zeros(1,number_sessions);
    sum_all = zeros(number_sessions,3);
    %loop through different sessions
    for n_session = 1:number_sessions

        phase_tmp = cell2mat(phase_tuning(:,n_session)');

        T_MI = nnz(phase_tmp(:,1));
        T_GS = nnz(phase_tmp(:,2));
        T_CS = nnz(phase_tmp(:,3));
        %units tuned in all three tasks
        Common = find(sum(phase_tmp,2) == 3);
        phase_tmp(Common,:) = [];
        common_(n_session) = nnz(Common);

        MI = find(phase_tmp(:,1));
        GS = find(phase_tmp(:,2));
        CS = find(phase_tmp(:,3));

        %calculate intersecting tuned units
        MI_GS = intersect(MI,GS);
        MI_CS = intersect(MI,CS);
        GS_CS = intersect(GS,CS);

        %calculate units tuned only to one task
        MIt = nnz(MI) - length(MI_GS) - length(MI_CS);
        GSt = nnz(GS) - length(MI_GS) - length(GS_CS);
        CSt = nnz(CS) - length(MI_CS) - length(GS_CS);

        MI_(n_session)= MIt;
        GS_(n_session) = GSt;
        CS_(n_session) = CSt;

        I_MI_GS(n_session) = length(MI_GS);
        I_MI_CS(n_session) = length(MI_CS);
        I_GS_CS(n_session) = length(GS_CS);

        sum_all(n_session,1) = MIt + length(MI_GS) + length(MI_CS) + length(Common);
        sum_all(n_session,2) = GSt + length(MI_GS) + length(GS_CS) + length(Common);
        sum_all(n_session,3) = CSt + length(GS_CS) + length(MI_CS) + length(Common);

    end 
    %plot venn diagram 
    figure()
    util.vennX([sum(MI_), sum(I_MI_GS), sum(GS_), sum(I_GS_CS),sum(CS_), sum(I_MI_CS), sum(common_)], 0.01);

    phase_names = {'ITI','Cue', 'Delay', 'Action'};
    title(phase_names{phase_to_evaluate});
    
end
