%% compute data for cross-task classiifcation (figure 4E)

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all 
close all

%% changeable parameters
flagShuffleTest = false; % 
%if true - shuffle labels to compute shuffled test
flagSaveData = false; % saves the data results 
SavedData = [pwd '\Data\NeuralData\SpeechDataset2\'];
%% variables
unit_region = 'SMG'; 
data_names = {'MotorImagery', 'Grasps', 'Colors'};
number_phases = 4;
number_tasks = length(data_names);
Sessions_tasks = {'20190911','20190923','20190924','20190930','20191016'}; 
data_per_phase = cell(number_tasks, number_phases, length(Sessions_tasks));
flagSelectTimeBins = false; 

%prepare data for cross-task classification
for n_task = 1:number_tasks 
    
    %load data
    data_name = data_names{n_task};
    Data = load(fullfile(SavedData,['Table_sorting_aligned_thr_-4.5_Speaking_' data_name 'PerModality.mat']));
    Go_data = Data.Go_data;
    
    Sessions = unique(Go_data.session_date);
    %Only keeping the sessions that have all three tasks 
    Sessions = Sessions(ismember(Sessions, Sessions_tasks));
    number_sessions = length(Sessions);

    %predeclare variables
    labels = Go_data.GoLabels;

    for n_session = 1:number_sessions       
        idx_this_session = ismember(Go_data.session_date, Sessions(n_session));

        if strcmp('SMG', unit_region) 
            SessionData = Go_data.SMG_Go(idx_this_session,:);
        elseif strcmp('PMV', unit_region) 
            SessionData = Go_data.PMV_Go(idx_this_session,:);
        elseif strcmp('S1X', unit_region)
            SessionData = Go_data.S1X_Go(idx_this_session,:);
        else
            error([unit_region ' does not exist '])
        end

        labels_per_session = labels(idx_this_session);
        time_phase_labels = Go_data.time_phase_labels(idx_this_session);
        
        %Select shorter cue phase window for color cue phase
        if n_task == 3
            flagSelectTimeBins = true;       
        end 
        
        if flagSelectTimeBins
            time_phase_labels_old = time_phase_labels{1};
            time_phase_labels_adapted = time_phase_labels_old;
            timeTrial = 0.05*(1:length(time_phase_labels_adapted));
            time_phase_labels_adapted(timeTrial > 4 & timeTrial < 6.2) = 0;
            tmp = cell(size(time_phase_labels));
            [tmp{:}] = deal(time_phase_labels_adapted);
            time_phase_labels = tmp;
        end
        
        for n_phase = 1:number_phases
            %extract data per phase
            data_per_phase_tmp = cell2mat(arrayfun(@(x,y) nanmean(x{1,1}(time_phase_labels{1}== n_phase,:),1),SessionData, 'UniformOutput', false));
            data_per_phase{n_task,n_phase,n_session} = cell2mat(arrayfun(@(x) data_per_phase_tmp(labels_per_session == x,:), unique(labels_per_session), 'UniformOutput', false)); 
        end               
    end
end

%%
%data per task
MotorImageryDataPerPhase = squeeze(data_per_phase(1,:,:))';
GraspSpeechDataPerPhase = squeeze(data_per_phase(2,:,:))';
ColorSpeechDataPerPhase = squeeze(data_per_phase(3,:,:))';
tasks = {MotorImageryDataPerPhase,GraspSpeechDataPerPhase,ColorSpeechDataPerPhase};

data_names = {'Motor Imagery', 'Spoken Grasps', 'Spoken Colors'};
%preallocate space 
errorsAll = cell(1, number_tasks);
DataNamesAll = cell(1, number_tasks);

%perform cross-task classification, using each time data from a different
%task to train the model. Test on all models

for n_task = 1:length(tasks)
    data_per_phase_all1 = tasks{n_task};

    flagPCA = true; 

    if isequal(data_per_phase_all1,MotorImageryDataPerPhase)

        data_name1 = data_names{1};
        data_per_phase_all2 = GraspSpeechDataPerPhase;    
        data_name2 = data_names{2};
        data_per_phase_all3 = ColorSpeechDataPerPhase;
        data_name3 = data_names{3};

    elseif isequal(data_per_phase_all1,GraspSpeechDataPerPhase) 

        data_name1 = data_names{2};
        data_per_phase_all2 = MotorImageryDataPerPhase;
        data_name2 = data_names{1};
        data_per_phase_all3 = ColorSpeechDataPerPhase;
        data_name3 = data_names{3};

    elseif isequal(data_per_phase_all1,ColorSpeechDataPerPhase)

        data_name1 =  data_names{3};
        data_per_phase_all2 = GraspSpeechDataPerPhase;
        data_name2 =  data_names{2};
        data_per_phase_all3 = MotorImageryDataPerPhase;
        data_name3 = data_names{1};    

    end 

    k = 40;
    classifier = 'diaglinear'; 
    number_repetition = 1; 
    %predeclare variables 
    err_test1 = zeros(number_repetition,k,number_phases); 
    err_test2 = zeros(number_repetition,k,number_phases);
    err_test3 = zeros(number_repetition,k,number_phases);

    labels = sort(labels_per_session);

    color_cues = util.get_color_rgb_codes({data_name1, data_name2, data_name3});
    %loop through sessions
    for n_session = 1:number_sessions 
        %loop through phases
        for n_phase = 1:number_phases
            %extract data for each phase
            data_per_phase1 = data_per_phase_all1{n_session,n_phase};
            data_per_phase2 = data_per_phase_all2{n_session,n_phase};
            data_per_phase3 = data_per_phase_all3{n_session,n_phase};
            cv = cvpartition(size(data_per_phase1,1), 'Leaveout');
            %loop through cross validation sets
            for cvRun = 1:cv.NumTestSets 
                %define training and testing data and labels
                trIdx = find(cv.training(cvRun));
                teIdx = find(cv.test(cvRun));

                data_training1 = data_per_phase1(trIdx,:);
                labels_train = labels(trIdx); 

                data_testing1 = data_per_phase1(teIdx,:); 
                labels_test = labels(teIdx);

                data_testing2 = data_per_phase2(teIdx,:); 
                data_testing3 = data_per_phase3(teIdx,:); 
                %perform PCA
                if flagPCA
                    [coeff, score, latent,~, explained] = pca(data_training1);
                    variance_c = cumsum(explained);
                    %find # PC explaining 90% of the variance
                    idx_90 = find(variance_c>90,1);
                    %project data onto PC
                    data_training1 = data_training1*coeff(:,1:idx_90);
                    data_testing1 = data_testing1*coeff(:,1:idx_90);
                    data_testing2 = data_testing2*coeff(:,1:idx_90);
                    data_testing3 = data_testing3*coeff(:,1:idx_90);

                end 
                %train classification model
                model = fitcdiscr(data_training1, labels_train, 'DiscrimType', classifier);
                %calculate classification error
                err_test1(n_session,cvRun,n_phase) = analysis.classerror(labels_test, predict(model, data_testing1)); 
                err_test2(n_session,cvRun,n_phase) = analysis.classerror(labels_test, predict(model, data_testing2));
                err_test3(n_session,cvRun,n_phase) = analysis.classerror(labels_test, predict(model, data_testing3)); 
            end
        end
    end
  
err_test_tmp{1} = err_test1;
err_test_tmp{2} = err_test2;
err_test_tmp{3} = err_test3;
labels_all = {data_name1, data_name2, data_name3};

%save data
errorsAll{n_task} = err_test_tmp;
DataNamesAll{n_task} = {labels_all};

end 

if flagSaveData
    if flagShuffleTest
        save([pwd '\Data\CrossModalityClassification\ShuffledData','DataNamesAll','errorsAll']);
    else
        save([pwd '\Data\CrossModalityClassification\RealData'], 'DataNamesAll','errorsAll');
    end 
end 

