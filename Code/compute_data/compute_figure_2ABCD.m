%% computes data for motor imagery classification (figure 2)

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all 
close all

%% changeable parameters
flagGoData = false;
%if true  - computes classification for Go trials 
%if false - computes classification for NoGo trials 

flagShuffleTest = false; % 
%if true - shuffle labels to compute shuffled test
flagSaveData = true; % saves the data results 

SavedData = [pwd '\Data\NeuralData\MotorImagery\'];
unit_region = 'PMV'; %SMG, PMV or S1

%%
if ~flagShuffleTest
    number_repetitions = 1; 
else
    number_repetitions = 1000; 
end

%load data
cueType = ['Table_sorting_aligned_thr_-4.5_ImageCue_perTrials.mat'];
Data = load(fullfile(SavedData, cueType));

% use Go or NoGo data
if flagGoData
    Go_data = Data.Go_data;
else
    Go_data = Data.NoGo_data;
end 

Sessions = unique(Go_data.session_date);
number_sessions = length(Sessions);
labels = Go_data.GoLabels;
classifier = 'diaglinear'; 
number_phases = 4; 
number_channels = 96;
font_size = 12;

%predeclare variables 
k = 40;
errTrain = ones(number_repetitions,k,number_phases,number_sessions); 
errTest = ones(number_repetitions,k,number_phases,number_sessions);

keepTestingLabels = cell(number_repetitions,number_phases,number_sessions); 
keepPredictedLabels = cell(number_repetitions,number_phases,number_sessions);

session_idx_to_remove = [];
tic 
for n_session = 1:number_sessions
    disp(['Classifying session ' Sessions{n_session} ]);
    idx_this_session = ismember(Go_data.session_date, Sessions(n_session));
    
    if strcmp('SMG', unit_region)
        SessionData = Go_data.SMG_Go(idx_this_session,:);
    elseif strcmp('PMV', unit_region)
        SessionData = Go_data.PMV_Go(idx_this_session,:);
    elseif strcmp('S1', unit_region)
        SessionData = Go_data.S1X_Go(idx_this_session,:);
    end
     
    %remove session days that have no sorted units or less than 5 sorted units 
    if nnz(cell2mat(cellfun(@(x) isempty(x), SessionData, 'UniformOutput', false))) >1 || unique(cell2mat(cellfun(@(x) size(x,2), SessionData, 'UniformOutput', false))) < 5
        session_idx_to_remove = [session_idx_to_remove, n_session];
        disp(['Removing session ' Sessions{n_session}])
        continue
    end
   
    labels_per_session = labels(idx_this_session);
    time_phase_labels = Go_data.time_phase_labels(idx_this_session);
    %loop through number of repetitions
    for n_rep = 1:number_repetitions 
        %loop through number of task phases
        for n_phase = 1:number_phases
            LabelsPredRepRun = [];
            LabelsRepRun = [];
            %extract averaged firing rate over task phases
            data_per_phase = arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:)),SessionData,time_phase_labels, 'UniformOutput', false);       

            %randomize labels to generate shuffle distribution
            if flagShuffleTest
                labels_per_session = labels_per_session(randperm(length(labels_per_session)));
            end 

            labels_adapted = arrayfun(@(x,y) ones(size(x{1,1},1),1)*y ,data_per_phase,labels_per_session, 'UniformOutput', false);
            %shuffle partition
            cv = cvpartition(length(data_per_phase), 'LeaveOut');
            %loop through cross validation sets
            for cvRun = 1:cv.NumTestSets %
                %define training and testing data and labels
                trIdx = find(cv.training(cvRun));
                teIdx = find(cv.test(cvRun));
                data_training = cell2mat(data_per_phase(trIdx));
                LabelsTrain = cell2mat(labels_adapted(trIdx)); 
                data_testing = cell2mat(data_per_phase(teIdx)); 
                LabelsTest = cell2mat(labels_adapted(teIdx)); 

                %compute PCA
                [coeff, score, latent,~, explained] = pca(data_training); 
                variance_c = cumsum(explained);
                %find #PC's explaining 90% of explained variance
                idx_90 = find(variance_c>90);
                %project train and test data on coefficients
                PCA_data_train = data_training*coeff;
                PCA_data_test = data_testing*coeff;
                PCA_dimensions = idx_90(1); 

                %if there are less then 10 PC's that explain 90% of
                %the variance, we keep all features
                if PCA_dimensions < 10
                    if idx_90(end) >= 10
                        PCA_dimensions = 10;
                    else
                        PCA_dimensions = idx_90(end);
                    end
                end 

                data_training = PCA_data_train(:,1:PCA_dimensions);    
                data_testing = PCA_data_test(:,1:PCA_dimensions);

                %train model
                model = fitcdiscr(data_training, LabelsTrain, 'DiscrimType', classifier); 
                %model predictions
                predicted_training = predict(model, data_training);
                predicted_testing = predict(model, data_testing);
                %keep labels and predictions
                LabelsRepRun = [LabelsRepRun; LabelsTest];
                LabelsPredRepRun =  [LabelsPredRepRun; predicted_testing];
                %calculate classification error
                errTrain(n_rep,cvRun,n_phase,n_session) = analysis.classerror(LabelsTrain, predicted_training); 
                errTest(n_rep,cvRun,n_phase,n_session) = analysis.classerror(LabelsTest, predicted_testing); 
            end 
               keepTestingLabels{n_rep,n_phase, n_session} = LabelsRepRun;
               keepPredictedLabels{n_rep,n_phase, n_session} = LabelsPredRepRun;
               clear LabelsRepRun;
               clear LabelsPredRepRun;
        end 
    end 
end

%save data
shuffle_label = {'', 'Shuffle'};
GoNoGo = {'_NoGo',''};
if flagSaveData
    save([pwd '\Data\ClassificationMotorImagery\' unit_region '_Error'   shuffle_label{flagShuffleTest+1}  GoNoGo{flagGoData +1}], 'errTrain', 'errTest', 'keepTestingLabels', 'keepPredictedLabels','Sessions','number_sessions')
end 
